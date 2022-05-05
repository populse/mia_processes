"""The reporting library of the mia_processes package.

The purpose of this module is to provide the reporting bricks necessary to
generate automatic report at the end of a pipeline calculation.

:Contains:
    :Class:
        - MRIQC_report

    :Function:
        -
"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# nibabel import
import nibabel as nib

# nipype import
from nipype.interfaces.base import (OutputMultiPath, InputMultiPath, File,
                                    traits, TraitListObject, Undefined,
                                    DictStrStr, Str)
from nipype.interfaces.spm.base import ImageFileSPM

#mia_processes import:
from mia_processes.utils import recupCover, PageNumCanvas, ReportLine

# populse_mia import
from populse_mia import sources_images
from populse_mia.data_manager.project import COLLECTION_CURRENT
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA



# soma-base imports
from soma.qt_gui.qt_backend.Qt import QMessageBox


# Other import
from datetime import datetime
import os
from sys import exit, path, version
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.lib.pagesizes import A4, landscape, portrait
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                 Table, TableStyle, PageBreak)
import tempfile
import platform
import nibabel as nib
from nibabel.processing import resample_from_to
import matplotlib.pyplot as plt
import scipy.ndimage as ndi



class MRIQC_report(ProcessMIA):
    """
        * Generates the report for MRIQC pipeline *

        Please, see the complete documentation for the `MRIQC_report brick in the populse.mia_processes website
        https://populse.github.io/mia_processes/documentation/bricks/reports/MRIQC_report.html

        """

    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(MRIQC_report, self).__init__()

        # Third party software required for the execution of the brick
        self.requirement = []

        # Inputs description
        anat_desc = ('An existing, uncompressed anatomic image file (valid '
                       'extensions: .nii)')

        # Outputs description
        report_desc = 'The generated report (pdf)'

        # Inputs traits
        self.add_trait("anat",
                       ImageFileSPM(copyfile=False,
                                    output=False,
                                    optional=False,
                                    desc=anat_desc))
        # Outputs traits
        self.add_trait("report",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=report_desc))

        self.add_trait("dict4runtime",
                       traits.Dict(output=False,
                                   optional=False,
                                   userlevel=1))
        self.dict4runtime['patient_name'] = Undefined
        self.dict4runtime['study_name'] = Undefined
        self.dict4runtime['acquisition_date'] = Undefined
        self.dict4runtime['sex'] = Undefined
        self.dict4runtime['site'] = Undefined
        self.dict4runtime['mri_scanner'] = Undefined
        self.dict4runtime['age'] = Undefined


        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(MRIQC_report, self).list_outputs()

        # As we do not have access to the database at the runtime (see #272),
        # we prepare here the data that the run_process_mia method will need
        # via dict4runtime dictionary
        # patient_name
        file_position = (self.anat.find(self.project.getName())
                         + len(self.project.getName()) + 1)
        database_filename = self.anat[file_position:]


        if ('PatientName' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['patient_name'] = self.project.session.get_value(COLLECTION_CURRENT,
                                                  database_filename,
                                                  ('PatientName'))

        if self.dict4runtime['patient_name'] in ('', Undefined):
            self.dict4runtime['patient_name'] = "Undefined_name_ref"

        if self.anat and self.anat not in ["<undefined>", traits.Undefined]:
            self.outputs['report'] = os.path.join(self.output_directory, self.dict4runtime['patient_name'] + '_mriqcReport_' + datetime.now().strftime('%Y_%m_%d_%H_%M_%S_%f')[:22] +'.pdf')



        # study_name
        if ('StudyName' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['study_name'] = self.project.session.get_value(COLLECTION_CURRENT,
                                                  database_filename,
                                                  ('StudyName'))
        # acquisition_date
        if ('AcquisitionDate' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['acquisition_date'] = str(self.project.session.get_value(COLLECTION_CURRENT,
                                                  database_filename,
                                                  ('AcquisitionDate')))

        # sex
        if ('Sex' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['sex'] = self.project.session.get_value(COLLECTION_CURRENT,
                                                  database_filename,
                                                  ('Sex'))
        # site
        if ('Site' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['site'] = self.project.session.get_value(COLLECTION_CURRENT,
                                                  database_filename,
                                                  ('Site'))
        # mri_scanner
        if ('Spectro' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['mri_scanner'] = self.project.session.get_value(COLLECTION_CURRENT,
                                                  database_filename,
                                                  ('Spectro'))
        # age
        if ('Age' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['age'] = self.project.session.get_value(COLLECTION_CURRENT,
                                                  database_filename,
                                                  ('Age'))

        # FIXME: Currently the following data is hard-coded. A solution should
        #        be found to retrieve them automatically or to put them in the
        #        input parameters of the brick:
        # Site
        if self.dict4runtime['site'] in ('', Undefined):
            self.dict4runtime['site'] = 'Grenoble University Hospital - CLUNI'

        # MriScanner
        if self.dict4runtime['mri_scanner'] in ('', Undefined):
            self.dict4runtime['mri_scanner'] = 'Philips Achieva 3.0T TX'


        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(MRIQC_report, self).run_process_mia()

        header_image = os.path.join(sources_images.__path__[0],
                                    'Logo_populse_mia_HR.jpeg')
        header_image = Image(header_image, 40.8 * mm, 30.0 * mm)  # 2505px x 1843px

        header_title = '<font size=30><b>MRI&nbsp Q</b></font>  <font size=11>uality</font>   <font size=30><b>&nbsp C</b></font>   <font size=11>ontrol</font>'

        if 'site' in self.dict4runtime and self.dict4runtime['site'] not in ('', Undefined):
            site = self.dict4runtime['site']

        else:
            site = "*Undefined site*"

        if 'study_name' in self.dict4runtime and self.dict4runtime['study_name'] not in ('', Undefined):
            study_name = self.dict4runtime['study_name']

        else:
            study_name = "*Undefined study name*"

        if 'acquisition_date' in self.dict4runtime and self.dict4runtime['acquisition_date'] not in ('', Undefined):
            acqu_date = self.dict4runtime['acquisition_date']

        else:
            acqu_date = "*Undefined study name*"

        if 'patient_name' in self.dict4runtime and self.dict4runtime['patient_name'] not in ('', Undefined):
            patient_ref = self.dict4runtime['patient_name']

        else:
            patient_ref = "*Undefined patient reference*"

        if 'sex' in self.dict4runtime and self.dict4runtime['sex'] not in ('', Undefined):
            patient_sex = self.dict4runtime['sex']

        else:
            patient_sex = "*Undefined patient sex*"

        if 'age' in self.dict4runtime and self.dict4runtime['age'] not in ('', Undefined):
            patient_age = self.dict4runtime['age']

        else:
            patient_age = "*Undefined patient age*"

        with tempfile.NamedTemporaryFile() as temp_file:
            temp_file.write(bytes("SITE:{}".format(self.dict4runtime['site']), encoding='utf8'))
            temp_file.write(bytes("\nMRI SCANNER: {}".format(self.dict4runtime['mri_scanner']), encoding='utf8'))
            temp_file.write(bytes("\nSTUDY NAME: {}".format(study_name), encoding='utf8'))
            temp_file.write(bytes("\nEXAMINATION DATE: {}".format(acqu_date), encoding='utf8'))
            temp_file.write(bytes("\nMRIQC CALCULATION DATE: {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')), encoding='utf8'))
            temp_file.write(bytes("\nPATIENT REFERENCE: {}".format(patient_ref), encoding='utf8'))
            temp_file.write(bytes("\nPATIENT SEX: {}".format(patient_sex), encoding='utf8'))
            temp_file.write(bytes("\nPATIENT AGE: {}".format(patient_age), encoding='utf8'))
            temp_file.write(bytes("\nSOFTWARES:Python {}".format(version.split()[0]), encoding='utf8'))
            from populse_mia import info as mia_info
            temp_file.write(bytes("\n : populse_mia {}".format(mia_info.__version__), encoding='utf8'))
            from capsul import info as capsul_info
            temp_file.write(bytes("\n : capsul {}".format(capsul_info.__version__), encoding='utf8'))
            from nipype import info as nipype_info
            temp_file.write(bytes("\n : nipype {}".format(nipype_info.__version__), encoding='utf8'))
            from mia_processes import info as mia_processes_info
            temp_file.write(bytes("\n : mia_processes {}".format(mia_processes_info.__version__), encoding='utf8'))

            if platform.system() == 'Darwin':
                os_sys = 'Mac OS X'

            else:
                os_sys = platform.system()

            temp_file.write(bytes("\n : Operating System {0} {1}".format(os_sys, platform.release()), encoding='utf8'))
            temp_file.flush()

            cover_data = recupCover(temp_file.name)

        cover_data = Table(cover_data, [60*mm, 97*mm])  # colWidths, rowHeights
        output = self.report

        page = SimpleDocTemplate(output,
                                 pagesize=portrait(A4),
                                 rightMargin=20*mm,
                                 leftMargin=20*mm,
                                 topMargin=10*mm,
                                 bottomMargin=10*mm)                          # Output document template definition for page; margins and pages size.

        styles = getSampleStyleSheet()                                          # Initialises stylesheet with few basic heading and text styles, return a stylesheet object.
        styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))      # ParagraphStyle gives all the attributes available for formatting paragraphs.
        styles.add(ParagraphStyle(name='Center2', alignment=TA_CENTER))
        styles['Center2'].leading = 24  # If more than 1 line.
        styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
        styles.add(ParagraphStyle(name='Right', alignment=TA_RIGHT))
        styles.add(ParagraphStyle(name='Left', alignment=TA_LEFT))
        styles.add(ParagraphStyle(name='Bullet1',
                                  leftIndent=30,
                                  bulletOffsetY=2,
                                  bulletIndent=20,
                                  bulletFontSize=6,
                                  bulletColor='black',
                                  bulletText=u'●')
                                 )
        styles.add(ParagraphStyle(name='Bullet2',
                                  leftIndent=60,
                                  bulletOffsetY=1,
                                  bulletIndent=50,
                                  bulletFontSize=6,
                                  bulletColor='black',
                                  bulletText=u'❍')
                                 )

        title = ('<font size=18> <b>Image-Quality Metrics summary report</b><br/> <br/></font>')


        textDisclaimer = ("<font size=7 ><i><b>DISCLAIMER</b><br/>Mia software,"
                          " from the Populse project, is executed in a "
                          "research environment on anonymized data. The "
                          "conclusions obtained with Mia software are an help "
                          "to diagnose and prognosticate. They do not "
                          "substitute themselves to the clinical care of the "
                          "physicians and remain under their responsibilities. "
                          "Consequently, Populse team is not responsible for "
                          "any direct or indirect damages resulting from the "
                          "use of data, informations, or results stemming from "
                          "the Mia software. The user recognizes to use these "
                          "informations under his sole and exclusive "
                          "responsibility.<br/> <br/> <b>DECHARGE DE "
                          "RESPONSABILITE</b><br/>Le logiciel Mia, provenant du"
                          " projet Populse, est exécuté dans un environnement "
                          "de recherche sur des données anonymisées. Les "
                          "conclusions obtenues grâce au logiciel Mia sont une "
                          "aide au diagnostic et au pronostic. Elles ne se "
                          "substituent pas à la prise en charge médicale des "
                          "médecins et demeurent sous leurs responsabilités. "
                          "Par conséquent, l'équipe Populse ne peut être tenu "
                          "responsable de dommage direct ou indirect résultant "
                          "de l'utilisation des données, des informations ou "
                          "des résultats issus du logiciel Mia. L'utilisateur "
                          "reconnaît utiliser ces informations sous sa seule et"
                          " entière responsabilité.</i> </font>")

        #Canvas creation
        report = []

        # First page - cover

        header_image.hAlign = 'CENTER'
        report.append(header_image)

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        report.append(Paragraph(header_title,
                                 styles['Center']))

        report.append(Spacer(0 * mm, 12 * mm))  # (width, height)  ## In Linux
        # report.append(Spacer(0*mm, 4*mm)) # (width, height)  ## In Mac

        line = ReportLine(250)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 20 * mm))

        report.append(Paragraph(title,
                                 styles['Center']))

        report.append(Spacer(0 * mm, 10 * mm))

        cover_data.setStyle(TableStyle([
            ('LINEABOVE', (0, 0), (-1, 0), 2, colors.black),
            ('LINEBELOW', (0, -1), (-1, -1), 2, colors.black),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE')
        ]))  # Carefull: Table use (col,raw) and 0-base for top left start as usual OR -1-base for lower right start.

        cover_data.hAlign = 'CENTER'
        report.append(cover_data)

        report.append(Spacer(0 * mm, 10 * mm))

        report.append(Paragraph(textDisclaimer,
                                 styles["Justify"]))

        report.append(PageBreak())

        report.append(Paragraph("<font size = 18 > <b> MRI axial slice planes display <br/> </b> </font>",
                                     styles['Center']))

        report.append(Spacer(0*mm, 4*mm)) # (width, height) ## In Linux

        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0*mm, 25*mm)) # (width, height)

        fig_rows = 5
        fig_cols = 5
        n_subplots = fig_rows * fig_cols
        brain_img = nib.load(self.anat)
        brain_voxels = brain_img.get_fdata()

        # axial
        n_slice = brain_voxels.shape[1]
        step_size = n_slice // n_subplots
        plot_range = n_subplots * step_size
        start_stop = int((n_slice - plot_range) / 2)

        tmpdir = tempfile.TemporaryDirectory()

        fig, axs = plt.subplots(fig_rows, fig_cols, figsize=[10, 10])

        for idx, img in enumerate(range(start_stop, plot_range, step_size)):
            #axs.flat[idx].imshow(ndi.rotate(t1_voxels[:, img, :], 270),  cmap='gray')
            axs.flat[idx].imshow(brain_voxels[:, img, :],  cmap='gray')
            axs.flat[idx].axis('off')

        plt.tight_layout()
        #plt.show()
        plt.savefig(os.path.join(tmpdir.name, 'slices.png'), bbox_inches='tight')

        #slices_image = os.path.join(sources_images.__path__[0],
        #                            'Logo_populse_mia_HR.jpeg')
        slices_image = Image(os.path.join(tmpdir.name, 'slices.png'), 177.4 * mm, 222.0 * mm)  #791x990
        slices_image.hAlign = 'CENTER'
        report.append(slices_image)


        page.build(report, canvasmaker=PageNumCanvas)
        tmpdir.cleanup()