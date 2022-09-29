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
from mia_processes.utils import recupCover, PageNumCanvas, ReportLine, slice_planes_plot

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
from reportlab.lib.units import mm, inch
from reportlab.lib import colors
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                 Table, TableStyle, PageBreak)
import tempfile
import platform
import json
#import numpy as np
#from nibabel.processing import resample_from_to
#import matplotlib.pyplot as plt
#import scipy.ndimage as ndi



class MRIQC_report(ProcessMIA):
    """
        * Generates the report for MRIQC pipeline *

        Please, see the complete documentation for the `MRIQC_report brick in the populse.mia_processes website
        https://populse.github.io/mia_processes/documentation/bricks/reports/MRIQC_report.html

        """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(MRIQC_report, self).__init__()

        # Third party software required for the execution of the brick
        self.requirement = []

        # Inputs description

        IQMs_file_desc = 'A .JSON file containing the IQMs'

        anat_desc = ('An existing, uncompressed anatomic image file (valid '
                       'extensions: .nii)')

        fig_rows_desc = 'The number of lines for the slice planes plot'

        fig_cols_desc = 'The number of columns for the slice planes plot'

        inf_slice_start_desc = 'The first index displayed in slice planes plot'

        slices_gap_desc = 'Gap between slices in slice planes plot'


        # Outputs description
        report_desc = 'The generated report (pdf)'

        # Inputs traits
        self.add_trait("IQMs_file",
                       File(copyfile=False,
                            output=False,
                            optional=False,
                            desc=IQMs_file_desc))

        self.add_trait("anat",
                       ImageFileSPM(copyfile=False,
                                    output=False,
                                    optional=False,
                                    desc=anat_desc))


        self.add_trait("fig_rows", traits.Int(5, output=False, optional=False, desc=fig_rows_desc))
        self.add_trait("fig_cols", traits.Int(5, output=False, optional=False, desc=fig_cols_desc))
        self.add_trait("inf_slice_start", traits.Either(Undefined, traits.Int, output=False, optional=True, desc=inf_slice_start_desc))
        self.add_trait("slices_gap", traits.Either(Undefined, traits.Int, output=False, optional=True, desc=slices_gap_desc))

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

        header_image_1 = os.path.join(sources_images.__path__[0],
                                    'Logo_populse_square.jpg')
        header_image_1 = Image(header_image_1, 43.0 * mm, 43.0 * mm)  # 2156px x 2156px

        header_image_2 = os.path.join(sources_images.__path__[0],
                                    'Logo_populse_mia_HR.jpeg')
        header_image_2 = Image(header_image_2, 54.4 * mm, 40.0 * mm)  # 2505px x 1843px

        #header_title = '<font size = 20><sup>$</sup></font><font size=30><b>MRI&nbsp Q</b></font>  <font size=11>uality</font>   <font size=30><b>&nbsp C</b></font>   <font size=11>ontrol</font>'
        header_title = '<font size = 20><sup>$</sup></font><font size=30><b>MRIQ</b></font><font size=11>uality</font><font size=30><b>C</b></font><font size=11>ontrol</font>'

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
            temp_file.write(bytes("SITE:{}".format(site), encoding='utf8'))
            temp_file.write(bytes("\nMRI SCANNER: {}".format(self.dict4runtime['mri_scanner']), encoding='utf8'))
            temp_file.write(bytes("\nSTUDY NAME: {}".format(study_name), encoding='utf8'))
            temp_file.write(bytes("\nEXAMINATION DATE: {}".format(acqu_date), encoding='utf8'))
            temp_file.write(bytes("\nMRIQC CALCULATION DATE: {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')), encoding='utf8'))
            temp_file.write(bytes("\nNAME OF THE INPUT DATA: {}".format(self.anat), encoding='utf8'))
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
        image_cov = Table([[header_image_1, header_image_2]], [70*mm, 70*mm])
        image_cov.setStyle(TableStyle([
                                  ('ALIGN', (0, 0), (0, 0), 'LEFT'),
                                  ('ALIGN', (-1, 0), (-1, 0), 'RIGHT'),
                                  ('VALIGN', (0, 0), (-1, 0), 'MIDDLE')
                                      ]))

        report.append(image_cov)
        report.append(Spacer(0 * mm, 8 * mm))  # (width, height)

        report.append(Paragraph(header_title,
                                 styles['Center']))

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        #line = ReportLine(250)
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 10 * mm))

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

        report.append(Spacer(0 * mm, 5 * mm))

        report.append(Paragraph(textDisclaimer,
                                 styles["Justify"]))

        report.append(Spacer(0 * mm, 5 * mm))  # (width, height)

        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        report.append(Paragraph(
            "<font size = 8><sup>$</sup>Esteban O et al., <i>MRIQC: Advancing the Automatic Prediction of Image Quality in MRI from Unseen Sites</i>, PLOS ONE 12(9):e0184661.</font>",
            styles['Left']))

        report.append(PageBreak())

        # Second page - IQMs

        f = open(self.IQMs_file)
        data = json.load(f)
        f.close()

        report.append(Paragraph("<font size = 18 > <b> Image parameters <br/> </b> </font>",
                                     styles['Center']))

        report.append(Spacer(0 * mm, 4 * mm))  # (width, height)

        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 20 * mm))

        # Spatial resolution
        report.append(Paragraph(
            "<font size = 15 > <b>SPATIAL RESOLUTION</b> </font> length - spacing:",
            styles['Bullet1']))

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        # size_x
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Voxel size in X </b> </font>(size_x): ' + \
                str(round(data.get('size_x'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Voxel size in X </b> </font>(size_x): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # size_y
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Voxel size in Y </b> </font>(size_y): ' + \
                str(round(data.get('size_y'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Voxel size in Y </b> </font>(size_y): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # size_z
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Voxel size in Z </b> </font>(size_z): ' + \
                str(round(data.get('size_z'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Voxel size in Z </b> </font>(size_z): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # spacing_x
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Spacing in X (mm)</b> </font>(spacing_x): ' + \
                str(round(data.get('spacing_x'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Spacing in X (mm)</b> </font>(spacing_x): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # spacing_y
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Spacing in Y (mm)</b> </font>(spacing_y): ' + \
                str(round(data.get('spacing_y'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Spacing in Y (mm)</b> </font>(spacing_y): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # spacing_z
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Spacing in Z (mm)</b> </font>(spacing_z): ' + \
                str(round(data.get('spacing_z'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Spacing in Z (mm)</b> </font>(spacing_z): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 20 * mm))

        # Signal-to-noise
        report.append(Paragraph(
            "<font size = 15 > <b>NOISE</b> </font> impact of noise and/or evaluation of the fitness of a noise model:",
            styles['Bullet1']))

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        # snr_csf
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Signal-to-noise ratio (SNR) for cerebrospinal fluid </b> </font>(snr_csf): ' + \
                str(round(data.get('snr_csf'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Signal-to-noise ratio (SNR) for cerebrospinal fluid </b> </font>(snr_csf): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # snr_wm
        try:
            report.append(Paragraph(
                '<font size = 11> <b> SNR for white matter </b> </font>(snr_wm): ' + \
                str(round(data.get('snr_wm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> SNR for white matter </b> </font>(snr_wm): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # snr_gm
        try:
            report.append(Paragraph(
                '<font size = 11> <b> SNR for gray matter </b> </font>(snr_gm): ' + \
                str(round(data.get('snr_gm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> SNR for gray matter </b> </font>(snr_gm): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # snr_total
        try:
            report.append(Paragraph(
                '<font size = 11> <b> SNR for brain parenchyma </b> </font>(snr_total): ' + \
                str(round(data.get('snr_total'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> SNR for brain parenchyma </b> </font>(snr_total): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # snrd_csf
        try:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for cerebrospinal fluid</b></font> (snrd_csf): ' + \
                str(round(data.get('snrd_csf'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for cerebrospinal fluid</b></font> (snrd_csf): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # snrd_wm
        try:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for white matter</b></font> (snrd_wm): ' + \
                str(round(data.get('snrd_wm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for white matter</b></font> (snrd_wm): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # snrd_gm
        try:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for gray matter</b></font> (snrd_gm): ' + \
                str(round(data.get('snrd_gm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for gray matter</b></font> (snrd_gm): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # snrd_total
        try:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for brain parenchyma</b></font> (snrd_total): ' + \
                str(round(data.get('snrd_total'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for brain parenchyma</b></font> (snrd_total): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # cnr
        try:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>Contrast-to-noise ratio</b></font> (cnr): ' + \
                str(round(data.get('cnr'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>Contrast-to-noise ratio</b></font> (cnr): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # qi_2
        try:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11><b>Mortamet’s quality index 2</b></font> (qi_2): ' + \
                '{:.2e}'.format(data['qi_2']),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11><b>Mortamet’s quality index 2</b></font> (qi_2): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # cjv
        try:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Coefficient of joint variation</b></font> (cjv): ' + \
                str(round(data.get('cjv'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Coefficient of joint variation</b></font> (cjv): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 65 * mm))  # (width, height)

        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        report.append(Paragraph(
            "<font size = 8><sup>$</sup>Dietrich et al., <i>Measurement of SNRs in MR images: influence of multichannel coils, parallel imaging and reconstruction filters</i>, JMRI 26(2):375–385, 2007. Higher values are better.</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>#</sup>Magnotta, VA., & Friedman, L., <i>Measurement of signal-to-noise and contrast-to-noise in the fBIRN multicenter imaging study</i>, J Dig Imag 19(2):140-147, 2006. Higher values are better.</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>&</sup>Mortamet B et al., <i>Automatic quality assessment in structural brain magnetic resonance imaging</i>, Mag Res Med 62(2):365-372, 2009. Lower values are better.</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>%</sup>Ganzetti et al., <i>Intensity inhomogeneity correction of structural MR images: a data-driven approach to define input algorithm parameters</i>, Front Neuroinform 10:10, 2016. Lower values are better.</font>",
            styles['Left']))

        report.append(PageBreak())

        # third page - IQMs

        report.append(Paragraph("<font size = 18 > <b> Image parameters <br/> </b> </font>",
                                     styles['Center']))

        report.append(Spacer(0 * mm, 4 * mm))  # (width, height)

        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 20 * mm))

        # Spatial distribution
        report.append(Paragraph(
            "<font size = 15 > <b>SPATIAL DISTRIBUTION</b> </font> information theory to evaluate the spatial distribution of information:",
            styles['Bullet1']))

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        # fber
        try:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>FBER</b></font> (fber): ' + \
                str(round(data.get('fber'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>FBER</b></font> (fber): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # efc
        try:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>EFC</b></font> (efc): ' + \
                str(round(data.get('efc'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>EFC</b></font> (efc): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 20 * mm))

        # Artifacts
        report.append(Paragraph(
            "<font size = 15 > <b>ARTIFACTS</b> </font> look for the presence and impact of particular artifacts (specifically, the INU artifact), and the signal leakage due to rapid motion (e.g. eyes motion or blood vessel pulsation):",
            styles['Bullet1']))

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        # wm2max
        try:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11> <b> White-matter to maximum intensity ratio </b> </font>(wm2max): ' + \
                str(round(data.get('wm2max'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11> <b> White-matter to maximum intensity ratio </b> </font>(wm2max): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # qi_1
        try:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11> <b>Mortamet’s quality index 1</b> </font>(qi_1): ' + \
                str(round(data.get('qi_1'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11> <b>Mortamet’s quality index 1</b> </font>(qi_1): ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)


        report.append(Spacer(0 * mm, 65 * mm))  # (width, height)

        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        report.append(Paragraph(
            "<font size = 8><sup>%</sup>Shehzad Z et al., <i>The Preprocessed Connectomes Project Quality Assessment Protocol - a resource for measuring the quality of MRI data</i>, Front. Neurosci. Conference Abstract: Neuroinformatics 2015. Higher values are better.</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>#</sup>Atkinson et al., <i>Automatic correction of motion artifacts in magnetic resonance images using an entropy focus criterion</i>, IEEE Trans Med Imag 16(6):903-910, 1997. Lower values are better.</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>&</sup>The median intensity within the white matter mask over the 95% percentile of the full intensity distribution. Values should be around the interval [0.6, 0.8].</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>$</sup>Mortamet B et al., <i>Automatic quality assessment in structural brain magnetic resonance imaging</i>, Mag Res Med 62(2):365-372, 2009. Lower values are better.</font>",
            styles['Left']))

        report.append(PageBreak())

        #fourth page

        report.append(Paragraph("<font size = 18 > <b> MRI axial slice planes display <br/> </b> </font>",
                                     styles['Center']))

        report.append(Spacer(0*mm, 4*mm)) # (width, height) ## In Linux

        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0*mm, 15*mm)) # (width, height)
        report.append(Paragraph('<font size = 9 > <i> "Neurological" convention, the left side of the image corresponds to the left side of the brain. </i> <br/> </font>',
                                     styles['Center']))
        report.append(Spacer(0*mm, 1*mm)) # (width, height)

        # fig_rows = 5
        # fig_cols = 5
        # n_subplots = fig_rows * fig_cols
        # brain_img = nib.load(self.anat)
        # brain_voxels = brain_img.get_fdata()
        #
        # # axial
        # n_slice = brain_voxels.shape[1]
        # step_size = n_slice // n_subplots
        # plot_range = n_subplots * step_size
        # start_stop = int((n_slice - plot_range) / 2)

        tmpdir = tempfile.TemporaryDirectory()

        # fig, axs = plt.subplots(fig_rows, fig_cols, figsize=[10, 10])
        #
        # for idx, img in enumerate(range(start_stop, plot_range, step_size)):
        #     #axs.flat[idx].imshow(ndi.rotate(t1_voxels[:, img, :], 270),  cmap='gray')
        #     axs.flat[idx].imshow(brain_voxels[:, img, :],  cmap='gray')
        #     axs.flat[idx].axis('off')
        #
        # plt.tight_layout()
        #plt.show()
        #plt.savefig(os.path.join(tmpdir.name, 'slices.png'), bbox_inches='tight')

        #slices_image = os.path.join(sources_images.__path__[0],
        #                            'Logo_populse_mia_HR.jpeg')

        #slices_image = slice_planes_plot(self.anat, fig_rows=5, fig_cols=5, inf_slice_start=None, slices_gap=None, cmap="Greys_r", out_dir=tmpdir.name)
        slices_image = slice_planes_plot(self.anat, self.fig_rows, self.fig_cols, inf_slice_start=self.inf_slice_start, slices_gap=self.slices_gap, cmap="Greys_r", out_dir=tmpdir.name)
        # reminder: A4 == 210mmx297mm
        #slices_image = Image(slices_image, width=177.4 * mm, height=222.0 * mm)  #791x990
        #slices_image = Image(slices_image, width=190 * mm, height=230 * mm)
        slices_image = Image(slices_image, width=7.4803 * inch, height=9.0551 * inch)
        slices_image.hAlign = 'CENTER'
        report.append(slices_image)


        page.build(report, canvasmaker=PageNumCanvas)
        tmpdir.cleanup()
