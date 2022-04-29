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
from mia_processes.utils import recupCover

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



        if self.anat and self.anat not in ["<undefined>", traits.Undefined]:
            self.outputs['report'] = os.path.join(self.output_directory, 'report.pdf')


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
        # study_name
        if ('StudyName' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['study_name'] = self.project.session.get_value(COLLECTION_CURRENT,
                                                  database_filename,
                                                  ('StudyName'))
        # acquisition_date
        if ('AcquisitionDate' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['acquisition_date'] = self.project.session.get_value(COLLECTION_CURRENT,
                                                  database_filename,
                                                  ('AcquisitionDate'))
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
                                    'Logo_populse_mia_LR.jpeg')
        header_image = Image(header_image, 175.0 * mm,
                            30.0 * mm)  # 2508px × 430px

        header_title = '<font size=30><b>MRI &nbsp Q</b></font><font size=11>uality</font>\
                          <font size=30><b>&nbsp C</b></font><font size=11><ontrol/font>'


        if self.dict4runtime['site'] in ('', Undefined):
            site = "*Undefined site*"

        else:
            site = self.dict4runtime['site']

        if self.dict4runtime['study_name'] in ('', Undefined):
            study_name = "*Undefined study name*"

        else:
            study_name = self.dict4runtime['study_name']

        if self.dict4runtime['acquisition_date'] in ('', Undefined):
            acqu_date = "*Undefined study name*"

        else:
            acqu_date = self.dict4runtime['acquisition_date']

        if self.dict4runtime['patient_name'] in ('', Undefined):
            patient_ref = "*Undefined patient reference*"

        else:
            patient_ref = self.dict4runtime['patient_name']

        if self.dict4runtime['sex'] in ('', Undefined):
            patient_sex = "*Undefined patient sex*"

        else:
            patient_sex = self.dict4runtime['sex']

        if self.dict4runtime['age'] in ('', Undefined):
            patient_age = "*Undefined patient age*"

        else:
            patient_age = self.dict4runtime['age']


        run_date = datetime.now()

        with tempfile.TemporaryFile() as temp_file:
            temp_file.write("SITE:%s" % site)
            temp_file.write("\nMRI SCANNER:%s" % self.dict4runtime['mri_scanner'])
            temp_file.write("\nSTUDY NAME:%s" % study_name)
            temp_file.write("\nEXAMINATION DATE:%s" % acqu_date)
            temp_file.write("\nPATIENT REFERENCE:%s" % patient_ref)
            temp_file.write("\nPATIENT SEX:%s" % patient_sex)
            temp_file.write("\nPATIENT AGE:%s" % patient_age)
            temp_file.write("\nSOFTWARES:Python %s" % (version.split()[0]))
            from populse_mia import info as mia_info
            temp_file.write("\n : populse_mia %s" % mia_info.__version__)
            from capsul import info as capsul_info
            temp_file.write("\n : capsul %s" % capsul_info.__version__)
            from nipype import info as nipype_info
            temp_file.write("\n : nipype %s" % nipype_info.__version__)
            from mia_processes import info as mia_processes_info
            temp_file.write("\n : mia_processes %s" % mia_processes_info.__version__)

            if platform.system() == 'Darwin':
                os_sys = 'Mac OS X'

            else:
                os_sys = platform.system()

            temp_file.write("\n : Operating System %s %s" % (os_sys, platform.release()))

            cover_data = recupCover(temp_file)

        cover_data = Table(cover_data, [48*mm, 97*mm])  # colWidths, rowHeights

        output = self.report

        page = SimpleDocTemplate(output,
                                 pagesize = portrait(A4),
                                 rightMargin=20*mm,
                                 leftMargin=20*mm,
                                 topMargin=10*mm,
                                 bottomMargin = 10*mm)                          # Output document template definition for page; margins and pages size.

        styles = getSampleStyleSheet()                                          # Initialises stylesheet with few basic heading and text styles, return a stylesheet object.
        styles.add(ParagraphStyle(name = 'Center', alignment = TA_CENTER))      # ParagraphStyle gives all the attributes available for formatting paragraphs.
        styles.add(ParagraphStyle(name = 'Center2', alignment = TA_CENTER))
        styles['Center2'].leading=24 # If more than 1 line.
        styles.add(ParagraphStyle(name = 'Justify', alignment = TA_JUSTIFY))
        styles.add(ParagraphStyle(name = 'Right', alignment = TA_RIGHT))
        styles.add(ParagraphStyle(name = 'Left', alignment = TA_LEFT))
        styles.add(ParagraphStyle(name = 'Bullet1',
                                  leftIndent=30,
                                  bulletOffsetY=2,
                                  bulletIndent=20,
                                  bulletFontSize=6,
                                  bulletColor='black',
                                  bulletText=u'●')
                                 )
        styles.add(ParagraphStyle(name = 'Bullet2',
                                  leftIndent=60,
                                  bulletOffsetY=1,
                                  bulletIndent=50,
                                  bulletFontSize=6,
                                  bulletColor='black',
                                  bulletText=u'❍')
                                 )

        title = ('<font size=18> <b> {0} report: </b> {1} '
                '<br/> <br/></font>').format(patient_ref,
                                            run_date.strftime('%Y.%m.%d'))     # Title:  patient_ref + report - running date.


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
        rapport.append(header_image)

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        report.append(Paragraph(header_title,
                                 styles['Center']))

        report.append(Spacer(0 * mm, 4 * mm))  # (width, height)  ## In Linux
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

        page.build(rapport, canvasmaker = PageNumCanvas)