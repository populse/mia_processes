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
from nipype import info as nipype_info
from nipype.interfaces.base import (OutputMultiPath, InputMultiPath, File,
                                    traits, TraitListObject, Undefined,
                                    DictStrStr, Str)
from nipype.interfaces.spm.base import ImageFileSPM

# capsul import
from capsul import info as capsul_info

# mia_processes import:
from mia_processes import info as mia_processes_info
from mia_processes.utils import recupCover, PageNumCanvas, ReportLine, slice_planes_plot

# populse_mia import
from populse_mia import info as mia_info
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

        anat_desc = ('An existing, uncompressed anatomical image file (valid '
                     'extensions: .nii)')

        anat_fig_rows_desc = ('The number of lines for the anatomical slice '
                              'planes plot')

        anat_fig_cols_desc = ('The number of columns for the anatomical slice '
                              'planes plot')

        anat_inf_slice_start_desc = ('The first index displayed in anatomical '
                                     'slice planes plot')

        anat_slices_gap_desc = ('Gap between slices in anatomical slice planes '
                                'plot')

        norm_anat_desc = ('An existing, uncompressed normalised anatomical '
                          'image file (valid extensions: .nii)')

        norm_anat_fig_rows_desc = ('The number of lines for the normalised '
                                   'anatomical slice planes plot')

        norm_anat_fig_cols_desc = ('The number of columns for the normalised '
                                   'anatomical slice planes plot')

        norm_anat_inf_slice_start_desc = ('The first index displayed in '
                                          'normalised anatomical slice planes '
                                          'plot')

        norm_anat_slices_gap_desc = ('Gap between slices in normalised '
                                     'anatomical slice planes plot')


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

        self.add_trait("anat_fig_rows",
                       traits.Int(5,
                                  output=False,
                                  optional=False,
                                  desc=anat_fig_rows_desc))

        self.add_trait("anat_fig_cols",
                       traits.Int(5,
                                  output=False,
                                  optional=False,
                                  desc=anat_fig_cols_desc))

        self.add_trait("anat_inf_slice_start",
                       traits.Either(Undefined,
                                     traits.Int,
                                     output=False,
                                     optional=True,
                                     desc=anat_inf_slice_start_desc))

        self.add_trait("anat_slices_gap",
                       traits.Either(Undefined,
                                     traits.Int,
                                     output=False,
                                     optional=True,
                                     desc=anat_slices_gap_desc))

        self.add_trait("norm_anat",
                       ImageFileSPM(copyfile=False,
                                    output=False,
                                    optional=False,
                                    desc=norm_anat_desc))

        self.add_trait("norm_anat_fig_rows",
                       traits.Int(5,
                                  output=False,
                                  optional=False,
                                  desc=norm_anat_fig_rows_desc))

        self.add_trait("norm_anat_fig_cols",
                       traits.Int(5,
                                  output=False,
                                  optional=False,
                                  desc=norm_anat_fig_cols_desc))

        self.add_trait("norm_anat_inf_slice_start",
                       traits.Either(Undefined,
                                     traits.Int,
                                     output=False,
                                     optional=True,
                                     desc=norm_anat_inf_slice_start_desc))

        self.add_trait("norm_anat_slices_gap",
                       traits.Either(Undefined,
                                     traits.Int,
                                     output=False,
                                     optional=True,
                                     desc=norm_anat_slices_gap_desc))

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

        file_position = (self.anat.find(self.project.getName())
                         + len(self.project.getName()) + 1)
        database_filename = self.anat[file_position:]

        # As we do not have access to the database at the runtime (see #272),
        # we prepare here the data that the run_process_mia method will need
        # via dict4runtime dictionary:
        # patient_name
        if ('PatientName' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['patient_name'] = self.project.session.get_value(
                                                             COLLECTION_CURRENT,
                                                             database_filename,
                                                             'PatientName')

        if self.dict4runtime['patient_name'] in ('', Undefined):
            self.dict4runtime['patient_name'] = "Undefined_name_ref"

        # Generate an output name
        if self.anat and self.anat not in ["<undefined>", traits.Undefined]:
            self.outputs['report'] = os.path.join(
                                   self.output_directory,
                                   self.dict4runtime['patient_name'] +
                                   '_mriqcReport_' +
                                   datetime.now().strftime('%Y_%m_%d_'
                                                           '%H_%M_%S_%f')[:22] +
                                   '.pdf')

        # study_name
        if ('StudyName' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['study_name'] = self.project.session.get_value(
                                                             COLLECTION_CURRENT,
                                                             database_filename,
                                                             'StudyName')

        # acquisition_date
        if ('AcquisitionDate' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['acquisition_date'] = str(
                              self.project.session.get_value(COLLECTION_CURRENT,
                                                             database_filename,
                                                             'AcquisitionDate'))

        # sex
        if ('Sex' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['sex'] = self.project.session.get_value(
                                                             COLLECTION_CURRENT,
                                                             database_filename,
                                                             'Sex')
        # site
        if ('Site' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['site'] = self.project.session.get_value(
                                                             COLLECTION_CURRENT,
                                                             database_filename,
                                                             'Site')
        # mri_scanner
        if ('Spectro' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['mri_scanner'] = self.project.session.get_value(
                                                             COLLECTION_CURRENT,
                                                             database_filename,
                                                             'Spectro')
        # age
        if ('Age' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.dict4runtime['age'] = self.project.session.get_value(
                                                             COLLECTION_CURRENT,
                                                             database_filename,
                                                             'Age')

        # FIXME: Currently the following data is hard-coded. A solution should
        #        be found to retrieve them automatically or to put them in the
        #        input parameters of the brick:
        # Site
        if self.dict4runtime['site'] in ('', Undefined):
            self.dict4runtime['site'] = 'Grenoble University Hospital - CLUNI'

        # MriScanner
        if self.dict4runtime['mri_scanner'] in ('', Undefined):
            self.dict4runtime['mri_scanner'] = 'Philips Achieva 3.0T TX'

        # FIXME: Do we need tags inheritance ?

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(MRIQC_report, self).run_process_mia()

        # Logo populse
        header_image_1 = os.path.join(sources_images.__path__[0],
                                      'Logo_populse_square.jpg')
        # 2156px x 2156px
        header_image_1 = Image(header_image_1, 43.0 * mm, 43.0 * mm)
        # Logo Mia
        header_image_2 = os.path.join(sources_images.__path__[0],
                                      'Logo_populse_mia_HR.jpeg')
        # 2505px x 1843px
        header_image_2 = Image(header_image_2, 54.4 * mm, 40.0 * mm)
        header_title = ('<sup rise=20 size=9>$</sup>'
                        '<font size=30><b>MRIQ</b></font>'
                        '<font size=11>uality</font>'
                        '<font size=30><b>C</b></font>'
                        '<font size=11>ontrol</font>')

        if (('site' in self.dict4runtime) and
                            (self.dict4runtime['site'] not in ('', Undefined))):
            site = self.dict4runtime['site']

        else:
            site = "Undefined site"

        if (('study_name' in self.dict4runtime) and
                      (self.dict4runtime['study_name'] not in ('', Undefined))):
            study_name = self.dict4runtime['study_name']

        else:
            study_name = "Undefined study name"

        if (('acquisition_date' in self.dict4runtime) and
                (self.dict4runtime['acquisition_date'] not in ('', Undefined))):
            acqu_date = self.dict4runtime['acquisition_date']

        else:
            acqu_date = "Undefined study name"

        if (('patient_name' in self.dict4runtime) and
                    (self.dict4runtime['patient_name'] not in ('', Undefined))):
            patient_ref = self.dict4runtime['patient_name']

        else:
            patient_ref = "Undefined patient reference"

        if (('sex' in self.dict4runtime) and
                             (self.dict4runtime['sex'] not in ('', Undefined))):
            patient_sex = self.dict4runtime['sex']

        else:
            patient_sex = "Undefined patient sex"

        if (('age' in self.dict4runtime) and
                             (self.dict4runtime['age'] not in ('', Undefined))):
            patient_age = self.dict4runtime['age']

        else:
            patient_age = "Undefined patient age"

        with tempfile.NamedTemporaryFile() as temp_file:
            temp_file.write(bytes("SITE:{}".format(site), encoding='utf8'))
            temp_file.write(bytes("\nMRI SCANNER: {}".format(
                                              self.dict4runtime['mri_scanner']),
                                  encoding='utf8'))
            temp_file.write(bytes("\nSTUDY NAME: {}".format(study_name),
                                  encoding='utf8'))
            temp_file.write(bytes("\nEXAMINATION DATE: {}".format(acqu_date),
                                  encoding='utf8'))
            temp_file.write(bytes("\nMRIQC CALCULATION DATE: {}".format(
                                  datetime.now().strftime('%Y-%m-%d %H:%M:%S')),
                                  encoding='utf8'))
            temp_file.write(bytes("\nNAME OF THE INPUT DATA: {}".format(
                                                                     self.anat),
                                  encoding='utf8'))
            temp_file.write(bytes("\nPATIENT REFERENCE: {}".format(patient_ref),
                                  encoding='utf8'))
            temp_file.write(bytes("\nPATIENT SEX: {}".format(patient_sex),
                                  encoding='utf8'))
            temp_file.write(bytes("\nPATIENT AGE: {}".format(patient_age),
                                  encoding='utf8'))
            temp_file.write(bytes("\nSOFTWARES:Python {}".format(
                                                            version.split()[0]),
                                  encoding='utf8'))
            temp_file.write(bytes("\n : populse_mia {}".format(
                                                          mia_info.__version__),
                                  encoding='utf8'))
            temp_file.write(bytes("\n : capsul {}".format(
                                                       capsul_info.__version__),
                                  encoding='utf8'))
            temp_file.write(bytes("\n : nipype {}".format(
                                                       nipype_info.__version__),
                                  encoding='utf8'))
            temp_file.write(bytes("\n : mia_processes {}".format(
                                                mia_processes_info.__version__),
                                  encoding='utf8'))

            if platform.system() == 'Darwin':
                os_sys = 'Mac OS X'

            else:
                os_sys = platform.system()

            temp_file.write(bytes("\n : Operating System {0} {1}".format(
                                                            os_sys,
                                                            platform.release()),
                                  encoding='utf8'))
            temp_file.flush()
            cover_data = recupCover(temp_file.name)

        # colWidths, rowHeights
        cover_data = Table(cover_data, [60*mm, 97*mm])
        # Output document template definition for page; margins and pages size
        page = SimpleDocTemplate(self.report,
                                 pagesize=portrait(A4),
                                 rightMargin=20*mm,
                                 leftMargin=20*mm,
                                 topMargin=10*mm,
                                 bottomMargin=10*mm)
        # Initialises stylesheet with few basic heading and text styles,
        # return a stylesheet object
        styles = getSampleStyleSheet()
        # ParagraphStyle gives all the attributes available for formatting
        # paragraphs
        styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
        styles.add(ParagraphStyle(name='Center2', alignment=TA_CENTER))
        # If more than 1 line
        styles['Center2'].leading = 24
        styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
        styles.add(ParagraphStyle(name='Right', alignment=TA_RIGHT))
        styles.add(ParagraphStyle(name='Left', alignment=TA_LEFT))
        styles.add(ParagraphStyle(name='Left2', alignment=TA_LEFT))
        # For left indent 30
        styles['Left2'].leftIndent = 30
        styles.add(ParagraphStyle(name='Bullet1',
                                  leftIndent=30,
                                  bulletOffsetY=2,
                                  bulletIndent=20,
                                  bulletFontSize=6,
                                  bulletColor='black',
                                  bulletText=u'●'))
        styles.add(ParagraphStyle(name='Bullet2',
                                  leftIndent=60,
                                  bulletOffsetY=1,
                                  bulletIndent=50,
                                  bulletFontSize=6,
                                  bulletColor='black',
                                  bulletText=u'❍'))
        title = ('<font size=18> <b>Image-Quality Metrics '
                 'summary report</b></font>')


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
        image_cov.setStyle(TableStyle([('ALIGN', (0, 0), (0, 0), 'LEFT'),
                                       ('ALIGN', (-1, 0), (-1, 0), 'RIGHT'),
                                       ('VALIGN', (0, 0), (-1, 0), 'MIDDLE')]))
        report.append(image_cov)
        # width, height
        report.append(Spacer(0 * mm, 8 * mm))
        report.append(Paragraph(header_title, styles['Center']))
        report.append(Spacer(0 * mm, 10 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)
        report.append(Spacer(0 * mm, 10 * mm))
        report.append(Paragraph(title, styles['Center']))
        report.append(Spacer(0 * mm, 10 * mm))
        # Carefull: Table use (col, raw) and 0-base for top left start as usual
        # OR -1-base for lower right start
        cover_data.setStyle(TableStyle([
                              ('LINEABOVE', (0, 0), (-1, 0), 2, colors.black),
                              ('LINEBELOW', (0, -1), (-1, -1), 2, colors.black),
                              ('VALIGN', (0, 0), (-1, -1), 'MIDDLE')]))
        cover_data.hAlign = 'CENTER'
        report.append(cover_data)
        report.append(Spacer(0 * mm, 5 * mm))
        report.append(Paragraph(textDisclaimer, styles["Justify"]))
        report.append(Spacer(0 * mm, 6 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        report.append(line)
        report.append(Spacer(0 * mm, 2.5 * mm))
        report.append(Paragraph("<font size = 8><sup>$</sup>Esteban O et al., "
                                "<i>MRIQC: Advancing the Automatic Prediction "
                                "of Image Quality in MRI from Unseen Sites</i>"
                                ", PLOS ONE 12(9):e0184661.</font>",
                                styles['Left']))
        report.append(PageBreak())

        # Second page - IQMs
        f = open(self.IQMs_file)
        data = json.load(f)
        f.close()
        report.append(Paragraph("<font size = 18 ><b>Image parameters"
                                "</b></font>",
                                styles['Center']))
        report.append(Spacer(0 * mm, 4 * mm))

        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 20 * mm))

        ### Spatial resolution ################################################
        report.append(Paragraph(
            "<font size = 15 > <b>SPATIAL RESOLUTION</b> </font>",
            styles['Bullet1']))

        report.append(Spacer(0 * mm, 2 * mm))  # (width, height)

        report.append(Paragraph("Length - Spacing", styles['Left2']))

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        # size_x
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Voxel size in X </b> </font>: ' + \
                str(round(data.get('size_x'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Voxel size in X </b> </font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # size_y
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Voxel size in Y </b> </font>: ' + \
                str(round(data.get('size_y'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Voxel size in Y </b> </font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # size_z
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Voxel size in Z </b> </font>: ' + \
                str(round(data.get('size_z'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Voxel size in Z </b> </font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # spacing_x
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Spacing in X (mm)</b> </font>: ' + \
                str(round(data.get('spacing_x'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Spacing in X (mm)</b> </font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # spacing_y
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Spacing in Y (mm)</b> </font>: ' + \
                str(round(data.get('spacing_y'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Spacing in Y (mm)</b> </font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # spacing_z
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Spacing in Z (mm)</b> </font>: ' + \
                str(round(data.get('spacing_z'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Spacing in Z (mm)</b> </font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 20 * mm))

        ### Signal-to-noise ###################################################
        report.append(Paragraph(
            "<font size = 15 > <b>NOISE</b> </font>",
            styles['Bullet1']))

        report.append(Spacer(0 * mm, 2 * mm))  # (width, height)

        report.append(Paragraph(
            "Impact of noise and/or evaluation of the fitness of a noise model",
            styles['Left2']))

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        # snr_csf
        try:
            report.append(Paragraph(
                '<font size = 11> <b> Signal-to-noise ratio (SNR) for cerebrospinal fluid </b></font>: ' + \
                str(round(data.get('snr_csf'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> Signal-to-noise ratio (SNR) for cerebrospinal fluid </b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # snr_wm
        try:
            report.append(Paragraph(
                '<font size = 11> <b> SNR for white matter </b></font>: ' + \
                str(round(data.get('snr_wm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> SNR for white matter </b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # snr_gm
        try:
            report.append(Paragraph(
                '<font size = 11> <b> SNR for gray matter </b></font>: ' + \
                str(round(data.get('snr_gm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> SNR for gray matter </b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # snr_total
        try:
            report.append(Paragraph(
                '<font size = 11> <b> SNR for brain parenchyma </b></font>: ' + \
                str(round(data.get('snr_total'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11> <b> SNR for brain parenchyma </b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # snrd_csf
        try:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for cerebrospinal fluid</b></font>: ' + \
                str(round(data.get('snrd_csf'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for cerebrospinal fluid</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # snrd_wm
        try:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for white matter</b></font>: ' + \
                str(round(data.get('snrd_wm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for white matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # snrd_gm
        try:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for gray matter</b></font>: ' + \
                str(round(data.get('snrd_gm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # snrd_total
        try:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for brain parenchyma</b></font>: ' + \
                str(round(data.get('snrd_total'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Dietrich’s SNR for brain parenchyma</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # cnr
        try:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>Contrast-to-noise ratio</b></font>: ' + \
                str(round(data.get('cnr'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>Contrast-to-noise ratio</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # qi_2
        try:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11><b>Mortamet’s quality index 2</b></font>: ' + \
                '{:.2e}'.format(data['qi_2']),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11><b>Mortamet’s quality index 2</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # cjv
        try:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Coefficient of joint variation</b></font>: ' + \
                str(round(data.get('cjv'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Coefficient of joint variation</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 52 * mm))  # (width, height)

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

        ### Spatial distribution ##############################################
        report.append(Paragraph(
            "<font size = 15 > <b>SPATIAL DISTRIBUTION</b> </font>",
            styles['Bullet1']))

        report.append(Spacer(0 * mm, 2 * mm))  # (width, height)

        report.append(Paragraph(
            "Information theory to evaluate the spatial distribution of information",
            styles['Left2']))

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        # fber
        try:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>FBER</b></font>: ' + \
                str(round(data.get('fber'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>FBER</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # efc
        try:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>EFC</b></font>: ' + \
                str(round(data.get('efc'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>EFC</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 20 * mm))

        ### Artifacts #########################################################
        report.append(Paragraph(
            "<font size = 15 > <b>ARTIFACTS</b> </font>",
            styles['Bullet1']))

        report.append(Spacer(0 * mm, 2 * mm))  # (width, height)

        report.append(Paragraph(
            "Estimates artefacts and signal leakage due to rapid movement (e.g. eye movement or blood vessel pulsation)",
            styles['Left2']))

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        # wm2max
        try:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11><b>White-matter to maximum intensity ratio </b></font>: ' + \
                str(round(data.get('wm2max'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11><b>White-matter to maximum intensity ratio </b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # qi_1
        try:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Mortamet’s quality index 1</b></font>: ' + \
                '{:.2e}'.format(data['qi_1']),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>$</sup></font><font size = 11><b>Mortamet’s quality index 1</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # inu_range
        try:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Bias field range (95th percentile - 5th percentile)</b></font>: ' + \
                str(round(data.get('inu_range'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Bias field range (95th percentile - 5th percentile)</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # inu_med
        try:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Bias field median</b></font>: ' + \
                str(round(data.get('inu_med'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Bias field median</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 100 * mm))  # (width, height)

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
        report.append(Paragraph(
            "<font size = 8><sup>*</sup>Tustison NJ et al., <i>N4ITK: improved N3 bias correction</i>, IEEE Trans Med Imag, 29(6):1310-20, 2010. Median closer to 1 and range closer to 0 are better.</font>",
            styles['Left']))

        report.append(PageBreak())

        # fourth page - IQMs

        report.append(Paragraph(
            "<font size = 18 > <b> Image parameters <br/> </b> </font>",
            styles['Center']))

        report.append(Spacer(0 * mm, 4 * mm))  # (width, height)

        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 20 * mm))

        ### Tissues Quality ##############################################
        report.append(Paragraph(
            "<font size = 15 > <b>TISSUES QUALITY</b> </font>",
            styles['Bullet1']))

        report.append(Spacer(0 * mm, 2 * mm))  # (width, height)

        report.append(Paragraph(
            "Metrics that do not fall into the above categories: statistical properties of tissue distributions, volume overlap of tissues, image harpness/blurriness, etc.",
            styles['Left2']))

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        # summary_csf_mean
        try:
            report.append(Paragraph(
                '<font size = 11><b>Mean of the distribution of cerebrospinal fluid</b></font>: ' + \
                str(round(data.get('summary_csf_mean'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Mean of the distribution of cerebrospinal fluid</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_csf_stdv
        try:
            report.append(Paragraph(
                '<font size = 11><b>Standard deviation of the distribution of cerebrospinal fluid</b></font>: ' + \
                str(round(data.get('summary_csf_stdv'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Standard deviation of the distribution of cerebrospinal fluid</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_csf_median
        try:
            report.append(Paragraph(
                '<font size = 11><b>Median of the distribution of cerebrospinal fluid</b></font>: ' + \
                str(round(data.get('summary_csf_median'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Median of the distribution of cerebrospinal fluid</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_csf_mad
        try:
            report.append(Paragraph(
                '<font size = 11><b>Median absolute deviation of the distribution of cerebrospinal fluid</b></font>: ' + \
                str(round(data.get('summary_csf_mad'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Median absolute deviation of the distribution of cerebrospinal fluid</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_csf_p95
        try:
            report.append(Paragraph(
                '<font size = 11><b>95% percentile of the distribution of cerebrospinal fluid</b></font>: ' + \
                str(round(data.get('summary_csf_p95'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>95% percentile of the distribution of cerebrospinal fluid</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_csf_p05
        try:
            report.append(Paragraph(
                '<font size = 11><b>5% percentile of the distribution of cerebrospinal fluid</b></font>: ' + \
                str(round(data.get('summary_csf_p05'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>5% percentile of the distribution of cerebrospinal fluid</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_csf_k
        try:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis of the distribution of cerebrospinal fluid</b></font>: ' + \
                str(round(data.get('summary_csf_k'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis of the distribution of cerebrospinal fluid</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_csf_n
        try:
            report.append(Paragraph(
                '<font size = 11><b>Number of voxels in the distribution of cerebrospinal fluid</b></font>: ' + \
                '{:.2e}'.format(data['summary_csf_n']),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Number of voxels in the distribution of cerebrospinal fluid</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # summary_gm_mean
        try:
            report.append(Paragraph(
                '<font size = 11><b>Mean of the distribution of gray matter</b></font>: ' + \
                str(round(data.get('summary_gm_mean'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Mean of the distribution of gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_gm_stdv
        try:
            report.append(Paragraph(
                '<font size = 11><b>Standard deviation of the distribution of gray matter</b></font>: ' + \
                str(round(data.get('summary_gm_stdv'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Standard deviation of the distribution of gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_gm_median
        try:
            report.append(Paragraph(
                '<font size = 11><b>Median of the distribution of gray matter</b></font>: ' + \
                str(round(data.get('summary_gm_median'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Median of the distribution of gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_gm_mad
        try:
            report.append(Paragraph(
                '<font size = 11><b>Median absolute deviation of the distribution of gray matter</b></font>: ' + \
                str(round(data.get('summary_gm_mad'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Median absolute deviation of the distribution of gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_gm_p95
        try:
            report.append(Paragraph(
                '<font size = 11><b>95% percentile of the distribution of gray matter</b></font>: ' + \
                str(round(data.get('summary_gm_p95'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>95% percentile of the distribution of gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_gm_p05
        try:
            report.append(Paragraph(
                '<font size = 11><b>5% percentile of the distribution of gray matter</b></font>: ' + \
                str(round(data.get('summary_gm_p05'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>5% percentile of the distribution of gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_gm_k
        try:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis of the distribution of gray matter</b></font>: ' + \
                str(round(data.get('summary_gm_k'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis of the distribution of gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_gm_n
        try:
            report.append(Paragraph(
                '<font size = 11><b>Number of voxels in the distribution of gray matter</b></font>: ' + \
                '{:.2e}'.format(data['summary_gm_n']),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Number of voxels in the distribution of gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # summary_wm_mean
        try:
            report.append(Paragraph(
                '<font size = 11><b>Mean of the distribution of white matter</b></font>: ' + \
                str(round(data.get('summary_wm_mean'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Mean of the distribution of white matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_wm_stdv
        try:
            report.append(Paragraph(
                '<font size = 11><b>Standard deviation of the distribution of white matter</b></font>: ' + \
                str(round(data.get('summary_wm_stdv'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Standard deviation of the distribution of white matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_wm_median
        try:
            report.append(Paragraph(
                '<font size = 11><b>Median of the distribution of white matter</b></font>: ' + \
                str(round(data.get('summary_wm_median'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Median of the distribution of white matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_wm_mad
        try:
            report.append(Paragraph(
                '<font size = 11><b>Median absolute deviation of the distribution of white matter</b></font>: ' + \
                str(round(data.get('summary_wm_mad'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Median absolute deviation of the distribution of white matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_wm_p95
        try:
            report.append(Paragraph(
                '<font size = 11><b>95% percentile of the distribution of white matter</b></font>: ' + \
                str(round(data.get('summary_wm_p95'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>95% percentile of the distribution of white matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_wm_p05
        try:
            report.append(Paragraph(
                '<font size = 11><b>5% percentile of the distribution of white matter</b></font>: ' + \
                str(round(data.get('summary_wm_p05'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>5% percentile of the distribution of white matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_wm_k
        try:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis of the distribution of white matter</b></font>: ' + \
                str(round(data.get('summary_wm_k'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis of the distribution of white matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_wm_n
        try:
            report.append(Paragraph(
                '<font size = 11><b>Number of voxels in the distribution of white matter</b></font>: ' + \
                '{:.2e}'.format(data['summary_wm_n']),

                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Number of voxels in the distribution of white matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # summary_bg_mean
        try:
            report.append(Paragraph(
                '<font size = 11><b>Mean of the distribution of background</b></font>: ' + \
                str(round(data.get('summary_bg_mean'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Mean of the distribution of background</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_bg_stdv
        try:
            report.append(Paragraph(
                '<font size = 11><b>Standard deviation of the distribution of background</b></font>: ' + \
                str(round(data.get('summary_bg_stdv'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Standard deviation of the distribution of background</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_bg_median
        try:
            report.append(Paragraph(
                '<font size = 11><b>Median of the distribution of background</b></font>: ' + \
                str(round(data.get('summary_bg_median'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Median of the distribution of background</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_bg_mad
        try:
            report.append(Paragraph(
                '<font size = 11><b>Median absolute deviation of the distribution of background</b></font>: ' + \
                str(round(data.get('summary_bg_mad'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Median absolute deviation of the distribution of background</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_bg_p95
        try:
            report.append(Paragraph(
                '<font size = 11><b>95% percentile of the distribution of background</b></font>: ' + \
                str(round(data.get('summary_bg_p95'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>95% percentile of the distribution of background</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_bg_p05
        try:
            report.append(Paragraph(
                '<font size = 11><b>5% percentile of the distribution of background</b></font>: ' + \
                str(round(data.get('summary_bg_p05'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>5% percentile of the distribution of background</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_bg_k
        try:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis of the distribution of background</b></font>: ' + \
                str(round(data.get('summary_bg_k'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis of the distribution of background</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # summary_bg_n
        try:
            report.append(Paragraph(
                '<font size = 11><b>Number of voxels in the distribution of background</b></font>: ' + \
                '{:.2e}'.format(data['summary_bg_n']),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Number of voxels in the distribution of background</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # fwhm_x
        try:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of the distribution in units of voxels, along x dimension</b></font>: ' + \
                str(round(data.get('fwhm_x'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of the distribution in units of voxels, along x dimension</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # fwhm_y
        try:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of the distribution in units of voxels, along y dimension</b></font>: ' + \
                str(round(data.get('fwhm_y'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of the distribution in units of voxels, along y dimension</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # fwhm_z
        try:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of the distribution in units of voxels, along z dimension</b></font>: ' + \
                str(round(data.get('fwhm_z'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of the distribution in units of voxels, along z dimension</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # fwhm_avg
        try:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of the distribution in units of voxels, average</b></font>: ' + \
                str(round(data.get('fwhm_avg'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of the distribution in units of voxels, average</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 7 * mm))  # (width, height)

        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        report.append(Paragraph(
            "<font size = 8><sup>*</sup>K is always ≥ -2. If the distribution is Gaussian, K = 0. If a distribution has less weight on its center and tails compared to a Gaussian of the same variance, then K < 0. If the distribution has more weight on its center and tails, then K > 0.</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>#</sup>Forman SD et al., <i>Improved assessment of significant activation in functional magnetic resonance imaging(fMRI): use of a cluster - size threshold</i>, Magn.Reson.Med. 33(5), 636–647, 1995. Lower values are better, higher values indicate a blurrier image.</font>",
            styles['Left']))

        report.append(PageBreak())

        # Fifth page page - IQMs

        report.append(Paragraph(
            "<font size = 18 > <b> Image parameters <br/> </b> </font>",
            styles['Center']))

        report.append(Spacer(0 * mm, 4 * mm))  # (width, height)

        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 20 * mm))

        ### Tissues Quality ##############################################
        report.append(Paragraph(
            "<font size = 15 > <b>TISSUES QUALITY</b> </font>",
            styles['Bullet1']))

        report.append(Spacer(0 * mm, 2 * mm))  # (width, height)

        report.append(Paragraph(
            "Metrics that do not fall into the above categories: statistical properties of tissue distributions, volume overlap of tissues, image harpness/blurriness, etc.",
            styles['Left2']))

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        # icvs_csf
        try:
            report.append(Paragraph(
                '<font size = 11><b>Intracranial volume fraction of cerebrospinal fluid</b></font>: ' + \
                str(round(data.get('icvs_csf'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Intracranial volume fraction of cerebrospinal fluid</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # icvs_gm
        try:
            report.append(Paragraph(
                '<font size = 11><b>Intracranial volume fraction of gray matter</b></font>: ' + \
                str(round(data.get('icvs_gm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Intracranial volume fraction of gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # icvs_wm
        try:
            report.append(Paragraph(
                '<font size = 11><b>Intracranial volume fraction of white matter</b></font>: ' + \
                str(round(data.get('icvs_wm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 11><b>Intracranial volume fraction of white matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # rpve_csf
        try:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Residual partial voluming error of cerebrospinal fluid</b></font>: ' + \
                str(round(data.get('rpve_csf'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Residual partial voluming error of cerebrospinal fluid</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # rpve_gm
        try:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Residual partial voluming error of gray matter</b></font>: ' + \
                str(round(data.get('rpve_gm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Residual partial voluming error of gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # rpve_wm
        try:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Residual partial voluming error of white matter</b></font>: ' + \
                str(round(data.get('rpve_wm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Residual partial voluming error of gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        # tpm_overlap_csf
        try:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11><b>Overlap of the tissue probability map and the corresponding ICBM nonlinear-asymmetric 2009c template, for cerebrospinal fluid</b></font>: ' + \
                str(round(data.get('tpm_overlap_csf'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11><b>Overlap of the tissue probability map and the corresponding ICBM nonlinear-asymmetric 2009c template, for cerebrospinal fluid</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # tpm_overlap_gm
        try:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11><b>Overlap of the tissue probability map and the corresponding ICBM nonlinear-asymmetric 2009c template, for gray matter</b></font>: ' + \
                str(round(data.get('tpm_overlap_gm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11><b>Overlap of the tissue probability map and the corresponding ICBM nonlinear-asymmetric 2009c template, for gray matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))  # (width, height)

        # tpm_overlap_wm
        try:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11><b>Overlap of the tissue probability map and the corresponding ICBM nonlinear-asymmetric 2009c template, for  white matter</b></font>: ' + \
                str(round(data.get('tpm_overlap_wm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>&</sup></font><font size = 11><b>Overlap of the tissue probability map and the corresponding ICBM nonlinear-asymmetric 2009c template, for  white matter</b></font>: ' + \
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 147 * mm))  # (width, height)

        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 2.5 * mm))  # (width, height)

        report.append(Paragraph(
            "<font size = 8><sup>%</sup>Lower values are better.</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>&</sup>Higher values are better.</font>",
            styles['Left']))
        report.append(PageBreak())

        # Sixth page - slice planes display - Raw anatomic

        report.append(Paragraph("<font size = 18 > <b>MRI axial slice planes display</b> </font>",
                                     styles['Center']))

        report.append(Spacer(0*mm, 4*mm)) # (width, height) ## In Linux

        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        report.append(Paragraph(
            "<font size = 14 >Raw anatomic images</font>",
            styles['Center']))

        report.append(Spacer(0*mm, 10*mm))  # (width, height)
        report.append(Paragraph('<font size = 9 > <i> "Neurological" convention, the left side of the image corresponds to the left side of the brain. </i> <br/> </font>',
                                     styles['Center']))
        report.append(Spacer(0*mm, 1*mm))  # (width, height)

        tmpdir = tempfile.TemporaryDirectory()

        slices_image = slice_planes_plot(self.anat, self.anat_fig_rows, self.anat_fig_cols, inf_slice_start=self.anat_inf_slice_start, slices_gap=self.anat_slices_gap, cmap="Greys_r", out_dir=tmpdir.name)
        # reminder: A4 == 210mmx297mm
        slices_image = Image(slices_image, width=7.4803 * inch, height=9.0551 * inch)
        slices_image.hAlign = 'CENTER'
        report.append(slices_image)

        report.append(PageBreak())

        # Seventh page - slice planes display - Normalised anatomic (MNI)

        report.append(Paragraph("<font size = 18 > <b>MRI axial slice planes display</b> </font>",
                                     styles['Center']))

        report.append(Spacer(0*mm, 4*mm)) # (width, height) ## In Linux

        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)

        report.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        report.append(Paragraph(
            "<font size = 14 >Normalised anatomic (MNI)</font>",
            styles['Center']))

        report.append(Spacer(0*mm, 10*mm))  # (width, height)
        report.append(Paragraph('<font size = 9 > <i> "Neurological" convention, the left side of the image corresponds to the left side of the brain. </i> <br/> </font>',
                                     styles['Center']))
        report.append(Spacer(0*mm, 1*mm))  # (width, height)

        #tmpdir = tempfile.TemporaryDirectory()

        #slices_image = slice_planes_plot(self.anat, self.anat_fig_rows, self.anat_fig_cols, inf_slice_start=self.anat_inf_slice_start, slices_gap=self.anat_slices_gap, cmap="Greys_r", out_dir=tmpdir.name)
        # reminder: A4 == 210mmx297mm
        #slices_image = Image(slices_image, width=7.4803 * inch, height=9.0551 * inch)
        #slices_image.hAlign = 'CENTER'
        #report.append(slices_image)
        page.build(report, canvasmaker=PageNumCanvas)
        tmpdir.cleanup()
