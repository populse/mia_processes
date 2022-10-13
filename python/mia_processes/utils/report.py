"""
Module dedicated to report generation

:Contains:
    :Class:
        - Report

"""

##########################################################################
# Populse_mia - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# nibabel import
import nibabel as nib
#from nibabel.processing import resample_from_to

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
from mia_processes.utils import (recupCover, PageNumCanvas, ReportLine,
                                 slice_planes_plot, plot_qi2, dict4runtime_update)

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
import numpy as np

#import matplotlib.pyplot as plt
#import scipy.ndimage as ndi



class Report():
    """blabla"""

    def __init__(self, dict4runtime, IQMs_file, anat_file, norm_anat_file,
                 report_file):
        """blabla"""

        self.anat = anat_file
        self.norm_anat = norm_anat_file

        f = open(IQMs_file)
        self.data = json.load(f)
        f.close()
        # Initialises stylesheet with few basic heading and text styles,
        # return a stylesheet object
        self.styles = getSampleStyleSheet()
        # ParagraphStyle gives all the attributes available for formatting
        # paragraphs
        self.styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
        self.styles.add(ParagraphStyle(name='Center2', alignment=TA_CENTER))
        # If more than 1 line
        self.styles['Center2'].leading = 24
        self.styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
        self.styles.add(ParagraphStyle(name='Right', alignment=TA_RIGHT))
        self.styles.add(ParagraphStyle(name='Left', alignment=TA_LEFT))
        self.styles.add(ParagraphStyle(name='Left2', alignment=TA_LEFT))
        # For left indent 30
        self.styles['Left2'].leftIndent = 30
        self.styles.add(ParagraphStyle(name='Bullet1',
                                       leftIndent=30,
                                       bulletOffsetY=2,
                                       bulletIndent=20,
                                       bulletFontSize=6,
                                       bulletColor='black',
                                       bulletText=u'●'))
        self.styles.add(ParagraphStyle(name='Bullet2',
                                       leftIndent=60,
                                       bulletOffsetY=1,
                                       bulletIndent=50,
                                       bulletFontSize=6,
                                       bulletColor='black',
                                       bulletText=u'❍'))
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
        self.image_cov = Table([[header_image_1, header_image_2]],
                               [70 * mm, 70 * mm])
        self.image_cov.setStyle(TableStyle(
                                       [('ALIGN', (0, 0), (0, 0), 'LEFT'),
                                        ('ALIGN', (-1, 0), (-1, 0), 'RIGHT'),
                                        ('VALIGN', (0, 0), (-1, 0), 'MIDDLE')]))

        self.header_title = ('<sup rise=20 size=9>$</sup>'
                             '<font size=30><b>MRIQ</b></font>'
                             '<font size=11>uality</font>'
                             '<font size=30><b>C</b></font>'
                             '<font size=11>ontrol</font>')

        infos = [dict4runtime[i] for i in ('Site', 'Spectro', 'StudyName',
                                           'AcquisitionDate', 'PatientName',
                                           'Sex', 'Age')]
        infos.insert(4, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        infos.insert(5, self.anat)
        infos.extend((version.split()[0], mia_info.__version__,
                      capsul_info.__version__, nipype_info.__version__,
                      mia_processes_info.__version__,
                      "Operating System {0} {1}".format(
                                        'Mac OS X' if
                                            platform.system() == 'Darwin' else
                                            platform.system(),
                                        platform.release())))
        headers = ["SITE", "MRI SCANNER", "STUDY NAME", "EXAMINATION DATE",
                  "MRIQC CALCULATION DATE", "NAME OF THE INPUT DATA",
                  "PATIENT REFERENCE", "PATIENT SEX", "PATIENT AGE",
                  "SOFTWARES", ]
        headers.extend([" "] *  5)

        cover_data = [[0, 0]]

        for header, info in zip(headers, infos):

            if cover_data == [[0, 0]]:

                if header == " ":
                    cover_data[0][0] = Paragraph("<para align=right>" +
                                                    header + "</para>",
                                                 self.styles["Normal"])

                else:
                    # Class reportlab.platypus.paragraph.Paragraph
                    # with XML-like markup
                    cover_data[0][0] = Paragraph("<para align=right><b>" +
                                                    header + "&nbsp&nbsp&nbsp :"
                                                             "</b></para>",
                                                 self.styles["BodyText"])

                cover_data[0][1] = Paragraph("<para align=justify>" +
                                                info + "</para>",
                                             self.styles["Normal"])

            else:

                if header == " ":
                    temp = [Paragraph("<para align=right>" + header + "</para>",
                                      self.styles["Normal"]),
                            Paragraph("<para align=justify>" + info + "</para>",
                                      self.styles["Normal"])]

                else:
                    temp = [Paragraph("<para align=right><b>" + header +
                                        "&nbsp&nbsp&nbsp :</b></para>",
                                      self.styles["BodyText"]),
                            Paragraph("<para align=justify>" + info + "</para>",
                                      self.styles["Normal"])]

                cover_data.append(temp)

        # colWidths, rowHeights
        self.cover_data = Table(cover_data, [60*mm, 97*mm])
        # Careful: Table use (col, raw) and 0-base for top left start as usual
        # OR -1-base for lower right start
        self.cover_data.setStyle(TableStyle([
                              ('LINEABOVE', (0, 0), (-1, 0), 2, colors.black),
                              ('LINEBELOW', (0, -1), (-1, -1), 2, colors.black),
                              ('VALIGN', (0, 0), (-1, -1), 'MIDDLE')]))
        self.cover_data.hAlign = 'CENTER'
        # Output document template definition for page; margins and pages size
        self.page = SimpleDocTemplate(report_file,
                                      pagesize=portrait(A4),
                                      rightMargin=20*mm,
                                      leftMargin=20*mm,
                                      topMargin=10*mm,
                                      bottomMargin=10*mm)

        self.title = ('<font size=18><b>Anatomical Image-Quality Metrics '
                      'summary report</b></font>')
        self.textDisclaimer = ("<font size=7 ><i><b>DISCLAIMER</b><br/>Mia "
                               "software, from the Populse project, is executed"
                               " in a research environment on anonymized data. "
                               "The conclusions obtained with Mia software are "
                               "an help to diagnose and prognosticate. They do "
                               "not substitute themselves to the clinical care "
                               "of the physicians and remain under their "
                               "responsibilities. Consequently, Populse team is"
                               " not responsible for any direct or indirect "
                               "damages resulting from the use of data, "
                               "informations, or results stemming from the Mia "
                               "software. The user recognizes to use these "
                               "informations under his sole and exclusive "
                               "responsibility.<br/> <br/> <b>DECHARGE DE "
                               "RESPONSABILITE</b><br/>Le logiciel Mia, "
                               "provenant du projet Populse, est exécuté dans "
                               "un environnement de recherche sur des données "
                               "anonymisées. Les conclusions obtenues grâce au "
                               "logiciel Mia sont une aide au diagnostic et au "
                               "pronostic. Elles ne se substituent pas à la "
                               "prise en charge médicale des médecins et "
                               "demeurent sous leurs responsabilités. Par "
                               "conséquent, l'équipe Populse ne peut être tenu "
                               "responsable de dommage direct ou indirect "
                               "résultant de l'utilisation des données, des "
                               "informations ou des résultats issus du logiciel"
                               "Mia. L'utilisateur reconnaît utiliser ces "
                               "informations sous sa seule et entière "
                               "responsabilité.</i> </font>")

        #Canvas creation
        self.report = []

    def make_report(self,):
        """blabla"""

       # First page - cover ##################################################
        #######################################################################

        self.report.append(self.image_cov)
        # width, height
        self.report.append(Spacer(0 * mm, 8 * mm))
        self.report.append(Paragraph(self.header_title, self.styles['Center']))
        self.report.append(Spacer(0 * mm, 10 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(Paragraph(self.title, self.styles['Center']))
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(self.cover_data)
        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(Paragraph(self.textDisclaimer,
                                     self.styles["Justify"]))
        self.report.append(Spacer(0 * mm, 6 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        self.report.append(Paragraph("<font size = 8><sup>$</sup>Esteban O et "
                                       "al., <i>MRIQC: Advancing the Automatic "
                                       "Prediction of Image Quality in MRI from"
                                       " Unseen Sites</i>, PLOS ONE 12(9)"
                                       ":e0184661.</font>",
                                     self.styles['Left']))
        self.report.append(PageBreak())

        # Second page - IQMs ##################################################
        #######################################################################
        self.report.append(Paragraph("<font size = 18 ><b>Image parameters"
                                        "</b></font>",
                                     self.styles['Center']))
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        ### Spatial resolution ################################################
        self.report.append(Paragraph(
                          "<font size = 15 > <b>SPATIAL RESOLUTION</b> </font>",
                          self.styles['Bullet1']))
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(Paragraph("Length - Spacing", styles['Left2']))
        self.report.append(Spacer(0 * mm, 10 * mm))

        # size_x
        self.report.append(self.get_iqms_data('size_x'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # size_y
        self.report.append(self.get_iqms_data('size_y'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # size_z
        self.report.append(self.get_iqms_data('size_z'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # spacing_x
        self.report.append(self.get_iqms_data('spacing_x'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # spacing_y
        self.report.append(self.get_iqms_data('spacing_y'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # spacing_z
        self.report.append(self.get_iqms_data('spacing_z'))
        self.report.append(Spacer(0 * mm, 20 * mm))
        ### Signal-to-noise ###################################################
        self.report.append(Paragraph("<font size = 15 > <b>NOISE</b> </font>",
                                     self.styles['Bullet1']))
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(Paragraph(
            "Impact of noise and/or evaluation of the fitness of a noise model",
            self.styles['Left2']))
        self.report.append(Spacer(0 * mm, 10 * mm))
        # snr_csf
        self.report.append(self.get_iqms_data('snr_csf'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # snr_wm
        self.report.append(self.get_iqms_data('snr_wm'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # snr_gm
        self.report.append(self.get_iqms_data('snr_gm'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # snr_total
        self.report.append(self.get_iqms_data('snr_total'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # snrd_csf
        self.report.append(self.get_iqms_data('snrd_csf'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # snrd_wm
        try:
            report.append(Paragraph(
                         '<font size = 9><sup>$</sup></font><font size = 11><b>'
                         'Dietrich’s SNR for white matter</b></font>: ' +
                         str(round(data.get('snrd_wm'), 2)),
                         styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                         '<font size = 9><sup>$</sup></font><font size = 11><b>'
                         'Dietrich’s SNR for white matter</b></font>: ' +
                         'Not determined',
                         styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # snrd_gm
        try:
            report.append(Paragraph(
                         '<font size = 9><sup>$</sup></font><font size = 11><b>'
                         'Dietrich’s SNR for gray matter</b></font>: ' +
                         str(round(data.get('snrd_gm'), 2)),
                         styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                         '<font size = 9><sup>$</sup></font><font size = 11><b>'
                        'Dietrich’s SNR for gray matter</b></font>: ' +
                        'Not determined',
                        styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # snrd_total
        try:
            report.append(Paragraph(
                         '<font size = 9><sup>$</sup></font><font size = 11><b>'
                         'Dietrich’s SNR for brain parenchyma</b></font>: ' +
                         str(round(data.get('snrd_total'), 2)),
                         styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                         '<font size = 9><sup>$</sup></font><font size = 11><b>'
                         'Dietrich’s SNR for brain parenchyma</b></font>: ' +
                         'Not determined',
                         styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))
        # cnr
        try:
            report.append(Paragraph(
                         '<font size = 9><sup>#</sup></font><font size = 11><b>'
                         'Contrast-to-noise ratio</b></font>: ' +
                         str(round(data.get('cnr'), 2)),
                         styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                         '<font size = 9><sup>#</sup></font><font size = 11><b>'
                        'Contrast-to-noise ratio</b></font>: ' +
                        'Not determined',
                        styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))
        # qi_2
        try:
            report.append(Paragraph(
                         '<font size = 9><sup>&</sup></font><font size = 11><b>'
                         'Mortamet’s quality index 2</b></font>: ' +
                         '{:.2e}'.format(data['qi_2']),
                         styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                         '<font size = 9><sup>&</sup></font><font size = 11><b>'
                         'Mortamet’s quality index 2</b></font>: ' +
                         'Not determined',
                         styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))
        # cjv
        try:
            report.append(Paragraph(
                         '<font size = 9><sup>%</sup></font><font size = 11><b>'
                         'Coefficient of joint variation</b></font>: ' +
                         str(round(data.get('cjv'), 2)),
                         styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                         '<font size = 9><sup>%</sup></font><font size = 11><b>'
                         'Coefficient of joint variation</b></font>: ' +
                         'Not determined',
                         styles['Bullet2']))

        report.append(Spacer(0 * mm, 52 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        report.append(line)
        report.append(Spacer(0 * mm, 2.5 * mm))
        report.append(Paragraph(
            "<font size = 8><sup>$</sup>Dietrich et al., <i>Measurement of SNRs"
            " in MR images: influence of multichannel coils, parallel imaging "
            "and reconstruction filters</i>, JMRI 26(2):375–385, 2007. "
            "Higher values are better.</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>#</sup>Magnotta, VA., & Friedman, L., <i>"
            "Measurement of signal-to-noise and contrast-to-noise in the fBIRN "
            "multicenter imaging study</i>, J Dig Imag 19(2):140-147, 2006. "
            "Higher values are better.</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>&</sup>Mortamet B et al., <i>Automatic "
            "quality assessment in structural brain magnetic resonance imaging"
            "</i>, Mag Res Med 62(2):365-372, 2009. Lower values are better."
            "</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>%</sup>Ganzetti et al., <i>Intensity "
            "inhomogeneity correction of structural MR images: a data-driven "
            "approach to define input algorithm parameters</i>, Front "
            "Neuroinform 10:10, 2016. Lower values are better.</font>",
            styles['Left']))
        report.append(PageBreak())

    def get_iqms_data(self, param):
        """blabla"""

        extra_info = {'size_x': ['Voxel size in X', 2, None],
                      'size_y': ['Voxel size in Y', 2, None],
                      'size_z': ['Voxel size in Z', 2, None],
                      'spacing_x': ['Spacing in X (mm)', 2, None],
                      'spacing_y': ['Spacing in Y (mm)', 2, None],
                      'spacing_z': ['Spacing in Z (mm)', 2, None],
                      'snr_csf': ['Signal-to-noise ratio (SNR) for '
                                    'cerebrospinal fluid',
                                  2, None],
                      'snr_wm': ['SNR for white matter', 2, None],
                      'snr_gm': ['SNR for gray matter', 2, None],
                      'snr_total': ['SNR for brain parenchyma', 2, None],
                      'snrd_csf': ["Dietrich’s SNR for cerebrospinal fluid",
                                   2, '$']
                      }

        if self.data.get(param) is None:

            if extra_info[param][2] is None:
                return Paragraph(
                         "<font size = 11> <b> " + extra_info[param][0] +
                            " </b> </font>: Not determined",
                         self.styles['Bullet2'])
            else:
                return Paragraph(
                        "<font size = 9><sup>" + extra_info[param][2] +
                            "</sup></font><font size = 11><b>" +
                            extra_info[param][0] +
                            "</b></font>: Not determined",
                        self.styles['Bullet2'])

        else:

            if extra_info[param][2] is None:
                return Paragraph(
                         "<font size = 11> <b> " + extra_info[param][0] +
                            " </b> </font>: " +
                            str(round(self.data.get(param),
                                      extra_info[param][1])),
                         self.styles['Bullet2'])

            else:
                return Paragraph(
                        "<font size = 9><sup>" + extra_info[param][2] +
                            "</sup></font><font size = 11><b>" +
                            extra_info[param][0] + "</b></font>: " +
                            str(round(self.data.get(param),
                                      extra_info[param][1])),
                        self.styles['Bullet2'])
