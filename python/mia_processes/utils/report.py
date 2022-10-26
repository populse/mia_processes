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
from mia_processes.utils import (PageNumCanvas, ReportLine,
                                 slice_planes_plot, plot_qi2,
                                 dict4runtime_update)

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

    def __init__(self, report_file, dict4runtime, **kwargs):
        """blabla"""

        self.dict4runtime = dict4runtime
        ref_exp = "Undefined"

        for key in kwargs:
            setattr(self, key, kwargs[key])

        if "IQMs_file" in kwargs:

            with open(self.IQMs_file) as f:
                self.iqms_data = json.load(f)

            if "anat" in kwargs:
                self.make_report = self.mriqc_anat_make_report
                ref_exp = self.anat
                self.title = ('<font size=18><b>Anatomical Image-Quality '
                              'Metrics summary report</b></font>')

            elif "func" in kwargs:
                self.make_report = self.mriqc_func_make_report
                ref_exp = self.func
                self.title = ('<font size=18><b>Functional Image-Quality '
                              'Metrics summary report</b></font>')

            self.header_title = ('<sup rise=20 size=9>$</sup>'
                                 '<font size=30><b>MRIQ</b></font>'
                                 '<font size=11>uality</font>'
                                 '<font size=30><b>C</b></font>'
                                 '<font size=11>ontrol</font>')

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

        infos = [self.dict4runtime[i] for i in (
                                           'Site', 'Spectro', 'StudyName',
                                           'AcquisitionDate', 'PatientName',
                                           'Sex', 'Age')]
        infos.insert(4, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        infos.insert(5, ref_exp)
        infos.extend(("Python " + version.split()[0],
                      "Populse_mia " + mia_info.__version__,
                      "Capsul " + capsul_info.__version__,
                      "Nipype " + nipype_info.__version__,
                      "Mia_processes " + mia_processes_info.__version__,
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
                            Paragraph("<para align=justify><font size = 7>" +
                                      info +
                                      "</font></para>" if
                                      header == "NAME OF THE INPUT DATA" else
                                      "<para align=justify>" + info + "</para>",
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

        # Output document template definition for page; margins and pages size
        self.page = SimpleDocTemplate(report_file,
                                      pagesize=portrait(A4),
                                      rightMargin=20*mm,
                                      leftMargin=20*mm,
                                      topMargin=10*mm,
                                      bottomMargin=10*mm)

        #Canvas creation
        self.report = []

    def mriqc_func_make_report(self,):
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
        self.report.append(Spacer(0 * mm, 6 * mm))
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
        ### Spatial and temporal resolution ###################################
        self.report.append(Paragraph(
            "<font size = 15 > <b>SPATIAL AND TEMPORAL RESOLUTION</b> </font>",
            self.styles['Bullet1']))
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(Paragraph("Length - Spacing", self.styles['Left2']))
        self.report.append(Spacer(0 * mm, 10 * mm))

        # size_x
        self.report.append(self.get_iqms_data('size_x'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # size_y
        self.report.append(self.get_iqms_data('size_y'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # size_z
        self.report.append(self.get_iqms_data('size_z'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # size_t
        self.report.append(self.get_iqms_data('size_t'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # spacing_x
        self.report.append(self.get_iqms_data('spacing_x'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # spacing_y
        self.report.append(self.get_iqms_data('spacing_y'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # spacing_z
        self.report.append(self.get_iqms_data('spacing_z'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # spacing_tr
        self.report.append(self.get_iqms_data('spacing_tr'))
        self.report.append(Spacer(0 * mm, 20 * mm))
        ### Measures for artifacts and other ##################################
        self.report.append(Paragraph("<font size = 15 > <b>MEASURES FOR "
                                         "ARTIFACTS AND OTHER</b> </font>",
                                     self.styles['Bullet1']))
        self.report.append(Spacer(0 * mm, 10 * mm))
        # gsr_x
        self.report.append(self.get_iqms_data('gsr_x'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # gsr_y
        self.report.append(self.get_iqms_data('gsr_y'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # fd_mean
        self.report.append(self.get_iqms_data('fd_mean'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fd_num
        self.report.append(self.get_iqms_data('fd_num'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fd_perc
        self.report.append(self.get_iqms_data('fd_perc'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # dummyTRs
        self.report.append(self.get_iqms_data('dummyTRs'))
        self.report.append(Spacer(0 * mm, 90 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("gsr_x", "fd_mean"):
            self.report.append(Paragraph(
                "<font size = 8><sup>" +
                self.dict4runtime['extra_info'][param][2] +
                "</sup>" + self.dict4runtime['extra_info'][param][3] +
                "</font>",
                self.styles['Left']))

        self.report.append(PageBreak())

        # Third page - IQMs ###################################################
        #######################################################################
        self.report.append(Paragraph("<font size = 18 ><b>Image parameters"
                                     "</b></font>",
                                     self.styles['Center']))
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        ### Measures for the spatial information ##############################
        self.report.append(Paragraph("<font size = 15 ><b>MEASURES FOR THE "
                                         "SPATIAL INFORMATION</b></font>",
                                     self.styles['Bullet1']))
        self.report.append(Spacer(0 * mm, 10 * mm))
        # summary_bg_mean
        self.report.append(self.get_iqms_data('summary_bg_mean'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_stdv
        self.report.append(self.get_iqms_data('summary_bg_stdv'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_median
        self.report.append(self.get_iqms_data('summary_bg_median'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_mad
        self.report.append(self.get_iqms_data('summary_bg_mad'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_p95
        self.report.append(self.get_iqms_data('summary_bg_p95'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_p05
        self.report.append(self.get_iqms_data('summary_bg_p05'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_k
        self.report.append(self.get_iqms_data('summary_bg_k'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_n
        self.report.append(self.get_iqms_data('summary_bg_n'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # summary_fg_mean
        self.report.append(self.get_iqms_data('summary_fg_mean'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_stdv
        self.report.append(self.get_iqms_data('summary_fg_stdv'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_median
        self.report.append(self.get_iqms_data('summary_fg_median'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_mad
        self.report.append(self.get_iqms_data('summary_fg_mad'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_p95
        self.report.append(self.get_iqms_data('summary_fg_p95'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_p05
        self.report.append(self.get_iqms_data('summary_fg_p05'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_k
        self.report.append(self.get_iqms_data('summary_fg_k'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_n
        self.report.append(self.get_iqms_data('summary_fg_n'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # snr
        self.report.append(self.get_iqms_data('snr'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # fber
        self.report.append(self.get_iqms_data('fber'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # efc
        self.report.append(self.get_iqms_data('efc'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # fwhm_x
        self.report.append(self.get_iqms_data('fwhm_x'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_y
        self.report.append(self.get_iqms_data('fwhm_y'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_z
        self.report.append(self.get_iqms_data('fwhm_z'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_avg
        self.report.append(self.get_iqms_data('fwhm_avg'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        self.report.append(Spacer(0 * mm, 66 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ('summary_bg_k', "fber", "efc", "fwhm_x"):
            self.report.append(Paragraph(
                "<font size = 8><sup>" +
                self.dict4runtime['extra_info'][param][2] +
                "</sup>" + self.dict4runtime['extra_info'][param][3] +
                "</font>",
                self.styles['Left']))

        self.report.append(PageBreak())

        # fourth page - IQMs ##################################################
        #######################################################################
        self.report.append(Paragraph("<font size = 18 ><b>Image parameters"
                                     "</b></font>",
                                     self.styles['Center']))
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        ### Measures for the temporal information #############################
        self.report.append(Paragraph("<font size = 15 ><b>MEASURES FOR THE "
                                         "TEMPORAL INFORMATION</b></font>",
                                     self.styles['Bullet1']))
        self.report.append(Spacer(0 * mm, 10 * mm))
        # dvars_nstd
        self.report.append(self.get_iqms_data('dvars_nstd'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # dvars_std
        self.report.append(self.get_iqms_data('dvars_std'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # dvars_vstd
        self.report.append(self.get_iqms_data('dvars_vstd'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # tsnr
        self.report.append(self.get_iqms_data('tsnr'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # gcor
        self.report.append(self.get_iqms_data('gcor'))
        self.report.append(Spacer(0 * mm, 170 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("dvars_nstd", "tsnr", "gcor"):
            self.report.append(Paragraph(
                "<font size = 8><sup>" +
                self.dict4runtime['extra_info'][param][2] +
                "</sup>" + self.dict4runtime['extra_info'][param][3] +
                "</font>",
                self.styles['Left']))

        self.report.append(PageBreak())

        # Fifth - slice planes display - Raw functional #######################
        #######################################################################
        self.report.append(Paragraph("<font size = 18 > <b>MRI axial slice "
                                     "planes display</b> </font>",
                                     self.styles['Center']))
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(Paragraph("<font size = 14 >Raw functional"
                                         "images (1<sup>st</sup> dynamic)"
                                         "</font>",
                                     self.styles['Center']))
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(Paragraph('<font size = 9 ><i>"Neurological" '
                                     'convention, the left side of the '
                                     'image corresponds to the left side of '
                                     'the brain.</i></font>',
                                     self.styles['Center']))
        self.report.append(Spacer(0 * mm, 1 * mm))
        tmpdir = tempfile.TemporaryDirectory()
        slices_image = slice_planes_plot(
            self.func, self.func_fig_rows,
            self.func_fig_cols,
            inf_slice_start=self.func_inf_slice_start,
            slices_gap=self.func_slices_gap,
            cmap="Greys_r", out_dir=tmpdir.name)
        # reminder: A4 == 210mmx297mm
        slices_image = Image(slices_image,
                             width=7.4803 * inch, height=9.0551 * inch)
        slices_image.hAlign = 'CENTER'
        self.report.append(slices_image)
        self.report.append(PageBreak())

        # Sixth page - slice planes display - Normalised functional (MNI) #####
        #######################################################################
        self.report.append(Paragraph("<font size = 18 > <b>MRI axial slice "
                                     "planes display</b> </font>",
                                     self.styles['Center']))
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(Paragraph("<font size = 14 >Normalised mean "
                                         "functional (MNI)</font>",
                                     self.styles['Center']))

        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(Paragraph('<font size = 9 > <i> "Neurological" '
                                     'convention, the left side of the '
                                     'image corresponds to the left side of '
                                     'the brain. </i> <br/> </font>',
                                     self.styles['Center']))
        self.report.append(Spacer(0 * mm, 1 * mm))
        slices_image = slice_planes_plot(
            self.norm_func,
            self.norm_func_fig_rows,
            self.norm_func_fig_cols,
            inf_slice_start=self.norm_func_inf_slice_start,
            slices_gap=self.norm_func_slices_gap,
            cmap="Greys_r", out_dir=tmpdir.name)
        # reminder: A4 == 210mmx297mm
        slices_image = Image(slices_image, width=7.4803 * inch,
                             height=9.0551 * inch)
        slices_image.hAlign = 'CENTER'
        self.report.append(slices_image)
        self.page.build(self.report, canvasmaker=PageNumCanvas)
        tmpdir.cleanup()

    def mriqc_anat_make_report(self,):
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
        self.report.append(Spacer(0 * mm, 6 * mm))
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
        self.report.append(Paragraph("Length - Spacing", self.styles['Left2']))
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
        self.report.append(self.get_iqms_data('snrd_wm'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # snrd_gm
        self.report.append(self.get_iqms_data('snrd_gm'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # snrd_total
        self.report.append(self.get_iqms_data('snrd_total'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # cnr
        self.report.append(self.get_iqms_data('cnr'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # qi_2
        self.report.append(self.get_iqms_data('qi_2'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # cjv
        self.report.append(self.get_iqms_data('cjv'))
        self.report.append(Spacer(0 * mm, 52 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("snrd_csf", "cnr", "qi_2", "cjv"):
            self.report.append(Paragraph(
                "<font size = 8><sup>" +
                    self.dict4runtime['extra_info'][param][2] +
                    "</sup>" + self.dict4runtime['extra_info'][param][3] +
                    "</font>",
                self.styles['Left']))

        self.report.append(PageBreak())

        # Third page - IQMs ###################################################
        #######################################################################
        self.report.append(Paragraph("<font size = 18 ><b>Image parameters"
                                        "</b></font>",
                                     self.styles['Center']))
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        ### Spatial distribution ##############################################
        self.report.append(Paragraph("<font size = 15 ><b>SPATIAL DISTRIBUTION"
                                        "</b></font>",
                                     self.styles['Bullet1']))
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(Paragraph("Information theory to evaluate the "
                                        "spatial distribution of information",
                                     self.styles['Left2']))
        self.report.append(Spacer(0 * mm, 10 * mm))
        # fber
        self.report.append(self.get_iqms_data('fber'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # efc
        self.report.append(self.get_iqms_data('efc'))
        self.report.append(Spacer(0 * mm, 20 * mm))
        ### Artifacts #########################################################
        self.report.append(Paragraph("<font size = 15 ><b>ARTIFACTS"
                                        "</b></font>",
                                     self.styles['Bullet1']))
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(
                      Paragraph("Estimates artefacts and signal leakage due to "
                                   "rapid movement (e.g. eye movement or blood "
                                   "vessel pulsation, etc.)",
                                self.styles['Left2']))

        self.report.append(Spacer(0 * mm, 10 * mm))
        # wm2max
        self.report.append(self.get_iqms_data('wm2max'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # qi_1
        self.report.append(self.get_iqms_data('qi_1'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # inu_range
        self.report.append(self.get_iqms_data('inu_range'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # inu_med
        self.report.append(self.get_iqms_data('inu_med'))
        self.report.append(Spacer(0 * mm, 100 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("fber", "efc", "wm2max", "qi_1", "inu_range"):
            self.report.append(Paragraph(
                "<font size = 8><sup>" +
                    self.dict4runtime['extra_info'][param][2] +
                    "</sup>" + self.dict4runtime['extra_info'][param][3] +
                    "</font>",
                self.styles['Left']))

        self.report.append(PageBreak())

        # fourth page - IQMs ##################################################
        #######################################################################
        self.report.append(Paragraph("<font size = 18 ><b>Image parameters"
                                        "</b></font>",
                                     self.styles['Center']))
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        ### Tissues Quality ###################################################
        self.report.append(Paragraph("<font size = 15 ><b>TISSUES QUALITY"
                                        "</b></font>",
                                     self.styles['Bullet1']))
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(Paragraph("Metrics that do not fall into the above "
                                        "categories: statistical properties of "
                                        "tissue distributions, volume overlap "
                                        "of tissues, image harpness/blurriness,"
                                        " etc.",
                                     self.styles['Left2']))
        self.report.append(Spacer(0 * mm, 10 * mm))
        # summary_csf_mean
        self.report.append(self.get_iqms_data('summary_csf_mean'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_stdv
        self.report.append(self.get_iqms_data('summary_csf_stdv'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_median
        self.report.append(self.get_iqms_data('summary_csf_median'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_mad
        self.report.append(self.get_iqms_data('summary_csf_mad'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_p95
        self.report.append(self.get_iqms_data('summary_csf_p95'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_p05
        self.report.append(self.get_iqms_data('summary_csf_p05'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_k
        self.report.append(self.get_iqms_data('summary_csf_k'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_n
        self.report.append(self.get_iqms_data('summary_csf_n'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # summary_gm_mean
        self.report.append(self.get_iqms_data('summary_gm_mean'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_stdv
        self.report.append(self.get_iqms_data('summary_gm_stdv'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_median
        self.report.append(self.get_iqms_data('summary_gm_median'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_mad
        self.report.append(self.get_iqms_data('summary_gm_mad'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_p95
        self.report.append(self.get_iqms_data('summary_gm_p95'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_p05
        self.report.append(self.get_iqms_data('summary_gm_p05'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_k
        self.report.append(self.get_iqms_data('summary_gm_k'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_n
        self.report.append(self.get_iqms_data('summary_gm_n'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # summary_wm_mean
        self.report.append(self.get_iqms_data('summary_wm_mean'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_stdv
        self.report.append(self.get_iqms_data('summary_wm_stdv'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_median
        self.report.append(self.get_iqms_data('summary_wm_median'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_mad
        self.report.append(self.get_iqms_data('summary_wm_mad'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_p95
        self.report.append(self.get_iqms_data('summary_wm_p95'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_p05
        self.report.append(self.get_iqms_data('summary_wm_p05'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_k
        self.report.append(self.get_iqms_data('summary_wm_k'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_n
        self.report.append(self.get_iqms_data('summary_wm_n'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # summary_bg_mean
        self.report.append(self.get_iqms_data('summary_bg_mean'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_stdv
        self.report.append(self.get_iqms_data('summary_bg_stdv'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_median
        self.report.append(self.get_iqms_data('summary_bg_median'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_mad
        self.report.append(self.get_iqms_data('summary_bg_mad'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_p95
        self.report.append(self.get_iqms_data('summary_bg_p95'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_p05
        self.report.append(self.get_iqms_data('summary_bg_p05'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_k
        self.report.append(self.get_iqms_data('summary_bg_k'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_n
        self.report.append(self.get_iqms_data('summary_bg_n'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # fwhm_x
        self.report.append(self.get_iqms_data('fwhm_x'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_y
        self.report.append(self.get_iqms_data('fwhm_y'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_z
        self.report.append(self.get_iqms_data('fwhm_z'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_avg
        self.report.append(self.get_iqms_data('fwhm_avg'))
        self.report.append(Spacer(0 * mm, 7 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("summary_csf_k", "fwhm_x"):
            self.report.append(Paragraph(
                "<font size = 8><sup>" +
                    self.dict4runtime['extra_info'][param][2] +
                    "</sup>" + self.dict4runtime['extra_info'][param][3] +
                    "</font>",
                self.styles['Left']))

        self.report.append(PageBreak())

       # Fifth page page - IQMs###############################################
        #######################################################################
        self.report.append(Paragraph("<font size = 18 ><b>Image parameters</b>"
                                        "</font>",
                                     self.styles['Center']))
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        ### Tissues Quality ###################################################
        self.report.append(Paragraph("<font size = 15 ><b>TISSUES QUALITY</b>"
                                        "</font>",
                                     self.styles['Bullet1']))
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(Paragraph(
              "Metrics that do not fall into the above categories: statistical "
                "properties of tissue distributions, volume overlap of "
                "tissues, image harpness/blurriness, etc.",
              self.styles['Left2']))
        self.report.append(Spacer(0 * mm, 10 * mm))
        # icvs_csf
        self.report.append(self.get_iqms_data('icvs_csf'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # icvs_gm
        self.report.append(self.get_iqms_data('icvs_gm'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # icvs_wm
        self.report.append(self.get_iqms_data('icvs_wm'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # rpve_csf
        self.report.append(self.get_iqms_data('rpve_csf'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # rpve_gm
        self.report.append(self.get_iqms_data('rpve_gm'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # rpve_wm
        self.report.append(self.get_iqms_data('rpve_wm'))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # tpm_overlap_csf
        self.report.append(self.get_iqms_data('tpm_overlap_csf'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # tpm_overlap_gm
        self.report.append(self.get_iqms_data('tpm_overlap_gm'))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # tpm_overlap_wm
        self.report.append(self.get_iqms_data('tpm_overlap_wm'))
        self.report.append(Spacer(0 * mm, 147 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("rpve_csf", "tpm_overlap_csf"):
            self.report.append(Paragraph(
                "<font size = 8><sup>" +
                self.dict4runtime['extra_info'][param][2] +
                "</sup>" + self.dict4runtime['extra_info'][param][3] +
                "</font>",
                self.styles['Left']))

        self.report.append(PageBreak())

        # Sixth page - slice planes display - Raw anatomic ####################
        #######################################################################
        self.report.append(Paragraph("<font size = 18 > <b>MRI axial slice "
                                        "planes display</b> </font>",
                                     self.styles['Center']))
        self.report.append(Spacer(0*mm, 4*mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(Paragraph("<font size = 14 >Raw anatomic "
                                        "images</font>",
                                     self.styles['Center']))
        self.report.append(Spacer(0*mm, 10*mm))
        self.report.append(Paragraph('<font size = 9 ><i>"Neurological" '
                                        'convention, the left side of the '
                                        'image corresponds to the left side of '
                                        'the brain.</i></font>',
                                     self.styles['Center']))
        self.report.append(Spacer(0*mm, 1*mm))
        tmpdir = tempfile.TemporaryDirectory()
        slices_image = slice_planes_plot(
                                      self.anat, self.anat_fig_rows,
                                      self.anat_fig_cols,
                                      inf_slice_start=self.anat_inf_slice_start,
                                      slices_gap=self.anat_slices_gap,
                                      cmap="Greys_r", out_dir=tmpdir.name)
        # reminder: A4 == 210mmx297mm
        slices_image = Image(slices_image,
                             width=7.4803 * inch, height=9.0551 * inch)
        slices_image.hAlign = 'CENTER'
        self.report.append(slices_image)
        self.report.append(PageBreak())

        # Seventh page - slice planes display - Normalised anatomic (MNI) #####
        #######################################################################
        self.report.append(Paragraph("<font size = 18 > <b>MRI axial slice "
                                        "planes display</b> </font>",
                                     self.styles['Center']))
        self.report.append(Spacer(0*mm, 4*mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(Paragraph("<font size = 14 >Normalised anatomic "
                                        "(MNI)</font>",
                                     self.styles['Center']))

        self.report.append(Spacer(0*mm, 10*mm))
        self.report.append(Paragraph('<font size = 9 > <i> "Neurological" '
                                        'convention, the left side of the '
                                        'image corresponds to the left side of '
                                        'the brain. </i> <br/> </font>',
                                     self.styles['Center']))
        self.report.append(Spacer(0*mm, 1*mm))
        slices_image = slice_planes_plot(
                                 self.norm_anat,
                                 self.norm_anat_fig_rows,
                                 self.norm_anat_fig_cols,
                                 inf_slice_start=self.norm_anat_inf_slice_start,
                                 slices_gap=self.norm_anat_slices_gap,
                                 cmap="Greys_r", out_dir=tmpdir.name)
        # reminder: A4 == 210mmx297mm
        slices_image = Image(slices_image, width=7.4803 * inch,
                             height=9.0551 * inch)
        slices_image.hAlign = 'CENTER'
        self.report.append(slices_image)


        # Currently, we make and save (in derived_data) the qi2 graph, but we
        # don't include it in the report because the result with the test data
        # looks strange.
        _ = plot_qi2(np.asarray(self.iqms_data['histogram_qi2_x_grid']),
                     np.asarray(self.iqms_data['histogram_qi2_ref_pdf']),
                     np.asarray(self.iqms_data['histogram_qi2_fit_pdf']),
                     np.asarray(self.iqms_data['histogram_qi2_ref_data']),
                     int(self.iqms_data['histogram_qi2_cutoff_idx']),
                     out_file=os.path.join(
                         os.path.split(os.path.split(self.anat)[0])[0],
                         'derived_data', 'qi2_plot.svg'))

        self.page.build(self.report, canvasmaker=PageNumCanvas)
        tmpdir.cleanup()

    def get_iqms_data(self, param):
        """blabla"""

        if self.iqms_data.get(param) is None:

            if self.dict4runtime['extra_info'][param][2] is None:
                return Paragraph(
                         "<font size = 11> <b> " +
                            self.dict4runtime['extra_info'][param][0] +
                            " </b> </font>: Not determined",
                         self.styles['Bullet2'])
            else:
                return Paragraph(
                        "<font size = 9><sup>" +
                            self.dict4runtime['extra_info'][param][2] +
                            "</sup></font><font size = 11><b>" +
                            self.dict4runtime['extra_info'][param][0] +
                            "</b></font>: Not determined",
                        self.styles['Bullet2'])

        else:

            if self.dict4runtime['extra_info'][param][2] is None:
                return Paragraph(
                    "<font size = 11> <b>{0}</b> </font>: {1}".format(
                    self.dict4runtime['extra_info'][param][0],
                    str(round(self.iqms_data.get(param),
                              self.dict4runtime['extra_info'][param][1])) if
                        isinstance(self.dict4runtime['extra_info'][param][1],
                                   int) else
                    self.dict4runtime['extra_info'][param][1].format(
                                                        self.iqms_data[param])),
                    self.styles['Bullet2'])

            else:
                return Paragraph(
                    "<font size = 9><sup>{0}</sup></font><font size = 11><b>"
                    "{1}</b></font>: {2}".format(
                    self.dict4runtime['extra_info'][param][2],
                    self.dict4runtime['extra_info'][param][0],
                    str(round(self.iqms_data.get(param),
                              self.dict4runtime['extra_info'][param][1])) if
                        isinstance(self.dict4runtime['extra_info'][param][1],
                                   int) else
                    self.dict4runtime['extra_info'][param][1].format(
                                                        self.iqms_data[param])),
                    self.styles['Bullet2'])


