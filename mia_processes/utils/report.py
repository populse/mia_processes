# -*- coding: utf-8 -*-
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

# Other import
import json
import os
import platform
import tempfile
import time
from datetime import datetime
from math import pi, sqrt
from sys import version

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

# capsul import
from capsul import info as capsul_info

# nipype import
from nipype import info as nipype_info

# populse_mia import
from populse_mia import info as mia_info
from populse_mia import sources_images
from reportlab.lib import colors

# reportlab import
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT, TA_RIGHT
from reportlab.lib.pagesizes import A4, portrait
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch, mm
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)
from scipy.interpolate import splev, splrep
from scipy.io import loadmat

# mia_processes import:
from mia_processes import info as mia_processes_info
from mia_processes.utils import (
    PageNumCanvas,
    ReportLine,
    plot_qi2,
    plot_segmentation,
    plot_slice_planes,
)

# import openpyxl


class Report:
    """Create pdf report

    IQMs_file --> mriqc individual report (with all IQMs)
    mriqc_group --> mriqc report group

    Methods:
      - get_iqms_data
      - co2_inhal_cvr_make_report
      - mriqc_anat_make_report
      - mriqc_func_make_report
      - mriqc_group_make_report

    """

    def __init__(self, report_file, dict4runtime, **kwargs):
        """Create Canvas , create cover and make report"""

        today_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
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
                self.title = (
                    "<font size=18><b>Anatomical Image-Quality "
                    "Metrics summary report</b></font>"
                )

            elif "func" in kwargs:
                self.make_report = self.mriqc_func_make_report
                ref_exp = self.func
                self.title = (
                    "<font size=18><b>Functional Image-Quality "
                    "Metrics summary report</b></font>"
                )

            self.header_title = (
                "<sup rise=20 size=9>$</sup>"
                "<font size=30><b>MRIQ</b></font>"
                "<font size=11>uality</font>"
                "<font size=30><b>C</b></font>"
                "<font size=11>ontrol</font>"
            )
            infos = [
                self.dict4runtime[i]
                for i in (
                    "Site",
                    "Spectro",
                    "StudyName",
                    "AcquisitionDate",
                    "PatientName",
                    "Sex",
                    "Age",
                )
            ]
            infos.insert(4, today_date)
            infos.insert(5, ref_exp)

            headers = [
                "SITE",
                "MRI SCANNER",
                "STUDY NAME",
                "EXAMINATION DATE",
                "MRIQC CALCULATION DATE",
                "NAME OF THE INPUT DATA",
                "PATIENT REFERENCE",
                "PATIENT SEX",
                "PATIENT AGE",
                "SOFTWARES",
            ]

        elif "mriqc_group" in kwargs:
            self.make_report = self.mriqc_group_make_report
            self.title = "<font size=18><b>MRIQc group report</b></font>"

            self.header_title = (
                "<sup rise=20 size=9>$</sup>"
                "<font size=30><b>MRIQ</b></font>"
                "<font size=11>uality</font>"
                "<font size=30><b>C</b></font>"
                "<font size=11>ontrol</font>"
            )
            infos = [today_date]
            headers = [
                "MRIQC GROUPS CALCULATION DATE",
                "SOFTWARES",
            ]

        elif "CVR" in kwargs:
            self.make_report = self.co2_inhal_cvr_make_report
            ref_exp = {
                "norm_anat": self.norm_anat,
                "norm_func": self.norm_func,
            }
            self.title = (
                "<font size=18><b>{0} report: </b>{1}"
                "</font>".format(
                    self.dict4runtime["norm_anat"]["PatientRef"],
                    today_date.split(" ")[0],
                )
            )

            self.header_title = (
                "<font size=30><b>C</b></font>"
                "<font size=11>erebro</font>"
                "<font size=30><b>V</b></font>"
                "<font size=11>ascular</font>"
                "<font size=30><b> R</b></font>"
                "<font size=11>reactivity <br/></font>"
                "<font size=7>with CO<sub>2</sub> inhalation "
                "as a physiological challenge</font>"
            )
            infos = [
                self.dict4runtime["norm_anat"][i]
                for i in (
                    "Site",
                    "Spectro",
                    "StudyName",
                    "AcquisitionDate",
                    "PatientRef",
                    "Sex",
                    "Age",
                    "Pathology",
                )
            ]
            infos.insert(4, today_date)
            infos.insert(5, ref_exp)
            infos.insert(10, "Undefined")

            headers = [
                "SITE",
                "MRI SCANNER",
                "STUDY NAME",
                "EXAMINATION DATE",
                "CVR CALCULATION DATE",
                "NAME OF THE INPUT DATA",
                "PATIENT REFERENCE",
                "PATIENT SEX",
                "PATIENT AGE",
                "PATHOLOGY",
                "REFERENCE GROUP",
                "SOFTWARES",
            ]

        infos.extend(
            (
                "Python " + version.split()[0],
                "Populse_mia " + mia_info.__version__,
                "Capsul " + capsul_info.__version__,
                "Nipype " + nipype_info.__version__,
                "Mia_processes " + mia_processes_info.__version__,
                "Operating System {0} {1}".format(
                    (
                        "Mac OS X"
                        if platform.system() == "Darwin"
                        else platform.system()
                    ),
                    platform.release(),
                ),
            )
        )

        headers.extend([" "] * 5)
        # Initialises stylesheet with few basic heading and text styles,
        # return a stylesheet object
        self.styles = getSampleStyleSheet()
        # ParagraphStyle gives all the attributes available for formatting
        # paragraphs
        self.styles.add(ParagraphStyle(name="Center", alignment=TA_CENTER))
        self.styles.add(ParagraphStyle(name="Center2", alignment=TA_CENTER))
        # If more than 1 line
        self.styles["Center2"].leading = 24
        self.styles.add(ParagraphStyle(name="Justify", alignment=TA_JUSTIFY))
        self.styles.add(ParagraphStyle(name="Right", alignment=TA_RIGHT))
        self.styles.add(ParagraphStyle(name="Left", alignment=TA_LEFT))
        self.styles.add(ParagraphStyle(name="Left2", alignment=TA_LEFT))
        # For left indent 30
        self.styles["Left2"].leftIndent = 30
        self.styles.add(
            ParagraphStyle(
                name="Bullet1",
                leftIndent=30,
                bulletOffsetY=2,
                bulletIndent=20,
                bulletFontSize=6,
                bulletColor="black",
                bulletText="●",
            )
        )
        self.styles.add(
            ParagraphStyle(
                name="Bullet2",
                leftIndent=60,
                bulletOffsetY=1,
                bulletIndent=50,
                bulletFontSize=6,
                bulletColor="black",
                bulletText="❍",
            )
        )
        # Logo populse
        header_image_1 = os.path.join(
            sources_images.__path__[0], "Logo_populse_square.jpg"
        )
        # 2156px x 2156px
        header_image_1 = Image(header_image_1, 43.0 * mm, 43.0 * mm)
        # Logo Mia
        header_image_2 = os.path.join(
            sources_images.__path__[0], "Logo_populse_mia_HR.jpeg"
        )
        # 2505px x 1843px
        header_image_2 = Image(header_image_2, 54.4 * mm, 40.0 * mm)
        self.image_cov = Table(
            [[header_image_1, header_image_2]], [70 * mm, 70 * mm]
        )
        self.image_cov.setStyle(
            TableStyle(
                [
                    ("ALIGN", (0, 0), (0, 0), "LEFT"),
                    ("ALIGN", (-1, 0), (-1, 0), "RIGHT"),
                    ("VALIGN", (0, 0), (-1, 0), "MIDDLE"),
                ]
            )
        )
        cover_data = [[0, 0]]

        for header, info in zip(headers, infos):
            if isinstance(info, dict):
                info = " AND ".join(str(val) for val in info.values())

            if cover_data == [[0, 0]]:
                if header == " ":
                    cover_data[0][0] = Paragraph(
                        "<para align=right>" + header + "</para>",
                        self.styles["Normal"],
                    )

                else:
                    # Class reportlab.platypus.paragraph.Paragraph
                    # with XML-like markup
                    cover_data[0][0] = Paragraph(
                        "<para align=right><b>" + header + "&nbsp&nbsp&nbsp :"
                        "</b></para>",
                        self.styles["BodyText"],
                    )

                cover_data[0][1] = Paragraph(
                    "<para align=justify>" + info + "</para>",
                    self.styles["Normal"],
                )

            else:
                if header == " ":
                    temp = [
                        Paragraph(
                            "<para align=right>" + header + "</para>",
                            self.styles["Normal"],
                        ),
                        Paragraph(
                            "<para align=justify>" + info + "</para>",
                            self.styles["Normal"],
                        ),
                    ]

                else:
                    temp = [
                        Paragraph(
                            "<para align=right><b>"
                            + header
                            + "&nbsp&nbsp&nbsp :</b></para>",
                            self.styles["BodyText"],
                        ),
                        Paragraph(
                            (
                                "<para align=justify><font size = 5>"
                                + info
                                + "</font></para>"
                                if header == "NAME OF THE INPUT DATA"
                                else "<para align=justify>"
                                + str(info)
                                + "</para>"
                            ),
                            self.styles["Normal"],
                        ),
                    ]

                cover_data.append(temp)

        # colWidths, rowHeights
        self.cover_data = Table(cover_data, [60 * mm, 97 * mm])
        # Careful: Table use (col, raw) and 0-base for top left start as usual
        # OR -1-base for lower right start
        self.cover_data.setStyle(
            TableStyle(
                [
                    ("LINEABOVE", (0, 0), (-1, 0), 2, colors.black),
                    ("LINEBELOW", (0, -1), (-1, -1), 2, colors.black),
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ]
            )
        )
        self.cover_data.hAlign = "CENTER"

        self.textDisclaimer = (
            "<font size=7 ><i><b>DISCLAIMER</b><br/>Mia "
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
            "responsabilité.</i> </font>"
        )

        # Output document template definition for page; margins and pages size
        self.page = SimpleDocTemplate(
            report_file,
            pagesize=portrait(A4),
            rightMargin=20 * mm,
            leftMargin=20 * mm,
            topMargin=10 * mm,
            bottomMargin=10 * mm,
        )

        # Canvas creation
        self.report = []

    def get_iqms_data(self, param):
        """Get iqms data"""

        if self.iqms_data.get(param) is None:
            if self.dict4runtime["extra_info"][param][2] is None:
                return Paragraph(
                    "<font size = 11> <b> "
                    + self.dict4runtime["extra_info"][param][0]
                    + " </b> </font>: Not determined",
                    self.styles["Bullet2"],
                )
            else:
                return Paragraph(
                    "<font size = 9><sup>"
                    + self.dict4runtime["extra_info"][param][2]
                    + "</sup></font><font size = 11><b>"
                    + self.dict4runtime["extra_info"][param][0]
                    + "</b></font>: Not determined",
                    self.styles["Bullet2"],
                )

        else:
            if self.dict4runtime["extra_info"][param][2] is None:
                return Paragraph(
                    "<font size = 11> <b>{0}</b> </font>: {1}".format(
                        self.dict4runtime["extra_info"][param][0],
                        (
                            str(
                                round(
                                    self.iqms_data.get(param),
                                    self.dict4runtime["extra_info"][param][1],
                                )
                            )
                            if isinstance(
                                self.dict4runtime["extra_info"][param][1], int
                            )
                            else self.dict4runtime["extra_info"][param][
                                1
                            ].format(self.iqms_data[param])
                        ),
                    ),
                    self.styles["Bullet2"],
                )

            else:
                return Paragraph(
                    "<font size = 9><sup>{0}</sup></font><font size = 11><b>"
                    "{1}</b></font>: {2}".format(
                        self.dict4runtime["extra_info"][param][2],
                        self.dict4runtime["extra_info"][param][0],
                        (
                            str(
                                round(
                                    self.iqms_data.get(param),
                                    self.dict4runtime["extra_info"][param][1],
                                )
                            )
                            if isinstance(
                                self.dict4runtime["extra_info"][param][1], int
                            )
                            else self.dict4runtime["extra_info"][param][
                                1
                            ].format(self.iqms_data[param])
                        ),
                    ),
                    self.styles["Bullet2"],
                )

    def co2_inhal_cvr_make_report(self):
        """Make CVR under CO2 challenge report"""

        # page 1 - cover ##################################################
        ###################################################################

        self.report.append(self.image_cov)
        # width, height
        self.report.append(Spacer(0 * mm, 8 * mm))
        self.report.append(Paragraph(self.header_title, self.styles["Center"]))
        self.report.append(Spacer(0 * mm, 10 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(Paragraph(self.title, self.styles["Center"]))
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(self.cover_data)
        self.report.append(Spacer(0 * mm, 4 * mm))
        self.report.append(
            Paragraph(self.textDisclaimer, self.styles["Justify"])
        )
        # self.report.append(Spacer(0 * mm, 6 * mm))
        self.report.append(PageBreak())

        # page 2 - Anatomical MRI - Acq & Post-pro parameters############
        #################################################################
        self.report.append(
            Paragraph(
                "<font size=18 ><b> Anatomical MRI <br/></b></font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                "<font size=15 ><b>Acquisition parameters:<br/></b></font>",
                self.styles["Left"],
            )
        )
        self.report.append(Spacer(0 * mm, 5 * mm))
        # TODO: We don't currently have the "Acquisition nbr" tag in the raw
        #       anat data. Would you like this to come from mri_conv?
        self.report.append(
            Paragraph(
                f"<font size=11> <b> Protocol name / Acquisition nr: </b> "
                f"</font> {self.dict4runtime['norm_anat']['ProtocolName']} / "
                f"{self.dict4runtime['norm_anat']['Acquisition nbr']}",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(
            Paragraph(
                "<font size=11> <b>Acquisition mode</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        self.report.append(
            Paragraph(
                f"<font size=11> <b> Sequence name: </b> "
                f"</font> {self.dict4runtime['norm_anat']['SequenceName']}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(
            Paragraph(
                "<font size=11> <b>Raw geometry parameters</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        d_dim = self.dict4runtime["norm_anat"][
            "Dataset dimensions (Count, X,Y,Z,T...)"
        ]
        self.report.append(
            Paragraph(
                f"<font size=11> <b> Dataset dimension:</b> "
                f"</font> {d_dim[1:]}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        sl_thick = self.dict4runtime["norm_anat"]["SliceThickness"]

        if not sl_thick == "Undefined":
            sl_thick = sl_thick[0]

        st_end_sl = self.dict4runtime["norm_anat"]["Start/end slice"]

        if not st_end_sl == "Undefined" and st_end_sl == [0, 0]:
            st_end_sl = 0.0

        self.report.append(
            Paragraph(
                f"<font size=11> <b> Slice thickness "
                f"/ Slice gap [mm]:</b> "
                f"</font> {sl_thick} / {st_end_sl}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        fov = self.dict4runtime["norm_anat"]["FOV"]

        if fov == "Undefined":
            fov = ["Undefined"] * 3

        if all(isinstance(elt, (int, float)) for elt in fov):
            fov = [round(elt, 1) for elt in fov]

        self.report.append(
            Paragraph(
                f"<font size=11> <b> FOV (ap / fh / rl) [mm]:</b> </font> "
                f"{fov[0]} / {fov[1]} / {fov[2]}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        # TODO: I don't know how to obtain the scan resolution parameter yet
        self.report.append(
            Paragraph(
                "<font size=11> <b> Scan resolution  "
                "(x / y):</b> </font> Undefined",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        vox_size = self.dict4runtime["norm_anat"][
            "Grid spacings (X,Y,Z,T,...)"
        ]

        if vox_size == "Undefined":
            vox_size = ["Undefined"] * 3

        if all(isinstance(elt, (int, float)) for elt in vox_size):
            vox_size = [round(elt, 1) for elt in vox_size]

        self.report.append(
            Paragraph(
                f"<font size=11> <b> Voxel size (x / y / z) [mm]:"
                f"</b> </font> {vox_size[0]} / {vox_size[1]} / {vox_size[2]}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(
            Paragraph(
                "<font size=11> <b>Other parameters" "</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        tr = self.dict4runtime["norm_anat"]["RepetitionTime"]
        te = self.dict4runtime["norm_anat"]["EchoTime"]
        flipang = self.dict4runtime["norm_anat"]["FlipAngle"]

        if flipang != "Undefined":
            flipang = round(flipang[0], 1)

        if te != "Undefined":
            te = round(te[0], 1)

        if tr != "Undefined":
            tr = round(tr[0], 1)

        self.report.append(
            Paragraph(
                f"<font size=11> <b> TR [ms] / TE [ms] / Image flip angle"
                f" [deg]:</b> </font> {tr} / {te} / {flipang}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))

        if len(d_dim) == 4:
            dyn_nb = 1

        elif len(d_dim) == 5:
            dyn_nb = d_dim[-1]

        self.report.append(
            Paragraph(
                f"<font size=11> <b> Number of dynamics:</b> "
                f"</font> {dyn_nb}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        # TODO: I don't know how to obtain the acquisition duration
        #       parameter yet
        self.report.append(
            Paragraph(
                "<font size=11> <b> Acquisition duration [s]:"
                "</b> </font> Undefined",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                "<font size=15 > <b> Post-processing:<br /> </b> </font>",
                self.styles["Left"],
            )
        )
        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(
            Paragraph(
                "<font size=11> <b>Spatial processing</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        # Hard-coded, we know the state of the normalization in the
        # CVR CO2 pipeline:
        self.report.append(
            Paragraph(
                "<font size=11> <b> Spatial " "normalization:</b> </font> Yes",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        brain_templ = self.dict4runtime["norm_anat"][
            "Affine regularization type"
        ]
        self.report.append(
            Paragraph(
                f"<font size=11> <b> Brain template:"
                f"</b> </font> {brain_templ}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        vox_size = self.dict4runtime["norm_anat"]["Voxel sizes"]

        if vox_size == "Undefined":
            vox_size = ["Undefined"] * 3

        if all(isinstance(elt, (int, float)) for elt in vox_size):
            vox_size = [round(elt, 1) for elt in vox_size]

        self.report.append(
            Paragraph(
                f"<font size=11> <b> Final voxel size "
                f"(x / y / z) [mm]:</b> "
                f"</font> {vox_size[0]} / {vox_size[1]} / {vox_size[2]}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        # Hard-coded, we know the state of the segmentation in the
        # CVR CO2 pipeline:
        self.report.append(
            Paragraph(
                "<font size=11> <b> Segmentation:</b>"
                " </font> Grey matter, white matter, "
                "cerebrospinal fluid, bone, soft tissues",
                self.styles["Bullet2"],
            )
        )
        self.report.append(PageBreak())

        # page 3 - Anatomical MRI - slice planes mosaic display ###########
        ###################################################################
        self.report.append(
            Paragraph(
                "<font size=15 > <b>MNI normalized axial anatomical "
                "images:</b> </font>",
                self.styles["Left"],
            )
        )
        self.report.append(Spacer(0 * mm, 20 * mm))
        self.report.append(
            Paragraph(
                '<font size=9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        tmpdir = tempfile.TemporaryDirectory()
        slices_image = plot_slice_planes(
            data_1=self.norm_anat,
            fig_rows=self.norm_anat_fig_rows,
            fig_cols=self.norm_anat_fig_cols,
            slice_start=self.norm_anat_inf_slice_start,
            slice_step=self.norm_anat_slices_gap,
            cmap_1=self.norm_anat_cmap,
            vmin_1=self.norm_anat_vmin,
            vmax_1=self.norm_anat_vmax,
            out_dir=tmpdir.name,
        )
        # remainder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)

        self.report.append(PageBreak())

        # page 4 - fmri-cvr MRI Acq & Post-pro parameters##############
        ###############################################################
        self.report.append(
            Paragraph(
                "<font size=18 ><b> fMRI under a vasoactive stimulus "
                "(hypercapnia)<br/></b></font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        self.report.append(
            Paragraph(
                "<font size=15 ><b>Acquisition parameters:<br/></b></font>",
                self.styles["Left"],
            )
        )
        self.report.append(Spacer(0 * mm, 5 * mm))
        # TODO: We don't currently have the "Acquisition nbr" tag in the raw
        #       anat data. Would you like this to come from mri_conv?
        self.report.append(
            Paragraph(
                f"<font size=11> <b> Protocol name / Acquisition nr: "
                f"</b> </font> "
                f"{self.dict4runtime['norm_func']['ProtocolName']} / "
                f"{self.dict4runtime['norm_func']['Acquisition nbr']}",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(
            Paragraph(
                "<font size=11> <b>Acquisition mode</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        self.report.append(
            Paragraph(
                f"<font size=11> <b> Sequence name: </b> </font> "
                f"{self.dict4runtime['norm_func']['SequenceName']}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(
            Paragraph(
                "<font size=11> <b>Raw geometry parameters</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        d_dim = self.dict4runtime["norm_func"][
            "Dataset dimensions (Count, X,Y,Z,T...)"
        ]
        self.report.append(
            Paragraph(
                f"<font size=11> <b> Dataset dimension:</b> "
                f"</font> {d_dim[1:-1]}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        sl_thick = self.dict4runtime["norm_func"]["SliceThickness"]

        if not sl_thick == "Undefined":
            sl_thick = sl_thick[0]

        st_end_sl = self.dict4runtime["norm_func"]["Start/end slice"]

        if not st_end_sl == "Undefined" and st_end_sl == [0, 0]:
            st_end_sl = 0.0

        self.report.append(
            Paragraph(
                f"<font size=11> <b> Slice thickness "
                f"/ Slice gap [mm]:</b> "
                f"</font> {sl_thick} / {st_end_sl}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        fov = self.dict4runtime["norm_func"]["FOV"]

        if fov == "Undefined":
            fov = ["Undefined"] * 3

        if all(isinstance(elt, (int, float)) for elt in fov):
            fov = [round(elt, 1) for elt in fov]

        self.report.append(
            Paragraph(
                f"<font size=11> <b> FOV (ap / fh / rl) [mm]:</b> </font> "
                f"{fov[0]} / {fov[1]} / {fov[2]}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        # TODO: I don't know how to obtain the scan resolution parameter yet
        self.report.append(
            Paragraph(
                "<font size=11> <b> Scan resolution  "
                "(x / y):</b> </font> Undefined",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        raw_func_vox_size = self.dict4runtime["norm_func"][
            "Grid spacings (X,Y,Z,T,...)"
        ]

        if raw_func_vox_size == "Undefined":
            raw_func_vox_size = ["Undefined"] * 3

        if all(isinstance(elt, (int, float)) for elt in raw_func_vox_size):
            raw_func_vox_size = [round(elt, 1) for elt in raw_func_vox_size]

        self.report.append(
            Paragraph(
                f"<font size=11> <b> Voxel size (x / y / z) [mm]:</b> "
                f"</font> {raw_func_vox_size[0]} "
                f"/ {raw_func_vox_size[1]} "
                f"/ {raw_func_vox_size[2]}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(
            Paragraph(
                "<font size=11> <b>Other parameters" "</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        tr = self.dict4runtime["norm_func"]["RepetitionTime"]
        te = self.dict4runtime["norm_func"]["EchoTime"]
        flipang = self.dict4runtime["norm_func"]["FlipAngle"]

        if flipang != "Undefined":
            flipang = round(flipang[0], 1)

        if te != "Undefined":
            te = round(te[0], 1)

        if tr != "Undefined":
            tr = round(tr[0], 1)

        self.report.append(
            Paragraph(
                f"<font size=11> <b> TR [ms] / TE [ms] / Image flip angle"
                f" [deg]:</b> </font> {tr} / {te} / {flipang}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))

        if len(d_dim) == 4:
            dyn_nb = 1

        elif len(d_dim) == 5:
            dyn_nb = d_dim[-1]

        self.report.append(
            Paragraph(
                f"<font size=11> <b> Number of dynamics:</b> "
                f"</font> {dyn_nb}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        # TODO: I don't know how to obtain the acquisition duration
        #       parameter yet
        self.report.append(
            Paragraph(
                "<font size=11> <b> Acquisition duration [s]:"
                "</b> </font> Undefined",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        # TODO: I don't know how to obtain the CVR regressor (Individual or
        #       Standard) parameter yet
        self.report.append(
            Paragraph(
                "<font size=11> <b>Stimulus and performance</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))  # (width, height)
        self.report.append(
            Paragraph(
                "<font size=11> <b>CVR regressor:</b> </font>"
                + "<font size = 9>Undefined</font>",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        self.report.append(
            Paragraph(
                f"<font size=11> <b>Stimulation device / Gas name:"
                f"</b> </font> {self.dict4runtime['norm_func']['GasAdmin']} "
                f"/ {self.dict4runtime['norm_func']['Gas']}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        # TODO: Stimulus is currently hard-coded at 8%. Should we take the
        #       value from the patient_info dictionary?
        self.report.append(
            Paragraph(
                "<font size=11> <b>Stimulus:</b> </font> 8% CO<sub>2</sub>",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                "<font size=15 > <b> Post-processing:<br /> </b> </font>",
                self.styles["Left"],
            )
        )
        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(
            Paragraph(
                "<font size=11> <b>Spatial processing</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        # Hard-coded, we know the state of the normalization in the
        # CVR CO2 pipeline:
        self.report.append(
            Paragraph(
                "<font size=11> <b> Spatial " "normalization:</b> </font> Yes",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        brain_templ = self.dict4runtime["norm_func"][
            "Affine regularization type"
        ]
        self.report.append(
            Paragraph(
                f"<font size=11> <b> Brain template:"
                f"</b> </font> {brain_templ}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        vox_size = self.dict4runtime["norm_func"]["Voxel sizes"]

        if vox_size == "Undefined":
            vox_size = ["Undefined"] * 3

        if all(isinstance(elt, (int, float)) for elt in vox_size):
            vox_size = [round(elt, 1) for elt in vox_size]

        self.report.append(
            Paragraph(
                f"<font size=11> <b> Voxel size after normalization "
                f"(x / z / y) [mm]:</b> "
                f"</font> {vox_size[0]} / {vox_size[1]} / {vox_size[2]}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        fwhm = self.dict4runtime["smooth_norm_func"][
            "FWHM (X, Y, Z) for Smooth"
        ]

        if fwhm == "Undefined":
            fwhm = ["Undefined"] * 3

        if all(isinstance(elt, (int, float)) for elt in fwhm):
            fwhm = [round(elt, 1) for elt in fwhm]

        self.report.append(
            Paragraph(
                f"<font size=11> <b>Spatial smoothing "
                f"(fwhm) [mm]:</b> </font> {fwhm[0]} "
                f"/ {fwhm[1]} / {fwhm[2]}",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        self.report.append(
            Paragraph(
                "<font size=11> <b>Head motion as "
                "nuisance regressor in GLM:"
                "</b> </font> Yes",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        # TODO: The pipeline doesn't currently use ART. Should we implement it?
        self.report.append(
            Paragraph(
                "<font size=11> <b>Artifact detection "
                "tools (ART):</b> </font> No",
                self.styles["Bullet2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        # Hard-coded, we know the state of the mask for the func in the
        # CVR CO2 pipeline:
        self.report.append(
            Paragraph(
                "<font size=11> <b> Mask:</b> </font> " "Grey matter",
                self.styles["Bullet2"],
            )
        )
        self.report.append(PageBreak())

        # page 5 - fmri-cvr MRI quality check: movements & EtCO2 regressor ####
        #######################################################################
        sources_images_dir = os.path.join(
            os.path.dirname(os.path.dirname(__file__)), "sources_images"
        )
        (
            rmsdTra,
            maxRmsdTra,
            maxTra,
            maxTraCheck,
            minTra,
            rmsdRot,
            maxRmsdRot,
            maxRot,
            maxRotCheck,
            minRot,
        ) = np.zeros(10)

        try:
            data_rp = np.loadtxt(self.realignment_parameters)
            data_rp[:, 3:] = data_rp[:, 3:] * 180 / pi  # Rad to deg conversion

        except IOError:
            print(
                "\n==> The "
                + self.realignment_parameters
                + " was not found ! <==\n"
            )
            data_rp = None

        try:
            matDataReg = loadmat(self.regressor_physio)

        except IOError:
            print(
                "\n==> The " + self.regressor_physio + " was not found ! <==\n"
            )
            matDataReg = None

        try:
            # TODO: We state that the file to take is "s" + self.norm_func.
            #       1- Since we want the BOLD time course in bold, perhaps
            #          we could take only the self.norm_func?
            #       2- Should we simply add "s" to self.norm_fun or add the
            #          data as an input to the brick?
            folder, file = os.path.split(self.norm_func)
            snorm_func = os.path.join(folder, "s" + file)
            brain_img_snorm_func = nib.load(snorm_func)
            snorm_func_data = brain_img_snorm_func.get_fdata()
            # TODO: We make the mask name from the self.norm_anat.
            #       Should we simply add the data as an input to the brick?
            folder, file = os.path.split(self.norm_anat)
            file, ext = os.path.splitext(file)
            grey_mat_mask = os.path.join(
                folder, "mask_swc1" + file[1:] + "_003" + ext
            )
            brain_img_grey_mat_mask = nib.load(grey_mat_mask)
            mask_data = brain_img_grey_mat_mask.get_fdata()
            binary_mask = np.where(mask_data > 0.1, 1, 0)
            masked_data = snorm_func_data * binary_mask[..., np.newaxis]
            bold_signal_tc = np.mean(masked_data, axis=(0, 1, 2))

        except Exception:
            print("\n==> The " + grey_mat_mask + " was not found ! <==\n")
            bold_signal_tc = None

        im_qualCheck = None
        im_qualCheckTra = None
        im_qualCheckRot = None
        im_qualCheckReg = None
        im_qualCheckBoldTC = None
        # figsize in inches
        fig = plt.figure(figsize=(12, 5.6), facecolor="white")
        # fig.set_size_inches(254 / 25.4, 142 / 25.4)
        ax = fig.add_subplot(111)
        fig.subplots_adjust(
            left=None,
            bottom=None,
            right=None,
            top=None,
            wspace=None,
            hspace=None,
        )
        xLim = None

        if data_rp is not None:
            xLim = len(data_rp)

            for i in range(3):

                for n in data_rp[:, i]:
                    rmsdTra = rmsdTra + n**2

                rmsdTra = sqrt(rmsdTra / len(data_rp[:, i]))

                if rmsdTra > maxRmsdTra:
                    # maxRmsdTra: Maximum root-mean-square deviation for
                    # x, y, z linear motion
                    maxRmsdTra = rmsdTra

                rmsdTra = 0

                if max(data_rp[:, i]) > maxTra:
                    # maxTra: Maximum x, y, z linear motion in
                    # POSITIVE direction
                    maxTra = max(data_rp[:, i])

                if min(data_rp[:, i]) < minTra:
                    # minTra: Maximum x, y, z linear motion in
                    # NEGATIVE direction
                    minTra = min(data_rp[:, i])

            maxTraCheck = maxTra

            if abs(minTra) > maxTra:
                maxTraCheck = abs(minTra)

            for i in range(3, 6):

                for n in data_rp[:, i]:
                    rmsdRot = rmsdRot + n**2

                rmsdRot = sqrt(rmsdRot / len(data_rp[:, i]))

                if rmsdRot > maxRmsdRot:
                    # maxRmsdRot: Maximum root-mean-square deviation for
                    # roll, pitch, yaw rotational motion
                    maxRmsdRot = rmsdRot

                rmsdRot = 0

                if max(data_rp[:, i]) > maxRot:
                    # maxRot: Maximum roll, pitch, yaw rotational motion in
                    # POSITIVE direction
                    maxRot = max(data_rp[:, i])

                if min(data_rp[:, i]) < minRot:
                    # minRot: Maximum roll, pitch, yaw rotational motion in
                    # NEGATIVE direction
                    minRot = min(data_rp[:, i])

            maxRotCheck = maxRot

            if abs(minRot) > maxRot:
                maxRotCheck = abs(minRot)

            if "Undefined" not in raw_func_vox_size:
                # Mean x,y,z voxel dimension calculation.
                traLim = (
                    (
                        raw_func_vox_size[0]
                        + raw_func_vox_size[1]
                        + raw_func_vox_size[2]
                    )
                    * 1
                    / 3
                )

                # maxRot for quality fixed to 2deg
                # maxTra(nslation) for quality fixed to
                # mean x,y,z voxel dim (traLim)
                if (maxTraCheck >= traLim) or (maxRotCheck >= 2):
                    # 912px × 892px
                    im_qualCheck = Image(
                        os.path.join(sources_images_dir, "No-check-mark.png"),
                        13.0 * mm,
                        12.7 * mm,
                    )

                else:
                    # 940px × 893px
                    im_qualCheck = Image(
                        os.path.join(sources_images_dir, "OK-check-mark.png"),
                        13.0 * mm,
                        12.4 * mm,
                    )

            ax.set_title("Linear head motion parameters", fontsize=20, y=1.03)
            ax.set_xlabel("Dynamic scans", fontsize=14)
            ax.set_ylabel("Linear motion -X, Y, Z- (mm)", fontsize=14)
            ax.plot(data_rp[:, :3], label=("X", "Y", "Z"))
            ax.legend(loc="best")
            ax.set_xlim(0, xLim)
            ax.set_ylim(minTra * 1.1, maxTra * 1.1)
            ax.yaxis.grid(
                True, linestyle="-", which="major", color="grey", alpha=0.5
            )
            ax.xaxis.grid(
                True, linestyle="-", which="major", color="grey", alpha=0.5
            )
            ax.get_yaxis().set_tick_params(direction="out")
            ax.get_xaxis().set_tick_params(direction="out")

            out_file_tra = os.path.join(
                tmpdir.name,
                self.dict4runtime["norm_anat"]["PatientRef"]
                + "_CVR_QualityControlMeasure_translation.png",
            )
            # High resolution: 2000px × 1118px.
            fig.savefig(out_file_tra, format="png", dpi=200)
            # im_qualCheckTra = Image(out_file_tra, 153.91 * mm, 86 * mm)
            im_qualCheckTra = Image(out_file_tra, 7.187 * inch, 3.354 * inch)
            ax.clear()
            ax.set_title(
                "Rotational head motion parameters", fontsize=20, y=1.03
            )
            ax.set_xlabel("Dynamic scans", fontsize=14)
            ax.set_ylabel(
                "Rotational motion -Roll, Pitch, Yaw- (°)", fontsize=14
            )
            ax.plot(data_rp[:, 3:], label=("Roll", "Pitch", "Yaw"))
            ax.legend(loc="best")
            ax.set_xlim(0, xLim)
            ax.set_ylim(minRot * 1.1, maxRot * 1.1)
            ax.yaxis.grid(
                True, linestyle="-", which="major", color="grey", alpha=0.5
            )
            ax.xaxis.grid(
                True, linestyle="-", which="major", color="grey", alpha=0.5
            )
            out_file_rot = os.path.join(
                tmpdir.name,
                self.dict4runtime["norm_anat"]["PatientRef"]
                + "_CVR_QualityControlMeasure_Rotation.png",
            )
            # High resolution: 2000px × 1118px.
            fig.savefig(out_file_rot, format="png", dpi=200)
            im_qualCheckRot = Image(out_file_rot, 7.187 * inch, 3.354 * inch)

        if matDataReg is not None:

            if xLim is None:
                xLim = len(matDataReg["R"])

            ax.clear()
            ax.set_title(
                r"$\mathsf{EtCO_2}$ variation regressor (Standard)",
                fontsize=20,
                y=1.03,
            )
            ax.set_xlabel("Dynamic scans", fontsize=14)
            ax.set_ylabel(r"$\mathsf{\Delta EtCO_2}$ ( mmHg)", fontsize=14)
            ax.plot(matDataReg["R"], label=r"$\mathsf{EtCO_2}$")
            # ax.legend(loc='best')
            ax.set_xlim(0, xLim)
            ax.set_ylim(
                matDataReg["R"].min() * 1.1, matDataReg["R"].max() * +1.1
            )
            ax.yaxis.grid(
                True, linestyle="-", which="major", color="grey", alpha=0.5
            )
            ax.xaxis.grid(
                True, linestyle="-", which="major", color="grey", alpha=0.5
            )
            out_file_reg = os.path.join(
                tmpdir.name,
                self.dict4runtime["norm_anat"]["PatientRef"]
                + "_CVR_QualityControlMeasure_EtCO2.png",
            )
            # High resolution: 2000px × 1118px.
            fig.savefig(out_file_reg, format="png", dpi=200)
            im_qualCheckReg = Image(out_file_reg, 7.187 * inch, 3.354 * inch)

        if bold_signal_tc is not None:

            if xLim is None:
                xLim = brain_img_snorm_func.shape[3]

            ax.clear()
            ax.set_title(
                "BOLD signal timecourse (smoothed)", fontsize=20, y=1.03
            )
            ax.set_xlabel("Dynamic scans", fontsize=14)
            ax.set_ylabel("BOLD response (A. U.)", fontsize=14)
            # With smoothing
            x_sparse = np.linspace(
                0, len(bold_signal_tc) - 1, len(bold_signal_tc)
            )
            x_dense = np.linspace(
                0, len(bold_signal_tc) - 1, len(bold_signal_tc) * 10
            )

            ix = np.mean(bold_signal_tc)
            #  Lagrange polynomial interpolation coefficient determined using
            #  the experimental average bold_signal_tc (ix) and the best s (y)
            #  giving a good visual result when smoothing the bold_signal_tc.
            s = (
                (1.64e-24) * ix**5
                - (6.729e-17) * ix**4
                + (2.49e-10) * ix**3
                + (6.301e-05) * ix**2
                + (0.03408) * ix
                - 1.249
            )
            tck = splrep(x_sparse, bold_signal_tc, s=s)
            ax.plot(x_dense, splev(x_dense, tck, der=0), label="BOLD response")
            # ax.legend(loc="best")
            ax.set_xlim(0, xLim)
            ax.set_ylim(
                bold_signal_tc.min() * 0.997, bold_signal_tc.max() * +1.003
            )
            ax.yaxis.grid(
                True, linestyle="-", which="major", color="grey", alpha=0.5
            )
            ax.xaxis.grid(
                True, linestyle="-", which="major", color="grey", alpha=0.5
            )
            out_file_rot = os.path.join(
                tmpdir.name,
                self.dict4runtime["norm_anat"]["PatientRef"]
                + "_CVR_QualityControlMeasure_BOLDtimeCourse.png",
            )
            # High resolution: 2000px × 1118px.
            fig.savefig(out_file_rot, format="png", dpi=200)
            im_qualCheckBoldTC = Image(
                out_file_rot, 7.187 * inch, 3.354 * inch
            )

        qualCheckMess = Paragraph(
            "<font size=14 > <b> fMRI quality "
            "check: movements </b> </font>",
            self.styles["Left"],
        )

        if im_qualCheck is None:
            im_qualCheck = Paragraph(
                "<font size=8 > Automatic evaluation not available </font>",
                self.styles["Center"],
            )

        im_qualCheck.hAlign = "CENTER"
        title = [[qualCheckMess, im_qualCheck]]
        t = Table(title, [110 * mm, 30 * mm])  # colWidths, rowHeight
        t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE")]))
        t.hAlign = "LEFT"
        self.report.append(t)

        if im_qualCheckTra is None:
            im_qualCheckTra = Paragraph(
                "<font size=14 > Linear head motion parameters not "
                "available </font>",
                self.styles["Center"],
            )
            self.report.append(Spacer(0 * mm, 45 * mm))

        im_qualCheckTra.hAlign = "CENTER"
        self.report.append(im_qualCheckTra)

        if im_qualCheckRot is None:
            im_qualCheckRot = Paragraph(
                "<font size=14 > Rotational head motion parameters "
                "not available </font>",
                self.styles["Center"],
            )
            self.report.append(Spacer(0 * mm, 45 * mm))

        im_qualCheckRot.hAlign = "CENTER"
        self.report.append(im_qualCheckRot)

        if data_rp is None:
            self.report.append(Spacer(0 * mm, 45 * mm))

        if im_qualCheckReg is None:
            im_qualCheckReg = Paragraph(
                "<font size=14 > EtCO<sub>2</sub> variation regressor "
                "parameters not available </font>",
                self.styles["Center"],
            )

            if data_rp is not None:
                self.report.append(Spacer(0 * mm, 35 * mm))

        im_qualCheckReg.hAlign = "CENTER"
        self.report.append(im_qualCheckReg)
        self.report.append(PageBreak())

        # page 6 - fmri-cvr MRI quality check: BOLD signal timecourse Vs EtCO2
        ######################################################################
        self.report.append(
            Paragraph(
                "<font size=14 > <b> fMRI quality "
                "check: BOLD signal timecourse Vs "
                "EtCO<sub>2</sub> model </b> </font>",
                self.styles["Left"],
            )
        )
        self.report.append(Spacer(0 * mm, 30 * mm))

        if im_qualCheckBoldTC is None:
            im_qualCheckBoldTC = Paragraph(
                "<font size=14 > BOLD signal time course not "
                "available </font>",
                self.styles["Center"],
            )
            self.report.append(Spacer(0 * mm, 20 * mm))

        im_qualCheckBoldTC.hAlign = "CENTER"
        self.report.append(im_qualCheckBoldTC)
        self.report.append(Spacer(0 * mm, 30 * mm))

        if matDataReg is None:
            self.report.append(Spacer(0 * mm, 20 * mm))

        self.report.append(im_qualCheckReg)
        self.report.append(PageBreak())

        # page 7 - BOLD: MNI normalized axial images, 1st dyn ################
        ######################################################################
        self.report.append(
            Paragraph(
                "<font size=14 > <b> BOLD: MNI "
                "normalized axial images (1<sup>st</sup> "
                "dynamic)</b> </font>",
                self.styles["Left"],
            )
        )
        self.report.append(Spacer(0 * mm, 20 * mm))
        self.report.append(
            Paragraph(
                '<font size=9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        slices_image = plot_slice_planes(
            data_1=self.norm_func,
            fig_rows=self.norm_func_fig_rows,
            fig_cols=self.norm_func_fig_cols,
            slice_start=self.norm_func_inf_slice_start,
            slice_step=self.norm_func_slices_gap,
            cmap_1=self.norm_func_cmap,
            vmin_1=self.norm_func_vmin,
            vmax_1=self.norm_func_vmax,
            out_dir=tmpdir.name,
        )
        # remainder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)
        self.report.append(PageBreak())

        # page 8 - parametric maps: beta / EtCO2 #############################
        ######################################################################
        self.report.append(
            Paragraph(
                "<font size=14><b> Parametric maps: </b></font>"
                "<font name='Symbol' size=14><b><greek>b</greek></b></font>"
                "<font size=14><b> weight values in </b></font>"
                "<font name='Symbol' size=14><b><greek>D</greek></b></font>"
                "<font size=14><b>(%BOLD) / EtCO<sub>2</sub> "
                "(mmHg)</b></font>",
                self.styles["Left"],
            )
        )

        self.report.append(Spacer(0 * mm, 20 * mm))
        self.report.append(
            Paragraph(
                '<font size=9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        slices_image = plot_slice_planes(
            data_1=self.norm_anat,
            data_2=self.beta_image,
            fig_rows=self.norm_anat_fig_rows,
            fig_cols=self.norm_anat_fig_cols,
            slice_start=self.norm_anat_inf_slice_start,
            slice_step=self.norm_anat_slices_gap,
            cmap_1=self.norm_anat_cmap,
            vmin_1=self.norm_anat_vmin,
            vmax_1=self.norm_anat_vmax,
            cmap_2=self.beta_cmap,
            vmin_2=self.beta_vmin,
            vmax_2=self.beta_vmax,
            out_dir=tmpdir.name,
        )

        # remainder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)
        self.report.append(PageBreak())

        # page 9 - parametric maps: SPMt #####################################
        ######################################################################
        self.report.append(
            Paragraph(
                "<font size=14><b> Parametric maps: statistic parametric "
                "t-Map (SPMt) </b></font>",
                self.styles["Left"],
            )
        )
        self.report.append(Spacer(0 * mm, 20 * mm))
        self.report.append(
            Paragraph(
                '<font size=9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        slices_image = plot_slice_planes(
            data_1=self.norm_anat,
            data_2=self.spmT_image,
            fig_rows=self.norm_anat_fig_rows,
            fig_cols=self.norm_anat_fig_cols,
            slice_start=self.norm_anat_inf_slice_start,
            slice_step=self.norm_anat_slices_gap,
            cmap_1=self.norm_anat_cmap,
            vmin_1=self.norm_anat_vmin,
            vmax_1=self.norm_anat_vmax,
            cmap_2=self.spmT_cmap,
            vmin_2=self.spmT_vmin,
            vmax_2=self.spmT_vmax,
            out_dir=tmpdir.name,
        )

        # remainder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)

        self.report.append(PageBreak())

        # page 10 - Laterality indexes #######################################
        ######################################################################
        self.report.append(
            Paragraph(
                f"<font size=16> <b> Laterality index values for "
                f"{self.dict4runtime['norm_anat']['PatientRef']} "
                f"against the box and whisker plot for a reference "
                f"population <br/> </b> </font>",
                self.styles["Center2"],
            )
        )
        # res_anal_data = os.path.join(
        #    self.output_directory,
        #    self.dict4runtime["norm_anat"]["PatientName"] + "_data",
        #    "results_aggregation",
        #    "BOLD_IL_mean_beta.xls",
        # )
        # CVR_ref_data = self.dict4runtime["CVR_ref_data"]

        self.report.append(Spacer(0 * mm, 0 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)

        self.report.append(PageBreak())

        # page 11 - Vascular territories used for IL #########################
        ######################################################################
        self.report.append(
            Paragraph(
                "<font size=14><b> Vascular territories of the cerebral "
                "arteries used for IL determination </b></font>",
                self.styles["Left"],
            )
        )
        self.report.append(Spacer(0 * mm, 15 * mm))
        self.report.append(
            Paragraph(
                '<font size=9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        conv_roi_dir = os.path.join(
            self.output_directory,
            self.dict4runtime["norm_anat"]["PatientName"] + "_data",
            "ROI_data",
            "convROI_BOLD2",
        )
        arter_terr = ["convACA", "convACM", "convACP", "convPICA", "convSCA"]
        arter_terr_files = [
            os.path.join(conv_roi_dir, f"{a_t}_{l_r}_2.nii")
            for a_t in arter_terr
            for l_r in ["L", "R"]
        ]

        max_timeout = 60  # Max timeout in seconds
        start_time = time.time()

        # TODO: This is a 2-cts hack in case arter_terr_files doesn't already
        #       exist. We can also add this data as input. In this case,
        #       the brick will wait until the files exist.
        #       This would be cleaner
        for elmt in arter_terr_files:
            i = 0

            while not os.path.exists(elmt):

                if time.time() - start_time > max_timeout:
                    print(
                        f"CVR make report: Max timeout "
                        f"reached ({max_timeout}s). The "
                        f"file does not exist."
                    )
                    break

                print(
                    f"{elmt} file does not exist yet. "
                    f"Waiting...{i} {time.time() - start_time}s"
                )
                time.sleep(1)
                i += 1

        slices_image = plot_slice_planes(
            data_1=self.norm_anat,
            data_2=arter_terr_files,
            fig_rows=self.norm_anat_fig_rows,
            fig_cols=self.norm_anat_fig_cols,
            slice_start=self.norm_anat_inf_slice_start,
            slice_step=self.norm_anat_slices_gap,
            cmap_1=self.norm_anat_cmap,
            vmin_1=self.norm_anat_vmin,
            vmax_1=self.norm_anat_vmax,
            cmap_2=[
                (255, 0, 0),  # red
                (0, 255, 0),  # lime (green)
                (0, 0, 255),  # blue
                (210, 105, 30),  # chocolate
                (255, 255, 0),  # yellow
            ],
            vmin_2=None,
            vmax_2=None,
            out_dir=tmpdir.name,
            out_name="arterialTerritories_axial",
        )

        # remainder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)
        footnote_artery = [
            [
                Paragraph(
                    "<b> "
                    '<font size = 12 color = "rgb(0, 64, 255)"> '
                    "Posterior cerebral artery </font> "
                    '<font size = 12 color = "grey"> - </font>'
                    '<font size = 12 color = "rgb(255, 64, 64)"> '
                    "Anterior cerebral artery </font> "
                    '<font size = 12 color = "grey"> - </font>'
                    '<font size = 12 color = "rgb(0, 255, 64)"> '
                    "Middle cerebral artery </font> <br/>"
                    '<font size = 12 color = "rgb(210, 105, 30)"> '
                    "Posterior inferior cerebellar artery </font> "
                    '<font size = 12 color = "grey"> - </font>'
                    '<font size = 12 color = "rgb(228, 228, 28)"> '
                    "Superior cerebellar artery </font> "
                    "</b>",
                    self.styles["Center"],
                )
            ]
        ]
        data = Table(footnote_artery, [7.4803 * inch])
        data.setStyle(
            TableStyle([("BACKGROUND", (0, 0), (-1, -1), colors.black)])
        )
        data.hAlign = "CENTER"
        self.report.append(data)

        self.report.append(PageBreak())

        # page 12 -  Lobes of the brain used for IL ##########################
        ######################################################################
        self.report.append(
            Paragraph(
                "<font size=16><b> Lobes of the brain used for IL "
                "determination </b></font>",
                self.styles["Left"],
            )
        )
        self.report.append(Spacer(0 * mm, 15 * mm))
        self.report.append(
            Paragraph(
                '<font size=9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        brain_lobes = [
            "convROI-CING",
            "convROI-FRON",
            "convROI-INSULA",
            "convROI-OCC",
            "convROI-PAR",
            "convROI-STR",
            "convROI-TEMP",
            "convROI-THA",
        ]
        brain_lobes_files = [
            os.path.join(conv_roi_dir, f"{b_l}_{l_r}_2.nii")
            for b_l in brain_lobes
            for l_r in ["L", "R"]
        ]
        slices_image = plot_slice_planes(
            data_1=self.norm_anat,
            data_2=brain_lobes_files,
            fig_rows=self.norm_anat_fig_rows,
            fig_cols=self.norm_anat_fig_cols,
            slice_start=self.norm_anat_inf_slice_start,
            slice_step=self.norm_anat_slices_gap,
            cmap_1=self.norm_anat_cmap,
            vmin_1=self.norm_anat_vmin,
            vmax_1=self.norm_anat_vmax,
            cmap_2=[
                (0, 255, 255),  # cyan (aqua)
                (255, 255, 0),  # yellow
                (210, 105, 30),  # chocolate
                (0, 0, 255),  # blue
                (255, 0, 0),  # red
                (255, 228, 181),  # moccasin
                (0, 255, 0),  # lime (green)
                (255, 0, 255),  # magenta (fuschia)
            ],
            # cmap_2=[
            #     (0, 255, 255),  # cyan (aqua)
            #     (255, 255, 0),  # yellow
            #     # (210, 105, 30),   # chocolate
            #     (128, 0, 128),  # purple
            #     (0, 0, 255),  # blue
            #     (255, 0, 0),  # red
            #     # (255, 228, 181),  # moccasin
            #     (255, 140, 0),  # orange
            #     (0, 255, 0),  # lime (green)
            #     (255, 0, 255),  # magenta (fuschia)
            # ],
            vmin_2=None,
            vmax_2=None,
            out_dir=tmpdir.name,
            out_name="roisTerritories_axial",
        )
        # remainder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)
        footnote_artery = [
            [
                Paragraph(
                    "<b> "
                    '<font size = 12 color = "rgb(0, 64, 255)"> '
                    "Occipital lobe </font> "
                    '<font size = 12 color = "grey"> - </font>'
                    '<font size = 12 color = "rgb(255, 64, 64)"> '
                    "Parietal lobe </font> "
                    '<font size = 12 color = "grey"> - </font>'
                    '<font size = 12 color = "rgb(0, 255, 64)"> '
                    "Temporal lobe </font> "
                    '<font size = 12 color = "grey"> - </font>'
                    '<font size = 12 color = "rgb(210, 105, 30)"> '
                    "Insular lobe </font> <br/>"
                    '<font size = 12 color = "rgb(228, 228, 28)"> '
                    "Frontal lobe </font> "
                    '<font size = 12 color = "grey"> - </font>'
                    '<font size = 12 color = "rgb(0, 255, 255)"> '
                    "Cingulate cortex </font> "
                    '<font size = 12 color = "grey"> - </font>'
                    '<font size = 12 color = "rgb(255, 64, 255)"> '
                    "Thalamus </font> "
                    '<font size = 12 color = "grey"> - </font>'
                    '<font size = 12 color = "rgb(255, 228, 181)"> '
                    "Striatum </font> "
                    "</b>",
                    self.styles["Center"],
                )
            ]
        ]
        data = Table(footnote_artery, [7.4803 * inch])
        data.setStyle(
            TableStyle([("BACKGROUND", (0, 0), (-1, -1), colors.black)])
        )
        data.hAlign = "CENTER"
        self.report.append(data)

        self.page.build(self.report, canvasmaker=PageNumCanvas)
        tmpdir.cleanup()

    def mriqc_anat_make_report(
        self,
    ):
        """Make mriqc anat individual report"""

        # First page - cover ##################################################
        #######################################################################

        self.report.append(self.image_cov)
        # width, height
        self.report.append(Spacer(0 * mm, 8 * mm))
        self.report.append(Paragraph(self.header_title, self.styles["Center"]))
        self.report.append(Spacer(0 * mm, 10 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(Paragraph(self.title, self.styles["Center"]))
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(self.cover_data)
        self.report.append(Spacer(0 * mm, 6 * mm))
        self.report.append(
            Paragraph(self.textDisclaimer, self.styles["Justify"])
        )
        self.report.append(Spacer(0 * mm, 6 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        self.report.append(
            Paragraph(
                "<font size = 8><sup>$</sup>Esteban O et "
                "al., <i>MRIQC: Advancing the Automatic "
                "Prediction of Image Quality in MRI from"
                " Unseen Sites</i>, PLOS ONE 12(9)"
                ":e0184661.</font>",
                self.styles["Left"],
            )
        )
        self.report.append(PageBreak())

        # Second page - IQMs ##################################################
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 ><b>Image parameters </b></font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        # Spatial resolution ################################################
        self.report.append(
            Paragraph(
                "<font size = 15 > <b>SPATIAL RESOLUTION</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(Paragraph("Length - Spacing", self.styles["Left2"]))
        self.report.append(Spacer(0 * mm, 10 * mm))

        # size_x
        self.report.append(self.get_iqms_data("size_x"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # size_y
        self.report.append(self.get_iqms_data("size_y"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # size_z
        self.report.append(self.get_iqms_data("size_z"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # spacing_x
        self.report.append(self.get_iqms_data("spacing_x"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # spacing_y
        self.report.append(self.get_iqms_data("spacing_y"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # spacing_z
        self.report.append(self.get_iqms_data("spacing_z"))
        self.report.append(Spacer(0 * mm, 20 * mm))
        # Signal-to-noise ###################################################
        self.report.append(
            Paragraph(
                "<font size = 15 > <b>NOISE</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(
            Paragraph(
                "Impact of noise and/or evaluation of "
                "the fitness of a noise model",
                self.styles["Left2"],
            )
        )
        self.report.append(Spacer(0 * mm, 10 * mm))
        # snr_csf
        self.report.append(self.get_iqms_data("snr_csf"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # snr_wm
        self.report.append(self.get_iqms_data("snr_wm"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # snr_gm
        self.report.append(self.get_iqms_data("snr_gm"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # snr_total
        self.report.append(self.get_iqms_data("snr_total"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # snrd_csf
        self.report.append(self.get_iqms_data("snrd_csf"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # snrd_wm
        self.report.append(self.get_iqms_data("snrd_wm"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # snrd_gm
        self.report.append(self.get_iqms_data("snrd_gm"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # snrd_total
        self.report.append(self.get_iqms_data("snrd_total"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # cnr
        self.report.append(self.get_iqms_data("cnr"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # qi_2
        self.report.append(self.get_iqms_data("qi_2"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # cjv
        self.report.append(self.get_iqms_data("cjv"))
        self.report.append(Spacer(0 * mm, 52 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("snrd_csf", "cnr", "qi_2", "cjv"):
            self.report.append(
                Paragraph(
                    "<font size = 8><sup>"
                    + self.dict4runtime["extra_info"][param][2]
                    + "</sup>"
                    + self.dict4runtime["extra_info"][param][3]
                    + "</font>",
                    self.styles["Left"],
                )
            )

        self.report.append(PageBreak())

        # Third page - IQMs ###################################################
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 ><b>Image parameters </b></font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        # Spatial distribution ##############################################
        self.report.append(
            Paragraph(
                "<font size = 15 ><b>SPATIAL DISTRIBUTION </b></font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(
            Paragraph(
                "Information theory to evaluate the "
                "spatial distribution of information",
                self.styles["Left2"],
            )
        )
        self.report.append(Spacer(0 * mm, 10 * mm))
        # fber
        self.report.append(self.get_iqms_data("fber"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # efc
        self.report.append(self.get_iqms_data("efc"))
        self.report.append(Spacer(0 * mm, 20 * mm))
        # Artifacts #########################################################
        self.report.append(
            Paragraph(
                "<font size = 15 ><b>ARTIFACTS </b></font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(
            Paragraph(
                "Estimates artefacts and signal leakage due to "
                "rapid movement (e.g. eye movement or blood "
                "vessel pulsation, etc.)",
                self.styles["Left2"],
            )
        )

        self.report.append(Spacer(0 * mm, 10 * mm))
        # wm2max
        self.report.append(self.get_iqms_data("wm2max"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # qi_1
        self.report.append(self.get_iqms_data("qi_1"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # inu_range
        self.report.append(self.get_iqms_data("inu_range"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # inu_med
        self.report.append(self.get_iqms_data("inu_med"))
        self.report.append(Spacer(0 * mm, 100 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("fber", "efc", "wm2max", "qi_1", "inu_range"):
            self.report.append(
                Paragraph(
                    "<font size = 8><sup>"
                    + self.dict4runtime["extra_info"][param][2]
                    + "</sup>"
                    + self.dict4runtime["extra_info"][param][3]
                    + "</font>",
                    self.styles["Left"],
                )
            )

        self.report.append(PageBreak())

        # fourth page - IQMs ##################################################
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 ><b>Image parameters </b></font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        # Tissues Quality ###################################################
        self.report.append(
            Paragraph(
                "<font size = 15 ><b>TISSUES QUALITY </b></font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(
            Paragraph(
                "Metrics that do not fall into the above "
                "categories: statistical properties of "
                "tissue distributions, volume overlap "
                "of tissues, image harpness/blurriness,"
                " etc.",
                self.styles["Left2"],
            )
        )
        self.report.append(Spacer(0 * mm, 10 * mm))
        # summary_csf_mean
        self.report.append(self.get_iqms_data("summary_csf_mean"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_stdv
        self.report.append(self.get_iqms_data("summary_csf_stdv"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_median
        self.report.append(self.get_iqms_data("summary_csf_median"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_mad
        self.report.append(self.get_iqms_data("summary_csf_mad"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_p95
        self.report.append(self.get_iqms_data("summary_csf_p95"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_p05
        self.report.append(self.get_iqms_data("summary_csf_p05"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_k
        self.report.append(self.get_iqms_data("summary_csf_k"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_n
        self.report.append(self.get_iqms_data("summary_csf_n"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # summary_gm_mean
        self.report.append(self.get_iqms_data("summary_gm_mean"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_stdv
        self.report.append(self.get_iqms_data("summary_gm_stdv"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_median
        self.report.append(self.get_iqms_data("summary_gm_median"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_mad
        self.report.append(self.get_iqms_data("summary_gm_mad"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_p95
        self.report.append(self.get_iqms_data("summary_gm_p95"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_p05
        self.report.append(self.get_iqms_data("summary_gm_p05"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_k
        self.report.append(self.get_iqms_data("summary_gm_k"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_n
        self.report.append(self.get_iqms_data("summary_gm_n"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # summary_wm_mean
        self.report.append(self.get_iqms_data("summary_wm_mean"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_stdv
        self.report.append(self.get_iqms_data("summary_wm_stdv"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_median
        self.report.append(self.get_iqms_data("summary_wm_median"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_mad
        self.report.append(self.get_iqms_data("summary_wm_mad"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_p95
        self.report.append(self.get_iqms_data("summary_wm_p95"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_p05
        self.report.append(self.get_iqms_data("summary_wm_p05"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_k
        self.report.append(self.get_iqms_data("summary_wm_k"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_n
        self.report.append(self.get_iqms_data("summary_wm_n"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # summary_bg_mean
        self.report.append(self.get_iqms_data("summary_bg_mean"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_stdv
        self.report.append(self.get_iqms_data("summary_bg_stdv"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_median
        self.report.append(self.get_iqms_data("summary_bg_median"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_mad
        self.report.append(self.get_iqms_data("summary_bg_mad"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_p95
        self.report.append(self.get_iqms_data("summary_bg_p95"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_p05
        self.report.append(self.get_iqms_data("summary_bg_p05"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_k
        self.report.append(self.get_iqms_data("summary_bg_k"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_n
        self.report.append(self.get_iqms_data("summary_bg_n"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # fwhm_x
        self.report.append(self.get_iqms_data("fwhm_x"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_y
        self.report.append(self.get_iqms_data("fwhm_y"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_z
        self.report.append(self.get_iqms_data("fwhm_z"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_avg
        self.report.append(self.get_iqms_data("fwhm_avg"))
        self.report.append(Spacer(0 * mm, 7 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("summary_csf_k", "fwhm_x"):
            self.report.append(
                Paragraph(
                    "<font size = 8><sup>"
                    + self.dict4runtime["extra_info"][param][2]
                    + "</sup>"
                    + self.dict4runtime["extra_info"][param][3]
                    + "</font>",
                    self.styles["Left"],
                )
            )

        self.report.append(PageBreak())

        # Fifth page page - IQMs###############################################
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 ><b>Image parameters</b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        # Tissues Quality ###################################################
        self.report.append(
            Paragraph(
                "<font size = 15 ><b>TISSUES QUALITY</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(
            Paragraph(
                "Metrics that do not fall into the above categories: "
                "statistical properties of tissue distributions, "
                "volume overlap of tissues, image harpness/blurriness, "
                "etc.",
                self.styles["Left2"],
            )
        )
        self.report.append(Spacer(0 * mm, 10 * mm))
        # icvs_csf
        self.report.append(self.get_iqms_data("icvs_csf"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # icvs_gm
        self.report.append(self.get_iqms_data("icvs_gm"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # icvs_wm
        self.report.append(self.get_iqms_data("icvs_wm"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # rpve_csf
        self.report.append(self.get_iqms_data("rpve_csf"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # rpve_gm
        self.report.append(self.get_iqms_data("rpve_gm"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # rpve_wm
        self.report.append(self.get_iqms_data("rpve_wm"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # tpm_overlap_csf
        self.report.append(self.get_iqms_data("tpm_overlap_csf"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # tpm_overlap_gm
        self.report.append(self.get_iqms_data("tpm_overlap_gm"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # tpm_overlap_wm
        self.report.append(self.get_iqms_data("tpm_overlap_wm"))
        self.report.append(Spacer(0 * mm, 147 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("rpve_csf", "tpm_overlap_csf"):
            self.report.append(
                Paragraph(
                    "<font size = 8><sup>"
                    + self.dict4runtime["extra_info"][param][2]
                    + "</sup>"
                    + self.dict4runtime["extra_info"][param][3]
                    + "</font>",
                    self.styles["Left"],
                )
            )

        self.report.append(PageBreak())

        # Sixth page - slice planes display - Raw anatomic ####################
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 > <b>MRI axial slice "
                "planes display</b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                "<font size = 14 >Raw anatomic images</font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                '<font size = 9 ><i>"Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain.</i></font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        tmpdir = tempfile.TemporaryDirectory()
        slices_image = plot_slice_planes(
            data_1=self.anat,
            fig_rows=self.anat_fig_rows,
            fig_cols=self.anat_fig_cols,
            slice_start=self.anat_inf_slice_start,
            slice_step=self.anat_slices_gap,
            cmap_1="Greys_r",
            out_dir=tmpdir.name,
        )
        # reminder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)
        self.report.append(PageBreak())

        # Seventh page - slice planes display - Normalised anatomic (MNI) #####
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 > <b>MRI axial slice "
                "planes display</b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                "<font size = 14 >Normalised anatomic (MNI)</font>",
                self.styles["Center"],
            )
        )

        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                '<font size = 9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        slices_image = plot_slice_planes(
            data_1=self.norm_anat,
            fig_rows=self.norm_anat_fig_rows,
            fig_cols=self.norm_anat_fig_cols,
            slice_start=self.norm_anat_inf_slice_start,
            slice_step=self.norm_anat_slices_gap,
            cmap_1="Greys_r",
            out_dir=tmpdir.name,
        )
        # reminder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)
        self.report.append(PageBreak())

        # Eighth page - slice planes display - Background #####
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 > <b>MRI axial slice "
                "planes display</b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                "<font size = 14 > Anatomical "
                "image with background enhancement</font>",
                self.styles["Center"],
            )
        )

        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                '<font size = 9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        slices_image = plot_slice_planes(
            data_1=self.anat,
            fig_rows=self.anat_fig_rows,
            fig_cols=self.anat_fig_cols,
            slice_start=self.anat_inf_slice_start,
            slice_step=self.anat_slices_gap,
            cmap_1="viridis_r",
            out_dir=tmpdir.name,
            only_noise=True,
            out_name="background",
        )
        # reminder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)
        self.report.append(PageBreak())

        # Nineth page - slice planes display - Segmentation #####
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 > <b>MRI display</b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(
            Paragraph(
                "<font size = 14 > Anatomical "
                "image with segmentation</font>",
                self.styles["Center"],
            )
        )

        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(
            Paragraph(
                '<font size = 9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(
            Paragraph(
                "Brain tissue segmentation",
                self.styles["Left2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        plot_seg = plot_segmentation(
            self.anat,
            self.segmentation,
            out_dir=tmpdir.name,
            name="PlotSegmentation",
            cut_coords=9,
            display_mode="z",
            levels=[0.5, 1.5, 2.5],
            colors=["r", "g", "b"],
        )

        image = Image(plot_seg, width=7.4803 * inch, height=1.5 * inch)
        image.hAlign = "CENTER"
        self.report.append(image)

        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(
            Paragraph(
                "Brain mask",
                self.styles["Left2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))

        plot_bmask = plot_segmentation(
            self.anat,
            self.brain_mask,
            out_dir=tmpdir.name,
            name="PlotBrainmask",
            cut_coords=9,
            display_mode="z",
            levels=[0.5],
            colors=["r"],
        )

        image = Image(plot_bmask, width=7.4803 * inch, height=1.5 * inch)
        image.hAlign = "CENTER"
        self.report.append(image)

        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(
            Paragraph(
                "Hat mask",
                self.styles["Left2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))

        plot_airmask = plot_segmentation(
            self.anat,
            self.air_mask,
            out_dir=tmpdir.name,
            name="PlotAirmask",
            cut_coords=6,
            display_mode="x",
            levels=[0.5],
            colors=["r"],
        )

        image = Image(plot_airmask, width=7.4803 * inch, height=1.5 * inch)
        image.hAlign = "CENTER"
        self.report.append(image)

        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(
            Paragraph(
                "Head outline",
                self.styles["Left2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))

        plot_headmask = plot_segmentation(
            self.anat,
            self.head_mask,
            out_dir=tmpdir.name,
            name="PlotHeadmask",
            cut_coords=6,
            display_mode="x",
            levels=[0.5],
            colors=["r"],
        )

        image = Image(plot_headmask, width=7.4803 * inch, height=1.5 * inch)
        image.hAlign = "CENTER"
        self.report.append(image)

        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(
            Paragraph(
                "Artifacts in background",
                self.styles["Left2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))

        plot_artmask = plot_segmentation(
            self.anat,
            self.art_mask,
            out_dir=tmpdir.name,
            name="PlotArtmask",
            cut_coords=9,
            display_mode="z",
            levels=[0.5],
            colors=["r"],
            saturate=True,
        )

        image = Image(plot_artmask, width=7.4803 * inch, height=1.5 * inch)
        image.hAlign = "CENTER"
        self.report.append(image)

        self.report.append(PageBreak())

        # Tenth page - qi2 plot #####
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 > <b>Distribution of the noise in "
                "the background</b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 10 * mm))

        qi2 = plot_qi2(
            np.asarray(self.iqms_data["histogram_qi2_x_grid"]),
            np.asarray(self.iqms_data["histogram_qi2_ref_pdf"]),
            np.asarray(self.iqms_data["histogram_qi2_fit_pdf"]),
            np.asarray(self.iqms_data["histogram_qi2_ref_data"]),
            int(self.iqms_data["histogram_qi2_cutoff_idx"]),
            out_file=os.path.join(
                tmpdir.name,
                "qi2_plot.png",
            ),
        )
        image = Image(qi2, width=7.4803 * inch, height=5 * inch)
        image.hAlign = "CENTER"
        self.report.append(image)

        self.page.build(self.report, canvasmaker=PageNumCanvas)
        tmpdir.cleanup()

    def mriqc_func_make_report(
        self,
    ):
        """Make mriqc functional individual report"""

        # First page - cover ##################################################
        #######################################################################

        self.report.append(self.image_cov)
        # width, height
        self.report.append(Spacer(0 * mm, 8 * mm))
        self.report.append(Paragraph(self.header_title, self.styles["Center"]))
        self.report.append(Spacer(0 * mm, 10 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(Paragraph(self.title, self.styles["Center"]))
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(self.cover_data)
        self.report.append(Spacer(0 * mm, 6 * mm))
        self.report.append(
            Paragraph(self.textDisclaimer, self.styles["Justify"])
        )
        self.report.append(Spacer(0 * mm, 6 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        self.report.append(
            Paragraph(
                "<font size = 8><sup>$</sup>Esteban O et "
                "al., <i>MRIQC: Advancing the Automatic "
                "Prediction of Image Quality in MRI from"
                " Unseen Sites</i>, PLOS ONE 12(9)"
                ":e0184661.</font>",
                self.styles["Left"],
            )
        )
        self.report.append(PageBreak())

        # Second page - IQMs ##################################################
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 ><b>Image parameters </b></font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        # Spatial and temporal resolution ###################################
        self.report.append(
            Paragraph(
                "<font size = 15 > <b>SPATIAL AND "
                "TEMPORAL RESOLUTION</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 2 * mm))
        self.report.append(Paragraph("Length - Spacing", self.styles["Left2"]))
        self.report.append(Spacer(0 * mm, 10 * mm))

        # size_x
        self.report.append(self.get_iqms_data("size_x"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # size_y
        self.report.append(self.get_iqms_data("size_y"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # size_z
        self.report.append(self.get_iqms_data("size_z"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # size_t
        self.report.append(self.get_iqms_data("size_t"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # spacing_x
        self.report.append(self.get_iqms_data("spacing_x"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # spacing_y
        self.report.append(self.get_iqms_data("spacing_y"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # spacing_z
        self.report.append(self.get_iqms_data("spacing_z"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # spacing_tr
        self.report.append(self.get_iqms_data("spacing_tr"))
        self.report.append(Spacer(0 * mm, 20 * mm))
        # Measures for artifacts and other ##################################
        self.report.append(
            Paragraph(
                "<font size = 15 > <b>MEASURES FOR "
                "ARTIFACTS AND OTHER</b> </font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 10 * mm))
        # aor
        self.report.append(self.get_iqms_data("aor"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # aqi
        self.report.append(self.get_iqms_data("aqi"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # gsr_x
        self.report.append(self.get_iqms_data("gsr_x"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # gsr_y
        self.report.append(self.get_iqms_data("gsr_y"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # fd_mean
        self.report.append(self.get_iqms_data("fd_mean"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fd_num
        self.report.append(self.get_iqms_data("fd_num"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fd_perc
        self.report.append(self.get_iqms_data("fd_perc"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # dummy_trs
        self.report.append(self.get_iqms_data("dummy_trs"))
        self.report.append(Spacer(0 * mm, 82 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("gsr_x", "fd_mean"):
            self.report.append(
                Paragraph(
                    "<font size = 8><sup>"
                    + self.dict4runtime["extra_info"][param][2]
                    + "</sup>"
                    + self.dict4runtime["extra_info"][param][3]
                    + "</font>",
                    self.styles["Left"],
                )
            )

        self.report.append(PageBreak())

        # Third page - IQMs ###################################################
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 ><b>Image parameters </b></font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        # Measures for the spatial information ##############################
        self.report.append(
            Paragraph(
                "<font size = 15 ><b>MEASURES FOR THE "
                "SPATIAL INFORMATION</b></font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 10 * mm))
        # summary_bg_mean
        self.report.append(self.get_iqms_data("summary_bg_mean"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_stdv
        self.report.append(self.get_iqms_data("summary_bg_stdv"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_median
        self.report.append(self.get_iqms_data("summary_bg_median"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_mad
        self.report.append(self.get_iqms_data("summary_bg_mad"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_p95
        self.report.append(self.get_iqms_data("summary_bg_p95"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_p05
        self.report.append(self.get_iqms_data("summary_bg_p05"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_k
        self.report.append(self.get_iqms_data("summary_bg_k"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_n
        self.report.append(self.get_iqms_data("summary_bg_n"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # summary_fg_mean
        self.report.append(self.get_iqms_data("summary_fg_mean"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_stdv
        self.report.append(self.get_iqms_data("summary_fg_stdv"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_median
        self.report.append(self.get_iqms_data("summary_fg_median"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_mad
        self.report.append(self.get_iqms_data("summary_fg_mad"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_p95
        self.report.append(self.get_iqms_data("summary_fg_p95"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_p05
        self.report.append(self.get_iqms_data("summary_fg_p05"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_k
        self.report.append(self.get_iqms_data("summary_fg_k"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # summary_fg_n
        self.report.append(self.get_iqms_data("summary_fg_n"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # snr
        self.report.append(self.get_iqms_data("snr"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # fber
        self.report.append(self.get_iqms_data("fber"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # efc
        self.report.append(self.get_iqms_data("efc"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # fwhm_x
        self.report.append(self.get_iqms_data("fwhm_x"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_y
        self.report.append(self.get_iqms_data("fwhm_y"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_z
        self.report.append(self.get_iqms_data("fwhm_z"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_avg
        self.report.append(self.get_iqms_data("fwhm_avg"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        self.report.append(Spacer(0 * mm, 66 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("summary_bg_k", "fber", "efc", "fwhm_x"):
            self.report.append(
                Paragraph(
                    "<font size = 8><sup>"
                    + self.dict4runtime["extra_info"][param][2]
                    + "</sup>"
                    + self.dict4runtime["extra_info"][param][3]
                    + "</font>",
                    self.styles["Left"],
                )
            )

        self.report.append(PageBreak())

        # fourth page - IQMs ##################################################
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 ><b>Image parameters </b></font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 20 * mm))
        # Measures for the temporal information #############################
        self.report.append(
            Paragraph(
                "<font size = 15 ><b>MEASURES FOR THE "
                "TEMPORAL INFORMATION</b></font>",
                self.styles["Bullet1"],
            )
        )
        self.report.append(Spacer(0 * mm, 10 * mm))
        # dvars_nstd
        self.report.append(self.get_iqms_data("dvars_nstd"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # dvars_std
        self.report.append(self.get_iqms_data("dvars_std"))
        self.report.append(Spacer(0 * mm, 1 * mm))
        # dvars_vstd
        self.report.append(self.get_iqms_data("dvars_vstd"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # tsnr
        self.report.append(self.get_iqms_data("tsnr"))
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        # gcor
        self.report.append(self.get_iqms_data("gcor"))
        self.report.append(Spacer(0 * mm, 170 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))

        for param in ("dvars_nstd", "tsnr", "gcor"):
            self.report.append(
                Paragraph(
                    "<font size = 8><sup>"
                    + self.dict4runtime["extra_info"][param][2]
                    + "</sup>"
                    + self.dict4runtime["extra_info"][param][3]
                    + "</font>",
                    self.styles["Left"],
                )
            )

        self.report.append(PageBreak())

        # Fifth - Carpet plot #######################
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 > <b>Summary plot </b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                "<font size = 9 > <i>Summary plot, "
                "showing the slice-wise signal intensity at the extremes for "
                "the identification of spikes, the outliers metric, "
                "the DVARS, the FD and the carpetplot. "
                "The carpet plot rows correspond to voxelwise time series, "
                "and are separated into regions: cortical gray matter, "
                "deep gray matter, white matter and cerebrospinal fluid, "
                "cerebellum and the brain-edge or “crown”. The crown "
                "corresponds to the voxels located on a closed "
                "band around the brain. Carpet plot done using the head "
                "motion corrected functional image. </i></font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        # reminder: A4 == 210mmx297mm
        image = Image(
            self.IQMS_plot, width=7.4803 * inch, height=9.0551 * inch
        )
        image.hAlign = "CENTER"
        self.report.append(image)
        self.report.append(PageBreak())

        # Sixth - slice planes display - Raw functional #######################
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 > <b>MRI axial slice "
                "planes display</b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                "<font size = 14 >Raw functional "
                "images (1<sup>st</sup> dynamic)"
                "</font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                '<font size = 9 ><i>"Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain.</i></font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        tmpdir = tempfile.TemporaryDirectory()
        slices_image = plot_slice_planes(
            data_1=self.func,
            fig_rows=self.func_fig_rows,
            fig_cols=self.func_fig_cols,
            slice_start=self.func_inf_slice_start,
            slice_step=self.func_slices_gap,
            cmap_1="Greys_r",
            out_dir=tmpdir.name,
        )
        # reminder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)
        self.report.append(PageBreak())

        # Seventh page - slice planes display - Normalised functional (MNI) ###
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 > <b>MRI axial slice "
                "planes display</b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                "<font size = 14 >Normalised mean functional (MNI)</font>",
                self.styles["Center"],
            )
        )

        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                '<font size = 9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        slices_image = plot_slice_planes(
            data_1=self.norm_func,
            fig_rows=self.norm_func_fig_rows,
            fig_cols=self.norm_func_fig_cols,
            slice_start=self.norm_func_inf_slice_start,
            slice_step=self.norm_func_slices_gap,
            cmap_1="Greys_r",
            out_dir=tmpdir.name,
        )
        # reminder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)
        self.report.append(PageBreak())

        # Eighth page - slice planes display - stddev #####
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 > <b>MRI axial slice "
                "planes display</b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                "<font size = 14 > Temporal standard deviation</font>",
                self.styles["Center"],
            )
        )

        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                '<font size = 9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        slices_image = plot_slice_planes(
            data_1=self.stddev,
            fig_rows=self.func_fig_rows,
            fig_cols=self.func_fig_cols,
            slice_start=self.func_inf_slice_start,
            slice_step=self.func_slices_gap,
            cmap_1="viridis",
            out_dir=tmpdir.name,
        )
        # reminder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)
        self.report.append(PageBreak())

        # Nineth page - slice planes display - Background #####
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 > <b>MRI axial slice "
                "planes display</b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                "<font size = 14 > Average functional "
                "image with background enhancement</font>",
                self.styles["Center"],
            )
        )

        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(
            Paragraph(
                '<font size = 9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))
        slices_image = plot_slice_planes(
            data_1=self.func_mean,
            fig_rows=self.func_fig_rows,
            fig_cols=self.func_fig_cols,
            slice_start=self.func_inf_slice_start,
            slice_step=self.func_slices_gap,
            cmap_1="viridis_r",
            out_dir=tmpdir.name,
            only_noise=True,
        )
        # reminder: A4 == 210mmx297mm
        slices_image = Image(
            slices_image, width=7.4803 * inch, height=9.0551 * inch
        )
        slices_image.hAlign = "CENTER"
        self.report.append(slices_image)
        self.report.append(PageBreak())

        # Ninth page - slice planes display - Segmentation #####
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 > <b>MRI display</b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(
            Paragraph(
                "<font size = 14 > Mean unctional "
                "image with segmentation</font>",
                self.styles["Center"],
            )
        )

        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(
            Paragraph(
                '<font size = 9 > <i> "Neurological" '
                "convention, the left side of the "
                "image corresponds to the left side of "
                "the brain. </i> <br/> </font>",
                self.styles["Center"],
            )
        )

        self.report.append(Spacer(0 * mm, 5 * mm))
        self.report.append(
            Paragraph(
                "Brain mask",
                self.styles["Left2"],
            )
        )
        self.report.append(Spacer(0 * mm, 1 * mm))

        plot_bmask = plot_segmentation(
            self.func_mean,
            self.brain_mask,
            out_dir=tmpdir.name,
            name="PlotBrainmask",
            cut_coords=9,
            display_mode="z",
            levels=[0.5],
            colors=["r"],
        )

        image = Image(plot_bmask, width=7.4803 * inch, height=1.5 * inch)
        image.hAlign = "CENTER"
        self.report.append(image)

        self.page.build(self.report, canvasmaker=PageNumCanvas)
        tmpdir.cleanup()

    def mriqc_group_make_report(self):
        """Make mriqc group report"""

        # First page - cover ##################################################
        #######################################################################

        self.report.append(self.image_cov)
        # width, height
        self.report.append(Spacer(0 * mm, 8 * mm))
        self.report.append(Paragraph(self.header_title, self.styles["Center"]))
        self.report.append(Spacer(0 * mm, 10 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(Paragraph(self.title, self.styles["Center"]))
        self.report.append(Spacer(0 * mm, 10 * mm))
        self.report.append(self.cover_data)
        self.report.append(Spacer(0 * mm, 6 * mm))
        self.report.append(
            Paragraph(self.textDisclaimer, self.styles["Justify"])
        )
        self.report.append(Spacer(0 * mm, 6 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 2.5 * mm))
        self.report.append(
            Paragraph(
                "<font size = 8><sup>$</sup>Esteban O et "
                "al., <i>MRIQC: Advancing the Automatic "
                "Prediction of Image Quality in MRI from"
                " Unseen Sites</i>, PLOS ONE 12(9)"
                ":e0184661.</font>",
                self.styles["Left"],
            )
        )
        self.report.append(PageBreak())

        # Second page - Subject #########################################
        #######################################################################
        self.report.append(
            Paragraph(
                "<font size = 18 > <b>Subjects</b> </font>",
                self.styles["Center"],
            )
        )
        self.report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = "CENTER"
        self.report.append(line)
        self.report.append(Spacer(0 * mm, 4 * mm))

        data = [["Input data", "Site", "MRI scanner"]]
        for key in list(self.dict4runtime.keys()):
            data_file = []
            for i in ("FileName", "Site", "Spectro"):
                tag = self.dict4runtime[key][i]
                if len(tag) > 30:
                    tag = f"{tag[:30]}\n{tag[30:]}"
                data_file.append(tag)
            data.append(data_file)

        subjects_table = Table(data)
        subjects_table.setStyle(
            TableStyle(
                [
                    ("FONT", (0, 0), (-1, 0), "Times-Bold"),
                    ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ]
            )
        )
        self.report.append(subjects_table)
        self.report.append(PageBreak())

        modality = self.mriqc_group["modality"]
        if modality == "anat":
            # Third page - Boxplots Noise ####################################
            ##################################################################
            self.report.append(
                Paragraph(
                    "<font size = 18 > <b>Noise</b> </font>",
                    self.styles["Center"],
                )
            )
            self.report.append(Spacer(0 * mm, 4 * mm))
            line = ReportLine(150)
            line.hAlign = "CENTER"
            self.report.append(line)
            self.report.append(Spacer(0 * mm, 4 * mm))

            boxplot_snr_image = Image(
                self.mriqc_group["SNR"], width=3 * inch, height=2.5 * inch
            )
            boxplot_snrd_image = Image(
                self.mriqc_group["SNRD"], width=3 * inch, height=2.5 * inch
            )
            boxplot_cnr_image = Image(
                self.mriqc_group["CNR"], width=3 * inch, height=2.5 * inch
            )
            boxplot_cjv_image = Image(
                self.mriqc_group["CJV"], width=3 * inch, height=2.5 * inch
            )
            boxplot_qi2_image = Image(
                self.mriqc_group["QI"], width=3 * inch, height=2.5 * inch
            )

            noise_data = [
                [boxplot_snr_image, boxplot_snrd_image],
                [boxplot_cnr_image, boxplot_cjv_image],
                [boxplot_qi2_image],
            ]
            noise_table = Table(noise_data, [90 * mm, 70 * mm])
            noise_table.setStyle(
                TableStyle(
                    [
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ]
                )
            )
            self.report.append(noise_table)
            self.report.append(PageBreak())

            # 4 page - Boxplots  #############################################
            ##################################################################
            # FBER, EFC
            self.report.append(
                Paragraph(
                    (
                        "<font size = 18 > <b>Spatial "
                        "distribution</b> </font>"
                    ),
                    self.styles["Center"],
                )
            )
            self.report.append(Spacer(0 * mm, 4 * mm))
            line = ReportLine(150)
            line.hAlign = "CENTER"
            self.report.append(line)
            self.report.append(Spacer(0 * mm, 4 * mm))
            boxplot_fber_image = Image(
                self.mriqc_group["FBER"], width=3 * inch, height=2.5 * inch
            )
            boxplot_efc_image = Image(
                self.mriqc_group["EFC"], width=3 * inch, height=2.5 * inch
            )
            spatial_distribution_data = [
                [boxplot_fber_image, boxplot_efc_image]
            ]
            spatial_distribution_table = Table(
                spatial_distribution_data, [90 * mm, 70 * mm]
            )
            spatial_distribution_table.setStyle(
                TableStyle(
                    [
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ]
                )
            )
            self.report.append(spatial_distribution_table)
            self.report.append(Spacer(0 * mm, 4 * mm))
            self.report.append(Spacer(0 * mm, 4 * mm))

            # wm2max, inu_range, inu_med
            self.report.append(
                Paragraph(
                    "<font size = 18 > <b>ARTIFACTS</b> </font>",
                    self.styles["Center"],
                )
            )
            self.report.append(Spacer(0 * mm, 4 * mm))
            line = ReportLine(150)
            line.hAlign = "CENTER"
            self.report.append(line)
            self.report.append(Spacer(0 * mm, 4 * mm))
            boxplot_w2max_image = Image(
                self.mriqc_group["WM2MAX"], width=3 * inch, height=2.5 * inch
            )
            boxplot_inu_image = Image(
                self.mriqc_group["INU"], width=3 * inch, height=2.5 * inch
            )
            artifact_data = [[boxplot_w2max_image, boxplot_inu_image]]
            artifact_table = Table(artifact_data, [90 * mm, 70 * mm])
            artifact_table.setStyle(
                TableStyle(
                    [
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ]
                )
            )
            self.report.append(artifact_table)
            self.report.append(PageBreak())

            # 4 page - Boxplots  #############################################
            ##################################################################
            self.report.append(
                Paragraph(
                    "<font size = 18 > <b>Tissues quality</b> </font>",
                    self.styles["Center"],
                )
            )
            self.report.append(Spacer(0 * mm, 4 * mm))
            line = ReportLine(150)
            line.hAlign = "CENTER"
            self.report.append(line)
            self.report.append(Spacer(0 * mm, 4 * mm))

            # FWHM, ICVS, RPVE, TMP_OVERLAP, SUMMARY
            boxplot_fwhm_image = Image(
                self.mriqc_group["FWHM"], width=3 * inch, height=2.5 * inch
            )
            boxplot_icvs_image = Image(
                self.mriqc_group["ICVS"], width=3 * inch, height=2.5 * inch
            )
            boxplot_rpve_image = Image(
                self.mriqc_group["RPVE"], width=3 * inch, height=2.5 * inch
            )
            boxplot_tmp_image = Image(
                self.mriqc_group["TMP_OVERLAP"],
                width=3 * inch,
                height=2.5 * inch,
            )
            boxplot_wm_image = Image(
                self.mriqc_group["SUMMARY_WM"],
                width=3 * inch,
                height=2.5 * inch,
            )
            boxplot_csf_image = Image(
                self.mriqc_group["SUMMARY_CSF"],
                width=3 * inch,
                height=2.5 * inch,
            )
            boxplot_gm_image = Image(
                self.mriqc_group["SUMMARY_GM"],
                width=3 * inch,
                height=2.5 * inch,
            )
            boxplot_bg_image = Image(
                self.mriqc_group["SUMMARY_BG"],
                width=3 * inch,
                height=2.5 * inch,
            )

            tissues_data = [
                [boxplot_fwhm_image, boxplot_icvs_image],
                [boxplot_rpve_image, boxplot_tmp_image],
                [boxplot_wm_image, boxplot_csf_image],
                [boxplot_gm_image, boxplot_bg_image],
            ]
            tissues_data_table = Table(tissues_data, [90 * mm, 70 * mm])
            tissues_data_table.setStyle(
                TableStyle(
                    [
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ]
                )
            )

            self.report.append(tissues_data_table)

        elif modality == "bold":
            # Third page - Boxplots ##########################################
            ##################################################################
            self.report.append(
                Paragraph(
                    "<font size = 18 > <b>Artifact</b> </font>",
                    self.styles["Center"],
                )
            )
            self.report.append(Spacer(0 * mm, 4 * mm))
            line = ReportLine(150)
            line.hAlign = "CENTER"
            self.report.append(line)
            self.report.append(Spacer(0 * mm, 4 * mm))
            boxplot_aor_image = Image(
                self.mriqc_group["AOR"], width=3 * inch, height=2.5 * inch
            )
            boxplot_aqi_image = Image(
                self.mriqc_group["AQI"], width=3 * inch, height=2.5 * inch
            )
            boxplot_dummy_image = Image(
                self.mriqc_group["DUMMY"], width=3 * inch, height=2.5 * inch
            )
            boxplot_gsr_image = Image(
                self.mriqc_group["GSR"], width=3 * inch, height=2.5 * inch
            )
            boxplot_fd_image = Image(
                self.mriqc_group["FD"], width=3 * inch, height=2.5 * inch
            )
            boxplot_fd_num_image = Image(
                self.mriqc_group["FD_NUM"], width=3 * inch, height=2.5 * inch
            )
            boxplot_fd_perc_image = Image(
                self.mriqc_group["FD_PERC"], width=3 * inch, height=2.5 * inch
            )
            data = [
                [boxplot_aor_image, boxplot_aqi_image],
                [boxplot_dummy_image, boxplot_gsr_image],
                [boxplot_fd_image, boxplot_fd_num_image],
                [boxplot_fd_perc_image],
            ]
            table = Table(data, [90 * mm, 70 * mm])
            table.setStyle(
                TableStyle(
                    [
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ]
                )
            )
            self.report.append(table)
            self.report.append(PageBreak())

            # 4 page - Boxplots ##############################################
            ##################################################################
            self.report.append(
                Paragraph(
                    (
                        "<font size = 18 > <b>Temporal "
                        "information</b> </font>"
                    ),
                    self.styles["Center"],
                )
            )
            self.report.append(Spacer(0 * mm, 4 * mm))
            line = ReportLine(150)
            line.hAlign = "CENTER"
            self.report.append(line)
            self.report.append(Spacer(0 * mm, 4 * mm))

            boxplot_dvars_image = Image(
                self.mriqc_group["DVARS"], width=3 * inch, height=2.5 * inch
            )
            boxplot_dvarsn_image = Image(
                self.mriqc_group["DVARSN"], width=3 * inch, height=2.5 * inch
            )
            boxplot_tsnr_num_image = Image(
                self.mriqc_group["TSNR"], width=3 * inch, height=2.5 * inch
            )
            boxplot_gcor_perc_image = Image(
                self.mriqc_group["GCOR"], width=3 * inch, height=2.5 * inch
            )

            data = [
                [boxplot_dvars_image, boxplot_dvarsn_image],
                [boxplot_tsnr_num_image, boxplot_gcor_perc_image],
            ]
            table = Table(data, [90 * mm, 70 * mm])
            table.setStyle(
                TableStyle(
                    [
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ]
                )
            )
            self.report.append(table)
            self.report.append(PageBreak())

            # 5 page - Boxplots ##############################################
            ##################################################################
            self.report.append(
                Paragraph(
                    (
                        "<font size = 18 > <b>Spatial "
                        "information</b> </font>"
                    ),
                    self.styles["Center"],
                )
            )
            self.report.append(Spacer(0 * mm, 4 * mm))
            line = ReportLine(150)
            line.hAlign = "CENTER"
            self.report.append(line)
            self.report.append(Spacer(0 * mm, 4 * mm))

            boxplot_fwhm_image = Image(
                self.mriqc_group["FWHM"], width=3 * inch, height=2.5 * inch
            )
            boxplot_icvs_image = Image(
                self.mriqc_group["SNR"], width=3 * inch, height=2.5 * inch
            )
            boxplot_rpve_image = Image(
                self.mriqc_group["FBER"], width=3 * inch, height=2.5 * inch
            )
            boxplot_tmp_image = Image(
                self.mriqc_group["EFC"], width=3 * inch, height=2.5 * inch
            )
            boxplot_gm_image = Image(
                self.mriqc_group["SUMMARY_FG"],
                width=3 * inch,
                height=2.5 * inch,
            )
            boxplot_bg_image = Image(
                self.mriqc_group["SUMMARY_BG"],
                width=3 * inch,
                height=2.5 * inch,
            )

            data = [
                [boxplot_fwhm_image, boxplot_icvs_image],
                [boxplot_rpve_image, boxplot_tmp_image],
                [boxplot_gm_image, boxplot_bg_image],
            ]
            table = Table(data, [90 * mm, 70 * mm])
            table.setStyle(
                TableStyle(
                    [
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ]
                )
            )
            self.report.append(table)
            self.report.append(PageBreak())

        self.report.append(PageBreak())

        self.page.build(self.report, canvasmaker=PageNumCanvas)
