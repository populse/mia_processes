# -*- coding: utf-8 -*-
"""The reporting library of the mia_processes package.

The purpose of this module is to provide the reporting bricks necessary to
generate automatic report at the end of a pipeline calculation.

:Contains:
    :Class:
        - ReportAnatMriqc
        - ReportFuncMriqc

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

# Other import
import os
from datetime import datetime
from sys import path

# nipype import
from nipype.interfaces.base import (
    File,
    OutputMultiPath,
    Str,
    Undefined,
    traits,
)
from nipype.interfaces.spm.base import ImageFileSPM

# populse_mia import
from populse_mia.data_manager.project import COLLECTION_CURRENT
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

from mia_processes.utils import Report, dict4runtime_update

# mia_processes import:


class ReportAnatMriqc(ProcessMIA):
    """
    * Generates the report for anatomical data in MRIQC pipeline.

    Please, see the complete documentation for the `MRIQC_report brick
    in the populse.mia_processes website
    https://populse.github.io/mia_processes/html/documentation/bricks/reports/MRIQC_report.html

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ReportAnatMriqc, self).__init__()

        # Third party software required for the execution of the brick
        self.requirement = []

        # Inputs description

        IQMs_file_desc = "A .JSON file containing the IQMs"

        anat_desc = (
            "An existing, uncompressed anatomical image file (valid "
            "extensions: .nii, .nii.gz)"
        )

        anat_fig_rows_desc = (
            "The number of lines for the anatomical slice " "planes plot"
        )

        anat_fig_cols_desc = (
            "The number of columns for the anatomical slice " "planes plot"
        )

        anat_inf_slice_start_desc = (
            "The first index displayed in anatomical " "slice planes plot"
        )

        anat_slices_gap_desc = (
            "Gap between slices in anatomical slice planes " "plot"
        )

        norm_anat_desc = (
            "An existing, uncompressed normalised anatomical "
            "image file (valid extensions: .nii)"
        )

        norm_anat_fig_rows_desc = (
            "The number of lines for the normalised "
            "anatomical slice planes plot"
        )

        norm_anat_fig_cols_desc = (
            "The number of columns for the normalised "
            "anatomical slice planes plot"
        )

        norm_anat_inf_slice_start_desc = (
            "The first index displayed in "
            "normalised anatomical slice planes "
            "plot"
        )

        norm_anat_slices_gap_desc = (
            "Gap between slices in normalised " "anatomical slice planes plot"
        )

        # Outputs description
        report_desc = "The generated report (pdf)"

        # Inputs traits
        self.add_trait(
            "IQMs_file",
            File(
                copyfile=False,
                output=False,
                optional=False,
                desc=IQMs_file_desc,
            ),
        )

        self.add_trait(
            "anat",
            ImageFileSPM(
                copyfile=False, output=False, optional=False, desc=anat_desc
            ),
        )

        self.add_trait(
            "anat_fig_rows",
            traits.Int(
                5, output=False, optional=True, desc=anat_fig_rows_desc
            ),
        )

        self.add_trait(
            "anat_fig_cols",
            traits.Int(
                5, output=False, optional=True, desc=anat_fig_cols_desc
            ),
        )

        self.add_trait(
            "anat_inf_slice_start",
            traits.Either(
                Undefined,
                traits.Int,
                default=Undefined,
                output=False,
                optional=True,
                desc=anat_inf_slice_start_desc,
            ),
        )

        self.add_trait(
            "anat_slices_gap",
            traits.Either(
                Undefined,
                traits.Int,
                default=Undefined,
                output=False,
                optional=True,
                desc=anat_slices_gap_desc,
            ),
        )

        self.add_trait(
            "norm_anat",
            ImageFileSPM(
                copyfile=False,
                output=False,
                optional=False,
                desc=norm_anat_desc,
            ),
        )

        self.add_trait(
            "norm_anat_fig_rows",
            traits.Int(
                5, output=False, optional=True, desc=norm_anat_fig_rows_desc
            ),
        )

        self.add_trait(
            "norm_anat_fig_cols",
            traits.Int(
                5, output=False, optional=True, desc=norm_anat_fig_cols_desc
            ),
        )

        self.add_trait(
            "norm_anat_inf_slice_start",
            traits.Either(
                Undefined,
                traits.Int,
                default=Undefined,
                output=False,
                optional=True,
                desc=norm_anat_inf_slice_start_desc,
            ),
        )

        self.add_trait(
            "norm_anat_slices_gap",
            traits.Either(
                Undefined,
                traits.Int,
                default=Undefined,
                output=False,
                optional=True,
                desc=norm_anat_slices_gap_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "report",
            OutputMultiPath(
                File(), output=True, optional=True, desc=report_desc
            ),
        )
        # Special parameter used as a messenger for the run_process_mia method
        self.add_trait(
            "dict4runtime",
            traits.Dict(output=False, optional=True, userlevel=1),
        )

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
        super(ReportAnatMriqc, self).list_outputs()

        file_position = (
            self.anat.find(self.project.getName())
            + len(self.project.getName())
            + 1
        )
        database_filename = self.anat[file_position:]

        # As we do not have access to the database at the runtime (see #272),
        # we prepare here the data that the run_process_mia method will need
        # via dict4runtime parameter:
        dict4runtime_update(
            self.dict4runtime,
            self.project.session,
            database_filename,
            "PatientName",
            "StudyName",
            "AcquisitionDate",
            "Sex",
            "Site",
            "Spectro",
            "Age",
        )

        # FIXME: Currently, Site and Spectro data is hard-coded. A solution
        #        should be found to retrieve them automatically or to put them
        #        in the input parameters of the brick:
        # Site
        if self.dict4runtime["Site"] in ("", "Undefined"):
            self.dict4runtime["Site"] = "Grenoble University Hospital - CLUNI"

        # MriScanner
        if self.dict4runtime["Spectro"] in ("", "Undefined"):
            self.dict4runtime["Spectro"] = "Philips Achieva 3.0T TX"

        # Generate an output name
        if self.anat and self.anat not in ["<undefined>", traits.Undefined]:
            self.outputs["report"] = os.path.join(
                self.output_directory,
                "{0}_anatomical_mriqcReport_{1}.pdf".format(
                    self.dict4runtime["PatientName"]
                    if self.dict4runtime["PatientName"] != "Undefined"
                    else "Undefined_name_ref",
                    datetime.now().strftime("%Y_%m_%d_" "%H_%M_%S_%f")[:22],
                ),
            )

        # Add additional information for report generation
        # {iqms parameter: [parameter header, parameter rounding,
        #                   wildcard, related publication]
        self.dict4runtime["extra_info"] = {
            "size_x": ["Voxel size in X", 2, None],
            "size_y": ["Voxel size in Y", 2, None],
            "size_z": ["Voxel size in Z", 2, None],
            "spacing_x": ["Spacing in X (mm)", 2, None],
            "spacing_y": ["Spacing in Y (mm)", 2, None],
            "spacing_z": ["Spacing in Z (mm)", 2, None],
            "snr_csf": [
                "Signal-to-noise ratio (SNR) for cerebrospinal fluid",
                2,
                None,
            ],
            "snr_wm": ["SNR for white matter", 2, None],
            "snr_gm": ["SNR for gray matter", 2, None],
            "snr_total": ["SNR for brain parenchyma", 2, None],
            "snrd_csf": [
                "Dietrich’s SNR for cerebrospinal fluid",
                2,
                "$",
                "Dietrich et al., <i>Measurement of SNRs in MR images:"
                " influence of multichannel coils, parallel "
                "imaging and reconstruction filters</i>, JMRI "
                "26(2):375–385, 2007. Higher values are better.",
            ],
            "snrd_wm": [
                "Dietrich’s SNR for white matter",
                2,
                "$",
                "Dietrich et al., <i>Measurement of SNRs in MR images:"
                " influence of multichannel coils, parallel "
                "imaging and reconstruction filters</i>, JMRI "
                "26(2):375–385, 2007. Higher values are better.",
            ],
            "snrd_gm": [
                "Dietrich’s SNR for gray matter",
                2,
                "$",
                "Dietrich et al., <i>Measurement of SNRs in MR images:"
                " influence of multichannel coils, parallel "
                "imaging and reconstruction filters</i>, JMRI "
                "26(2):375–385, 2007. Higher values are better.",
            ],
            "snrd_total": [
                "Dietrich’s SNR for brain parenchyma",
                2,
                "$",
                "Dietrich et al., <i>Measurement of SNRs in MR "
                "images: influence of multichannel coils, "
                "parallel imaging and reconstruction filters</i>,"
                "JMRI 26(2):375–385, 2007. Higher values are "
                "better.",
            ],
            "cnr": [
                "Contrast-to-noise ratio",
                2,
                "#",
                "Magnotta, VA., & Friedman, L., <i>Measurement of "
                "signal-to-noise and contrast-to-noise in the fBIRN "
                "multicenter imaging study</i>, J Dig Imag "
                "19(2):140-147, 2006. Higher values are better.",
            ],
            "qi_2": [
                "Mortamet’s quality index 2",
                "{:.2e}",
                "&",
                "Mortamet B et al., <i>Automatic quality assessment in "
                "structural brain magnetic resonance imaging</i>, Mag "
                "Res Med 62(2):365-372, 2009. Lower values are better.",
            ],
            "cjv": [
                "Coefficient of joint variation",
                2,
                "%",
                "Ganzetti et al., <i> Intensity inhomogeneity correction "
                "of structural MR images: a data-driven approach to "
                "define input algorithm parameters</i>, Front "
                "Neuroinform 10:10, 2016. Lower values are better.",
            ],
            "fber": [
                "FBER",
                "{:.2e}",
                "%",
                "Shehzad Z et al., <i>The Preprocessed Connectomes "
                "Project Quality Assessment Protocol - a resource for "
                "measuring the quality of MRI data</i>, Front. "
                "Neurosci. Conference Abstract: Neuroinformatics 2015. "
                "Higher values are better.",
            ],
            "efc": [
                "EFC",
                2,
                "#",
                "Atkinson et al., <i>Automatic correction of motion "
                "artifacts in magnetic resonance images using an "
                "entropy focus criterion</i>, IEEE Trans Med Imag "
                "16(6):903-910, 1997. Lower values are better.",
            ],
            "wm2max": [
                "White-matter to maximum intensity ratio",
                2,
                "&",
                "The median intensity within the white matter mask over "
                "the 95% percentile of the full intensity "
                "distribution. Values should be around the interval "
                "[0.6, 0.8].",
            ],
            "qi_1": [
                "Mortamet’s quality index 1",
                "{:.2e}",
                "$",
                "Mortamet B et al., <i>Automatic quality assessment in "
                "structural brain magnetic resonance imaging</i>, Mag "
                "Res Med 62(2):365-372, 2009. Lower values are better.",
            ],
            "inu_range": [
                "Bias field range (95th percentile - 5th percentile)",
                2,
                "*",
                "Tustison NJ et al., <i>N4ITK: improved N3 bias "
                "correction</i>, IEEE Trans Med Imag, "
                "29(6):1310-20, 2010. Median closer to 1 and range "
                "closer to 0 are better.",
            ],
            "inu_med": [
                "Bias field median",
                2,
                "*",
                "Tustison NJ et al., <i>N4ITK: improved N3 bias "
                "correction</i>, IEEE Trans Med Imag, "
                "29(6):1310-20, 2010. Median closer to 1 and range "
                "closer to 0 are better.",
            ],
            "summary_csf_mean": [
                "Mean of the distribution of cerebrospinal " "fluid",
                2,
                None,
            ],
            "summary_csf_stdv": [
                "Standard deviation of the distribution of "
                "cerebrospinal fluid",
                2,
                None,
            ],
            "summary_csf_median": [
                "Median of the distribution of " "cerebrospinal fluid",
                2,
                None,
            ],
            "summary_csf_mad": [
                "Median absolute deviation of the distribution"
                " of cerebrospinal fluid",
                2,
                None,
            ],
            "summary_csf_p95": [
                "95% percentile of the distribution of " "cerebrospinal fluid",
                2,
                None,
            ],
            "summary_csf_p05": [
                "5% percentile of the distribution of " "cerebrospinal fluid",
                2,
                None,
            ],
            "summary_csf_k": [
                "Kurtosis of the distribution of " "cerebrospinal fluid",
                2,
                "*",
                "K is always ≥ -2. If the distribution is "
                "Gaussian, K = 0. If a distribution has less "
                "weight on its center and tails compared to a "
                "Gaussian of the same variance, then K < 0. If "
                "the distribution has more weight on its "
                "center and tails, then K > 0.",
            ],
            "summary_csf_n": [
                "Number of voxels in the distribution "
                "of cerebrospinal fluid",
                "{:.2e}",
                None,
            ],
            "summary_gm_mean": [
                "Mean of the distribution of gray matter",
                2,
                None,
            ],
            "summary_gm_stdv": [
                "Standard deviation of the distribution of " "gray matter",
                2,
                None,
            ],
            "summary_gm_median": [
                "Median of the distribution of gray matter",
                2,
                None,
            ],
            "summary_gm_mad": [
                "Median absolute deviation of the distribution" " gray matter",
                2,
                None,
            ],
            "summary_gm_p95": [
                "95% percentile of the distribution of " "gray matter",
                2,
                None,
            ],
            "summary_gm_p05": [
                "5% percentile of the distribution of " "gray matter",
                2,
                None,
            ],
            "summary_gm_k": [
                "Kurtosis of the distribution of " "gray matter",
                2,
                "*",
                "K is always ≥ -2. If the distribution is "
                "Gaussian, K = 0. If a distribution has less "
                "weight on its center and tails compared to a "
                "Gaussian of the same variance, then K < 0. If "
                "the distribution has more weight on its "
                "center and tails, then K > 0.",
            ],
            "summary_gm_n": [
                "Number of voxels in the distribution of gray " "matter",
                "{:.2e}",
                None,
            ],
            "summary_wm_mean": [
                "Mean of the distribution of white matter",
                2,
                None,
            ],
            "summary_wm_stdv": [
                "Standard deviation of the distribution of " "white matter",
                2,
                None,
            ],
            "summary_wm_median": [
                "Median of the distribution of white matter",
                2,
                None,
            ],
            "summary_wm_mad": [
                "Median absolute deviation of the distribution"
                " white matter",
                2,
                None,
            ],
            "summary_wm_p95": [
                "95% percentile of the distribution of " "white matter",
                2,
                None,
            ],
            "summary_wm_p05": [
                "5% percentile of the distribution of white " "matter",
                2,
                None,
            ],
            "summary_wm_k": [
                "Kurtosis of the distribution of white matter",
                2,
                "*",
                "K is always ≥ -2. If the distribution is "
                "Gaussian, K = 0. If a distribution has less "
                "weight on its center and tails compared to a "
                "Gaussian of the same variance, then K < 0. If "
                "the distribution has more weight on its "
                "center and tails, then K > 0.",
            ],
            "summary_wm_n": [
                "Number of voxels in the distribution of white " "matter",
                "{:.2e}",
                None,
            ],
            "summary_bg_mean": [
                "Mean of the distribution of background",
                2,
                None,
            ],
            "summary_bg_stdv": [
                "Standard deviation of the distribution of " "background",
                2,
                None,
            ],
            "summary_bg_median": [
                "Median of the distribution of background",
                2,
                None,
            ],
            "summary_bg_mad": [
                "Median absolute deviation of the distribution" " background",
                2,
                None,
            ],
            "summary_bg_p95": [
                "95% percentile of the distribution of " "background",
                2,
                None,
            ],
            "summary_bg_p05": [
                "5% percentile of the distribution of " "background",
                2,
                None,
            ],
            "summary_bg_k": [
                "Kurtosis of the distribution of background",
                2,
                "*",
                "K is always ≥ -2. If the distribution is "
                "Gaussian, K = 0. If a distribution has less "
                "weight on its center and tails compared to a "
                "Gaussian of the same variance, then K < 0. If "
                "the distribution has more weight on its "
                "center and tails, then K > 0.",
            ],
            "summary_bg_n": [
                "Number of voxels in the distribution of " "background",
                "{:.2e}",
                None,
            ],
            "fwhm_x": [
                "FWHM of the distribution in units of voxels, along x "
                "dimension",
                2,
                "#",
                "Forman SD et al., <i>Improved assessment of "
                "significant activation in functional magnetic "
                "resonance imaging(fMRI): use of a cluster - size "
                "threshold</i>, Magn.Reson.Med. 33(5), 636–647, 1995. "
                "Lower values are better, higher values indicate a "
                "blurrier image.",
            ],
            "fwhm_y": [
                "FWHM of the distribution in units of voxels, along y "
                "dimension",
                2,
                "#",
                "Forman SD et al., <i>Improved assessment of "
                "significant activation in functional magnetic "
                "resonance imaging(fMRI): use of a cluster - size "
                "threshold</i>, Magn.Reson.Med. 33(5), 636–647, 1995. "
                "Lower values are better, higher values indicate a "
                "blurrier image.",
            ],
            "fwhm_z": [
                "FWHM of the distribution in units of voxels, along z "
                "dimension",
                2,
                "#",
                "Forman SD et al., <i>Improved assessment of "
                "significant activation in functional magnetic "
                "resonance imaging(fMRI): use of a cluster - size "
                "threshold</i>, Magn.Reson.Med. 33(5), 636–647, 1995. "
                "Lower values are better, higher values indicate a "
                "blurrier image.",
            ],
            "fwhm_avg": [
                "FWHM of the distribution in units of voxels, average",
                2,
                "#",
                "Forman SD et al., <i>Improved assessment of "
                "significant activation in functional magnetic "
                "resonance imaging(fMRI): use of a cluster - size "
                "threshold</i>, Magn.Reson.Med. 33(5), 636–647, 1995. "
                "Lower values are better, higher values indicate a "
                "blurrier image.",
            ],
            "icvs_csf": [
                "Intracranial volume fraction for cerebrospinal fluid",
                2,
                None,
            ],
            "icvs_gm": [
                "Intracranial volume fraction for gray matter",
                2,
                None,
            ],
            "icvs_wm": [
                "Intracranial volume fraction for white matter",
                2,
                None,
            ],
            "rpve_csf": [
                "Residual partial voluming error for cerebrospinal " "fluid",
                2,
                "%",
                "Lower values are better.",
            ],
            "rpve_gm": [
                "Residual partial voluming error for gray matter",
                2,
                "%",
                "Lower values are better.",
            ],
            "rpve_wm": [
                "Residual partial voluming error for white matter",
                2,
                "%",
                "Lower values are better.",
            ],
            "tpm_overlap_csf": [
                "Overlap of the tissue probability map and the "
                "corresponding ICBM nonlinear-asymmetric "
                "2009c template, for cerebrospinal fluid",
                2,
                "&",
                "Higher values are better.",
            ],
            "tpm_overlap_gm": [
                "Overlap of the tissue probability map and the "
                "corresponding ICBM nonlinear-asymmetric "
                "2009c template, for gray matter",
                2,
                "&",
                "Higher values are better.",
            ],
            "tpm_overlap_wm": [
                "Overlap of the tissue probability map and the "
                "corresponding ICBM nonlinear-asymmetric "
                "2009c template, for white matter",
                2,
                "&",
                "Higher values are better.",
            ],
        }

        # FIXME: Do we need tags inheritance ?

        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ReportAnatMriqc, self).run_process_mia()

        report = Report(
            self.report,
            self.dict4runtime,
            IQMs_file=self.IQMs_file,
            anat=self.anat,
            anat_fig_rows=self.anat_fig_rows,
            anat_fig_cols=self.anat_fig_cols,
            anat_inf_slice_start=self.anat_inf_slice_start,
            anat_slices_gap=self.anat_slices_gap,
            norm_anat=self.norm_anat,
            norm_anat_fig_rows=self.norm_anat_fig_rows,
            norm_anat_fig_cols=self.norm_anat_fig_cols,
            norm_anat_inf_slice_start=self.norm_anat_inf_slice_start,
            norm_anat_slices_gap=self.norm_anat_slices_gap,
        )

        report.make_report()


class ReportFuncMriqc(ProcessMIA):
    """
    * Generates the report for functional data in MRIQC pipeline.

    Please, see the complete documentation for the `MRIQC_report brick in the populse.mia_processes website
    https://populse.github.io/mia_processes/html/documentation/bricks/reports/MRIQC_report.html

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ReportFuncMriqc, self).__init__()

        # Third party software required for the execution of the brick
        self.requirement = []

        # Inputs description

        IQMs_file_desc = "A .JSON file containing the IQMs"

        func_desc = (
            "An existing, uncompressed functional image file (valid "
            "extensions: .nii)"
        )

        func_fig_rows_desc = (
            "The number of lines for the functional slice " "planes plot"
        )

        func_fig_cols_desc = (
            "The number of columns for the functional slice " "planes plot"
        )

        func_inf_slice_start_desc = (
            "The first index displayed in functional " "slice planes plot"
        )

        func_slices_gap_desc = (
            "Gap between slices in functional slice planes " "plot"
        )

        norm_func_desc = (
            "An existing, uncompressed normalised functional "
            "image file (valid extensions: .nii)"
        )

        norm_func_fig_rows_desc = (
            "The number of lines for the normalised "
            "functional slice planes plot"
        )

        norm_func_fig_cols_desc = (
            "The number of columns for the normalised "
            "functional slice planes plot"
        )

        norm_func_inf_slice_start_desc = (
            "The first index displayed in "
            "normalised functional slice planes "
            "plot"
        )

        norm_func_slices_gap_desc = (
            "Gap between slices in normalised " "functional slice planes plot"
        )

        # Outputs description
        report_desc = "The generated report (pdf)"

        # Inputs traits
        self.add_trait(
            "IQMs_file",
            File(
                copyfile=False,
                output=False,
                optional=False,
                desc=IQMs_file_desc,
            ),
        )

        self.add_trait(
            "func",
            ImageFileSPM(
                copyfile=False, output=False, optional=False, desc=func_desc
            ),
        )

        self.add_trait(
            "func_fig_rows",
            traits.Int(
                5, output=False, optional=True, desc=func_fig_rows_desc
            ),
        )

        self.add_trait(
            "func_fig_cols",
            traits.Int(
                5, output=False, optional=True, desc=func_fig_cols_desc
            ),
        )

        self.add_trait(
            "func_inf_slice_start",
            traits.Either(
                Undefined,
                traits.Int,
                default=Undefined,
                output=False,
                optional=True,
                desc=func_inf_slice_start_desc,
            ),
        )

        self.add_trait(
            "func_slices_gap",
            traits.Either(
                Undefined,
                traits.Int,
                default=Undefined,
                output=False,
                optional=True,
                desc=func_slices_gap_desc,
            ),
        )

        self.add_trait(
            "norm_func",
            ImageFileSPM(
                copyfile=False,
                output=False,
                optional=False,
                desc=norm_func_desc,
            ),
        )

        self.add_trait(
            "norm_func_fig_rows",
            traits.Int(
                5, output=False, optional=True, desc=norm_func_fig_rows_desc
            ),
        )

        self.add_trait(
            "norm_func_fig_cols",
            traits.Int(
                5, output=False, optional=True, desc=norm_func_fig_cols_desc
            ),
        )

        self.add_trait(
            "norm_func_inf_slice_start",
            traits.Either(
                Undefined,
                traits.Int,
                default=Undefined,
                output=False,
                optional=True,
                desc=norm_func_inf_slice_start_desc,
            ),
        )

        self.add_trait(
            "norm_func_slices_gap",
            traits.Either(
                Undefined,
                traits.Int,
                default=Undefined,
                output=False,
                optional=True,
                desc=norm_func_slices_gap_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "report",
            OutputMultiPath(
                File(), output=True, optional=True, desc=report_desc
            ),
        )
        # Special parameter used as a messenger for the run_process_mia method
        self.add_trait(
            "dict4runtime",
            traits.Dict(output=False, optional=True, userlevel=1),
        )

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
        super(ReportFuncMriqc, self).list_outputs()

        file_position = (
            self.func.find(self.project.getName())
            + len(self.project.getName())
            + 1
        )
        database_filename = self.func[file_position:]

        # As we do not have access to the database at the runtime (see #272),
        # we prepare here the data that the run_process_mia method will need
        # via dict4runtime parameter:
        dict4runtime_update(
            self.dict4runtime,
            self.project.session,
            database_filename,
            "PatientName",
            "StudyName",
            "AcquisitionDate",
            "Sex",
            "Site",
            "Spectro",
            "Age",
        )

        # FIXME: Currently, Site and Spectro data are hard-coded. A solution
        #        should be found to retrieve them automatically or to put them
        #        in the input parameters of the brick:
        # Site
        if self.dict4runtime["Site"] in ("", "Undefined", None):
            self.dict4runtime["Site"] = "Grenoble University Hospital - CLUNI"

        # MriScanner
        if self.dict4runtime["Spectro"] in ("", "Undefined", None):
            self.dict4runtime["Spectro"] = "Philips Achieva 3.0T TX"

        # Generate an output name
        if self.func and self.func not in ["<undefined>", traits.Undefined]:
            self.outputs["report"] = os.path.join(
                self.output_directory,
                "{0}_functional_mriqcReport_{1}.pdf".format(
                    self.dict4runtime["PatientName"]
                    if self.dict4runtime["PatientName"] != "Undefined"
                    else "Undefined_name_ref",
                    datetime.now().strftime("%Y_%m_%d_" "%H_%M_%S_%f")[:22],
                ),
            )

        else:
            return self.make_initResult()

        # Add additional information for report generation
        # {iqms parameter: [parameter header, parameter rounding,
        #                   wildcard, related publication]
        self.dict4runtime["extra_info"] = {
            "size_x": ["Voxel size in X", 2, None],
            "size_y": ["Voxel size in Y", 2, None],
            "size_z": ["Voxel size in Z", 2, None],
            "size_t": ["Number of dynamics", 0, None],
            "spacing_x": ["Spacing in X (mm)", 2, None],
            "spacing_y": ["Spacing in Y (mm)", 2, None],
            "spacing_z": ["Spacing in Z (mm)", 2, None],
            "spacing_tr": ["Dynamic spacing (s)", 2, None],
            "gsr_x": [
                "Ghost to Signal Ratio (GSR), calculated along x",
                "{:.2e}",
                "$",
                "Giannelli et al., <i>Characterization of Nyquist ghost "
                "in EPI-fMRI acquisition sequences implemented on two "
                "clinical 1.5 T MR scanner systems: effect of readout "
                "bandwidth and echo spacing</i>, J App Clin Med Phy, "
                "11(4), 2010.",
            ],
            "gsr_y": [
                "Ghost to Signal Ratio (GSR), calculated along y",
                "{:.2e}",
                "$",
                "Giannelli et al., <i>Characterization of Nyquist ghost "
                "in EPI-fMRI acquisition sequences implemented on two "
                "clinical 1.5 T MR scanner systems: effect of readout "
                "bandwidth and echo spacing</i>, J App Clin Med Phy, "
                "11(4), 2010.",
            ],
            "fd_mean": [
                "Framewise Displacement (FD), average (mm)",
                2,
                "#",
                "Jenkinson et al., <i>Improved Optimisation for the "
                "Robust and Accurate Linear Registration and Motion "
                "Correction of Brain Images</i>, NeuroImage, 17(2), "
                "825-841, 2002.",
            ],
            "fd_num": [
                "Number of timepoints above FD threshold (0.2mm)",
                2,
                "#",
                "Jenkinson et al., <i>Improved Optimisation for the "
                "Robust and Accurate Linear Registration and Motion "
                "Correction of Brain Images</i>, NeuroImage, 17(2), "
                "825-841, 2002.",
            ],
            "fd_perc": [
                "Percent of FDs above the FD threshold (0.2mm)",
                2,
                "#",
                "Jenkinson et al., <i>Improved Optimisation for the "
                "Robust and Accurate Linear Registration and Motion "
                "Correction of Brain Images</i>, NeuroImage, 17(2), "
                "825-841, 2002.",
            ],
            "dummy_trs": [
                "Number of volumes identified as non-steady state",
                2,
                None,
            ],
            "summary_bg_mean": [
                "Mean of the distribution of background",
                "{:.2e}",
                None,
            ],
            "summary_bg_stdv": [
                "Standard deviation of the distribution of " "background",
                "{:.2e}",
                None,
            ],
            "summary_bg_median": [
                "Median of the distribution of background",
                2,
                None,
            ],
            "summary_bg_mad": [
                "Median absolute deviation of the distribution"
                " of background",
                2,
                None,
            ],
            "summary_bg_p95": [
                "95% percentile of the distribution of " "background",
                "{:.2e}",
                None,
            ],
            "summary_bg_p05": [
                "5% percentile of the distribution of " "background",
                2,
                None,
            ],
            "summary_bg_k": [
                "Kurtosis of the distribution of background",
                2,
                "*",
                "K is always ≥ -2. If the distribution is "
                "Gaussian, K = 0. If a distribution has less "
                "weight on its center and tails compared to a "
                "Gaussian of the same variance, then K < 0. If "
                "the distribution has more weight on its "
                "center and tails, then K > 0.",
            ],
            "summary_bg_n": [
                "Number of voxels in the distribution " "of background",
                "{:.2e}",
                None,
            ],
            "summary_fg_mean": [
                "Mean of the distribution of foreground",
                "{:.2e}",
                None,
            ],
            "summary_fg_stdv": [
                "Standard deviation of the distribution of " "foreground",
                "{:.2e}",
                None,
            ],
            "summary_fg_median": [
                "Median of the distribution of foreground",
                "{:.2e}",
                None,
            ],
            "summary_fg_mad": [
                "Median absolute deviation of the distribution"
                " of foreground",
                "{:.2e}",
                None,
            ],
            "summary_fg_p95": [
                "95% percentile of the distribution of " "foreground",
                "{:.2e}",
                None,
            ],
            "summary_fg_p05": [
                "5% percentile of the distribution of " "foreground",
                "{:.2e}",
                None,
            ],
            "summary_fg_k": [
                "Kurtosis of the distribution of foreground",
                2,
                "*",
                "K is always ≥ -2. If the distribution is "
                "Gaussian, K = 0. If a distribution has less "
                "weight on its center and tails compared to a "
                "Gaussian of the same variance, then K < 0. If "
                "the distribution has more weight on its "
                "center and tails, then K > 0.",
            ],
            "summary_fg_n": [
                "Number of voxels in the distribution " "of foreground",
                "{:.2e}",
                None,
            ],
            "snr": ["SNR for brain parenchyma", 2, None],
            "fber": [
                "FBER",
                "{:.2e}",
                "%",
                "Shehzad Z et al., <i>The Preprocessed Connectomes "
                "Project Quality Assessment Protocol - a resource for "
                "measuring the quality of MRI data</i>, Front. "
                "Neurosci. Conference Abstract: Neuroinformatics 2015. "
                "Higher values are better.",
            ],
            "efc": [
                "EFC",
                2,
                "#",
                "Atkinson et al., <i>Automatic correction of motion "
                "artifacts in magnetic resonance images using an "
                "entropy focus criterion</i>, IEEE Trans Med Imag "
                "16(6):903-910, 1997. Lower values are better.",
            ],
            "fwhm_x": [
                "FWHM of the distribution in units of voxels, along x "
                "dimension",
                2,
                "&",
                "Forman SD et al., <i>Improved assessment of "
                "significant activation in functional magnetic "
                "resonance imaging(fMRI): use of a cluster - size "
                "threshold</i>, Magn.Reson.Med. 33(5), 636–647, 1995. "
                "Lower values are better, higher values indicate a "
                "blurrier image.",
            ],
            "fwhm_y": [
                "FWHM of the distribution in units of voxels, along y "
                "dimension",
                2,
                "&",
                "Forman SD et al., <i>Improved assessment of "
                "significant activation in functional magnetic "
                "resonance imaging(fMRI): use of a cluster - size "
                "threshold</i>, Magn.Reson.Med. 33(5), 636–647, 1995. "
                "Lower values are better, higher values indicate a "
                "blurrier image.",
            ],
            "fwhm_z": [
                "FWHM of the distribution in units of voxels, along z "
                "dimension",
                2,
                "&",
                "Forman SD et al., <i>Improved assessment of "
                "significant activation in functional magnetic "
                "resonance imaging(fMRI): use of a cluster - size "
                "threshold</i>, Magn.Reson.Med. 33(5), 636–647, 1995. "
                "Lower values are better, higher values indicate a "
                "blurrier image.",
            ],
            "fwhm_avg": [
                "FWHM of the distribution in units of voxels, average",
                2,
                "&",
                "Forman SD et al., <i>Improved assessment of "
                "significant activation in functional magnetic "
                "resonance imaging(fMRI): use of a cluster - size "
                "threshold</i>, Magn.Reson.Med. 33(5), 636–647, 1995. "
                "Lower values are better, higher values indicate a "
                "blurrier image.",
            ],
            "dvars_nstd": [
                "DVARS, per-image standard deviation of the temporal "
                "derivative of the BOLD signal, scaled to 1000",
                2,
                "*",
                "Power et al., <i>Spurious but systematic "
                "correlations in functional connectivity MRI "
                "networks arise from subject motion</i>, NeuroImage "
                "59(3):2142-2154, 2012.",
            ],
            "dvars_std": [
                "Normalization of DVARS with the standard deviation "
                "of the temporal difference time series",
                2,
                "*",
                "Power et al., <i>Spurious but systematic "
                "correlations in functional connectivity MRI "
                "networks arise from subject motion</i>, NeuroImage "
                "59(3):2142-2154, 2012.",
            ],
            "dvars_vstd": [
                "Voxel-wise standardization of DVARS",
                2,
                "*",
                "Power et al., <i>Spurious but systematic "
                "correlations in functional connectivity MRI "
                "networks arise from subject motion</i>, NeuroImage "
                "59(3):2142-2154, 2012.",
            ],
            "tsnr": [
                "Temporal SNR",
                2,
                "#",
                "Krüger et al., <i>Physiological noise in "
                "oxygenation-sensitive magnetic resonance imaging</i>, "
                "Magn. Reson. Med. 46(4):631-637, 2001. Higher values are"
                "better.",
            ],
            "gcor": [
                "Global Correlation",
                2,
                "&",
                "Saad et al., <i>Correcting Brain-Wide Correlation "
                "Differences in Resting-State FMRI</i>, Brain Conn "
                "3(4):339-352, 2013.",
            ],
            "aor": ["AFNI's outlier ratio", "{:.2e}", None],
            "aqi": ["AFNI's quality index", "{:.2e}", None],
        }

        # We force tags inheritance (see mia_processes issue #13)
        self.inheritance_dict[self.outputs["report"]] = self.func

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ReportFuncMriqc, self).run_process_mia()

        report = Report(
            self.report,
            self.dict4runtime,
            IQMs_file=self.IQMs_file,
            func=self.func,
            func_fig_rows=self.func_fig_rows,
            func_fig_cols=self.func_fig_cols,
            func_inf_slice_start=self.func_inf_slice_start,
            func_slices_gap=self.func_slices_gap,
            norm_func=self.norm_func,
            norm_func_fig_rows=self.norm_func_fig_rows,
            norm_func_fig_cols=self.norm_func_fig_cols,
            norm_func_inf_slice_start=self.norm_func_inf_slice_start,
            norm_func_slices_gap=self.norm_func_slices_gap,
        )

        report.make_report()
