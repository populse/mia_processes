# -*- coding: utf-8 -*-
"""The reporting library of the mia_processes package.

The purpose of this module is to provide the reporting bricks necessary to
generate automatic report at the end of a pipeline calculation.

:Contains:
    :Class:
        - ReportAnatMriqc
        - ReportCO2inhalCvr
        - ReportFuncMriqc
        - ReportGE2REC
        - ReportGroupMriqc
        - ReportPerfDsc

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
import tempfile
from datetime import datetime

# nipype import
from nipype.interfaces.base import File, OutputMultiPath, Undefined, traits
from nipype.interfaces.spm.base import ImageFileSPM

# populse_mia import
from populse_mia.software_properties import Config
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# mia_processes import:
from mia_processes.utils import (
    Report,
    dict4runtime_update,
    get_dbFieldValue,
    mriqc_get_all_run,
    mriqc_group_iqms_tsv,
    plot_boxplot_points,
)


class ReportAnatMriqc(ProcessMIA):
    """
    *Generates the report for anatomical data in MRIQC pipeline*

    Please, see the complete documentation for the `ReportAnatMriqc brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/ReportAnatMriqc.html>`_

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
            "The number of lines for the anatomical slice planes plot"
        )

        anat_fig_cols_desc = (
            "The number of columns for the anatomical slice planes plot"
        )

        anat_inf_slice_start_desc = (
            "The first index displayed in anatomical slice planes plot"
        )

        anat_slices_gap_desc = (
            "Gap between slices in anatomical slice planes plot"
        )

        air_mask_desc = "Air mask"

        art_mask_desc = "Art mask"

        brain_mask_desc = "Brain mask"

        head_mask_desc = "Head mask"

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
            "Gap between slices in normalised anatomical slice planes plot"
        )

        segmentation_desc = "Segmentation file"

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
            "brain_mask",
            ImageFileSPM(
                copyfile=False,
                output=False,
                optional=False,
                desc=brain_mask_desc,
            ),
        )

        self.add_trait(
            "segmentation",
            ImageFileSPM(
                copyfile=False,
                output=False,
                optional=False,
                desc=segmentation_desc,
            ),
        )

        self.add_trait(
            "head_mask",
            ImageFileSPM(
                copyfile=False,
                output=False,
                optional=False,
                desc=head_mask_desc,
            ),
        )

        self.add_trait(
            "art_mask",
            ImageFileSPM(
                copyfile=False,
                output=False,
                optional=False,
                desc=art_mask_desc,
            ),
        )

        self.add_trait(
            "air_mask",
            ImageFileSPM(
                copyfile=False,
                output=False,
                optional=False,
                desc=air_mask_desc,
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

    def list_outputs(self, is_plugged=None, iteration=False):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :param iteration: the state, iterative or not, of the process.
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
            "Institution",
            "Manufacturer",
            "Manufacturer's Model",
            "SoftwareVersions",
            "Age",
        )

        # Site
        if self.dict4runtime["Institution"] in ("", "Undefined"):
            self.dict4runtime["Institution"] = (
                "Grenoble University Hospital - CLUNI (default value)"
            )

        # MriScanner
        if self.dict4runtime["Manufacturer"] in ("", "Undefined"):
            self.dict4runtime["Manufacturer"] = "Philips (default value)"
        if self.dict4runtime["Manufacturer's Model"] in ("", "Undefined"):
            self.dict4runtime["Manufacturer's Model"] = (
                "Achieva dStream (default value)"
            )

        # Generate an output name
        if self.anat and self.anat not in ["<undefined>", traits.Undefined]:
            self.outputs["report"] = os.path.join(
                self.output_directory,
                "{0}_anatomical_mriqcReport_{1}.pdf".format(
                    (
                        self.dict4runtime["PatientName"]
                        if self.dict4runtime["PatientName"] != "Undefined"
                        else "Undefined_name_ref"
                    ),
                    datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f")[:22],
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
                "Mean of the distribution of cerebrospinal fluid",
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
                "Median of the distribution of cerebrospinal fluid",
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
                "95% percentile of the distribution of cerebrospinal fluid",
                2,
                None,
            ],
            "summary_csf_p05": [
                "5% percentile of the distribution of cerebrospinal fluid",
                2,
                None,
            ],
            "summary_csf_k": [
                "Kurtosis of the distribution of cerebrospinal fluid",
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
                "Standard deviation of the distribution of gray matter",
                2,
                None,
            ],
            "summary_gm_median": [
                "Median of the distribution of gray matter",
                2,
                None,
            ],
            "summary_gm_mad": [
                "Median absolute deviation of the distribution gray matter",
                2,
                None,
            ],
            "summary_gm_p95": [
                "95% percentile of the distribution of gray matter",
                2,
                None,
            ],
            "summary_gm_p05": [
                "5% percentile of the distribution of gray matter",
                2,
                None,
            ],
            "summary_gm_k": [
                "Kurtosis of the distribution of gray matter",
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
                "Number of voxels in the distribution of gray matter",
                "{:.2e}",
                None,
            ],
            "summary_wm_mean": [
                "Mean of the distribution of white matter",
                2,
                None,
            ],
            "summary_wm_stdv": [
                "Standard deviation of the distribution of white matter",
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
                "95% percentile of the distribution of white matter",
                2,
                None,
            ],
            "summary_wm_p05": [
                "5% percentile of the distribution of white matter",
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
                "Number of voxels in the distribution of white matter",
                "{:.2e}",
                None,
            ],
            "summary_bg_mean": [
                "Mean of the distribution of background",
                2,
                None,
            ],
            "summary_bg_stdv": [
                "Standard deviation of the distribution of background",
                2,
                None,
            ],
            "summary_bg_median": [
                "Median of the distribution of background",
                2,
                None,
            ],
            "summary_bg_mad": [
                "Median absolute deviation of the distribution background",
                2,
                None,
            ],
            "summary_bg_p95": [
                "95% percentile of the distribution of background",
                2,
                None,
            ],
            "summary_bg_p05": [
                "5% percentile of the distribution of background",
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
                "Number of voxels in the distribution of background",
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
                "Residual partial voluming error for cerebrospinal fluid",
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

        self.tags_inheritance(
            self.anat,
            self.outputs["report"],
        )

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
            brain_mask=self.brain_mask,
            segmentation=self.segmentation,
            air_mask=self.air_mask,
            art_mask=self.art_mask,
            head_mask=self.head_mask,
            norm_anat=self.norm_anat,
            norm_anat_fig_rows=self.norm_anat_fig_rows,
            norm_anat_fig_cols=self.norm_anat_fig_cols,
            norm_anat_inf_slice_start=self.norm_anat_inf_slice_start,
            norm_anat_slices_gap=self.norm_anat_slices_gap,
        )

        report.make_report()


class ReportCO2inhalCvr(ProcessMIA):
    """
    *Generates report for CVR study by CO2 inhalation task*

    Please, see the complete documentation for the `ReportCO2inhalCvr brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/ReportCO2inhalCvr.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ReportCO2inhalCvr, self).__init__()

        # Third party software required for the execution of the brick
        self.requirement = []

        # Inputs description
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
            "Gap between slices in normalised anatomical slice planes plot"
        )

        norm_anat_cmap_desc = (
            "Colormap name used for the normalised anatomical plot "
            "(a string, default: 'Greys_r')"
        )

        norm_anat_vmin_desc = (
            "Minimum value in the data range covered by the color map"
        )

        norm_anat_vmax_desc = (
            "Maximum value in the data range covered by the color map"
        )

        norm_func_desc = (
            "An existing, uncompressed normalised functional image file "
            "(valid extensions: .nii, .nii.gz)"
        )

        norm_func_fig_rows_desc = (
            "The number of lines for the normalised functional slice "
            "planes plot"
        )

        norm_func_fig_cols_desc = (
            "The number of columns for the normalised functional slice "
            "planes plot"
        )

        norm_func_inf_slice_start_desc = (
            "The first index displayed in the normalised functional slice "
            "planes plot"
        )

        norm_func_slices_gap_desc = (
            "Gap between slices in the normalised functional slice planes plot"
        )

        norm_func_cmap_desc = (
            "Colormap name used for the normalised functional plot "
            "(a string, default: 'Greys_r')"
        )

        norm_func_vmin_desc = (
            "Minimum value in the data range covered by the color map ("
        )

        norm_func_vmax_desc = (
            "Maximum value in the data range covered by the color map"
        )

        realignment_parameters_desc = (
            "Estimation of translation and rotation parameters when "
            "realigning functional data (a pathlike object or a"
            "string representing a file, or a list "
            "of pathlike objects or strings "
            "representing a file)"
        )
        regressor_physio_desc = (
            "Time series data that represent the "
            "physiological processes thought to "
            "contribute to the CVR in the fMRI signal "
            "(a pathlike object or a string representing "
            "a file, or a list of pathlike objects or "
            "strings representing a file)"
        )
        beta_image_desc = (
            "The 1st estimated parameter of the model, i.e. beta_0001.nii "
            "(valid extensions:.nii)"
        )

        beta_cmap_desc = (
            "Colormap name used for the first beta parameter plot "
            "(a string, default: 'rainbow')"
        )

        beta_vmin_desc = (
            "Minimum value in the data range covered by the color map ("
        )

        beta_vmax_desc = (
            "Maximum value in the data range covered by the color map"
        )

        spmT_image_desc = (
            "Stat images from a t-contrast, i.e. spmT_0001.nii, "
            "(valid extensions:.nii)"
        )

        spmT_cmap_desc = (
            "Colormap name used for the stat images plot "
            "(a string, default: 'rainbow')"
        )

        spmT_vmin_desc = (
            "Minimum value in the data range covered by the color map"
        )

        spmT_vmax_desc = (
            "Maximum value in the data range covered by the color map"
        )

        mask_003_desc = (
            "The grey matter mask at the resolution of the smoothed_func."
        )

        patient_info_desc = (
            "Optional dictionary with information about the patient "
            "(e.g. {"
            "'PatientName': 'sub-1', 'Pathology': 'ACMD', "
            "'Age': 64, 'Sex': 'M', 'MR': '3T', "
            "'Gas': 'BACTAL', 'GasAdmin': 'MASK'}"
        )

        # Outputs description
        report_desc = "The generated report (pdf)"

        # Inputs traits
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

        self.add_trait(
            "norm_anat_cmap",
            traits.String(
                "Greys_r",
                output=False,
                optional=True,
                desc=norm_anat_cmap_desc,
            ),
        )

        self.add_trait(
            "norm_anat_vmin",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=norm_anat_vmin_desc,
            ),
        )
        self.norm_anat_vmin = Undefined

        self.add_trait(
            "norm_anat_vmax",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=norm_anat_vmax_desc,
            ),
        )
        self.norm_anat_vmax = Undefined

        # TODO: We use the normalized functional. It seems that it's the
        #       normalized - smothed functional that's used in Amigo. We'll
        #       have to check and decide what we'll use here finally?
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

        self.add_trait(
            "norm_func_cmap",
            traits.String(
                "Greys_r",
                output=False,
                optional=True,
                desc=norm_func_cmap_desc,
            ),
        )

        self.add_trait(
            "norm_func_vmin",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=norm_func_vmin_desc,
            ),
        )
        self.norm_func_vmin = Undefined

        self.add_trait(
            "norm_func_vmax",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=norm_func_vmax_desc,
            ),
        )
        self.norm_func_vmax = Undefined

        self.add_trait(
            "realignment_parameters",
            File(
                output=False, optional=False, desc=realignment_parameters_desc
            ),
        )

        self.add_trait(
            "regressor_physio",
            File(output=False, optional=False, desc=regressor_physio_desc),
        )

        self.add_trait(
            "beta_image",
            File(output=False, optional=False, desc=beta_image_desc),
        )

        self.add_trait(
            "beta_cmap",
            traits.String(
                "rainbow", output=False, optional=True, desc=beta_cmap_desc
            ),
        )

        self.add_trait(
            "beta_vmin",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=beta_vmin_desc,
            ),
        )
        self.beta_vmin = Undefined

        self.add_trait(
            "beta_vmax",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=beta_vmax_desc,
            ),
        )
        self.beta_vmax = Undefined

        self.add_trait(
            "spmT_image",
            File(output=False, optional=False, desc=spmT_image_desc),
        )

        self.add_trait(
            "spmT_cmap",
            traits.String(
                "rainbow", output=False, optional=True, desc=spmT_cmap_desc
            ),
        )

        self.add_trait(
            "spmT_vmin",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=spmT_vmin_desc,
            ),
        )
        self.spmT_vmin = Undefined

        self.add_trait(
            "spmT_vmax",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=spmT_vmax_desc,
            ),
        )
        self.spmT_vmax = Undefined

        self.add_trait(
            "mask_003",
            File(output=False, optional=False, desc=mask_003_desc),
        )

        self.add_trait(
            "patient_info",
            traits.Dict(output=False, optional=True, desc=patient_info_desc),
        )
        self.patient_info = dict(
            PatientRef=Undefined,
            Pathology=Undefined,
            Age=Undefined,
            Sex=Undefined,
            MR=Undefined,
            Gas=Undefined,
            GasAdmin=Undefined,
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

    def list_outputs(self, is_plugged=None, iteration=False):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :param iteration: the state, iterative or not, of the process.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(ReportCO2inhalCvr, self).list_outputs()

        if iteration is True:
            self.patient_info = dict(
                PatientRef=Undefined,
                Pathology=Undefined,
                Age=Undefined,
                Sex=Undefined,
                MR=Undefined,
                Gas=Undefined,
                GasAdmin=Undefined,
            )

        file_position = (
            self.norm_anat.find(self.project.getName())
            + len(self.project.getName())
            + 1
        )
        db_file_norm_anat = self.norm_anat[file_position:]
        db_file_norm_func = self.norm_func[file_position:]
        # TODO: Do we need to explicitly take the smoothed func as input,
        #       or can we simply add the prefix 's' to the normalized func?
        folder, fil = os.path.split(db_file_norm_func)
        db_file_smooth_norm_func = os.path.join(folder, "s" + fil)
        db_file_regressor_physio = self.regressor_physio[file_position:]

        # As we do not have access to the database at the runtime (see #272),
        # we prepare here the data that the run_process_mia method will need
        # via dict4runtime parameter:
        self.dict4runtime["norm_anat"] = {}
        dict4runtime_update(
            self.dict4runtime["norm_anat"],
            self.project.session,
            db_file_norm_anat,
            "AcquisitionDuration",
            "AcquisitionDate",
            "AcquisitionNumber",
            "Affine regularization type",
            "Age",
            "Dataset dimensions (Count, X,Y,Z,T...)",
            "EchoTime",
            "FlipAngle",
            "FOV",
            "Grid spacings (X,Y,Z,T,...)",
            "Institution",
            "Manufacturer",
            "Manufacturer's Model",
            "MaxNumOfSlices",
            "Pathology",
            "PatientName",
            "PatientRef",
            "ProtocolName",
            "RepetitionTime",
            "SequenceName",
            "Sex",
            "ScanResolution",
            "SoftwareVersions",
            "SliceGap",
            "SliceThickness",
            "Start/end slice",
            "StudyName",
            "Voxel sizes",
        )
        self.dict4runtime["norm_func"] = {}
        dict4runtime_update(
            self.dict4runtime["norm_func"],
            self.project.session,
            db_file_norm_func,
            "AcquisitionDuration",
            "AcquisitionNumber",
            "Affine regularization type",
            "Dataset dimensions (Count, X,Y,Z,T...)",
            "EchoTime",
            "FlipAngle",
            "FOV",
            "Gas",
            "GasAdmin",
            "Grid spacings (X,Y,Z,T,...)",
            "MaxNumOfSlices",
            "ProtocolName",
            "RepetitionTime",
            "SequenceName",
            "ScanResolution",
            "SliceGap",
            "SliceThickness",
            "Start/end slice",
            "Voxel sizes",
        )
        self.dict4runtime["smooth_norm_func"] = {}
        dict4runtime_update(
            self.dict4runtime["smooth_norm_func"],
            self.project.session,
            db_file_smooth_norm_func,
            "FWHM (X, Y, Z) for Smooth",
        )
        self.dict4runtime["regressor_physio"] = {}
        dict4runtime_update(
            self.dict4runtime["regressor_physio"],
            self.project.session,
            db_file_regressor_physio,
            "Regressor state",
        )
        # FIXME: the data should be anonymized and we should use PatientRef
        #        instead of PatientName !
        if self.dict4runtime["norm_anat"]["PatientName"] == "Undefined":
            print(
                "\nReportCO2inhalCvr brick:\nThe tags PatientName was not "
                "found in the database for the {} file...\n The "
                "initialization is aborted...".format(db_file_norm_anat)
            )
            return self.make_initResult()

        if (
            self.patient_info.get("PatientRef") is None
            or self.patient_info["PatientRef"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["PatientRef"] == "Undefined":
                self.dict4runtime["norm_anat"]["PatientRef"] = (
                    self.dict4runtime["norm_anat"]["PatientName"]
                )
                self.patient_info["PatientRef"] = self.dict4runtime[
                    "norm_anat"
                ]["PatientName"]

            else:
                self.patient_info["PatientRef"] = self.dict4runtime[
                    "norm_anat"
                ]["PatientRef"]

        else:
            self.dict4runtime["norm_anat"]["PatientRef"] = (
                self.patient_info.get("PatientRef")
            )

        if (
            self.patient_info.get("Pathology") is None
            or self.patient_info["Pathology"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Pathology"] != "Undefined":
                self.patient_info["Pathology"] = self.dict4runtime[
                    "norm_anat"
                ]["Pathology"]

        else:
            self.dict4runtime["norm_anat"]["Pathology"] = (
                self.patient_info.get("Pathology")
            )

        if (
            self.patient_info.get("Age") is None
            or self.patient_info["Age"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Age"] != "Undefined":
                self.patient_info["Age"] = self.dict4runtime["norm_anat"][
                    "Age"
                ]

        else:
            self.dict4runtime["norm_anat"]["Age"] = self.patient_info.get(
                "Age"
            )

        if (
            self.patient_info.get("Sex") is None
            or self.patient_info["Sex"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Sex"] != "Undefined":
                self.patient_info["Sex"] = self.dict4runtime["norm_anat"][
                    "Sex"
                ]

        else:
            self.dict4runtime["norm_anat"]["Sex"] = self.patient_info.get(
                "Sex"
            )

        if (
            self.patient_info.get("MR") is None
            or self.patient_info["MR"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Manufacturer"] != "Undefined":
                self.patient_info["MR"] = self.dict4runtime["norm_anat"][
                    "Manufacturer"
                ]

        else:
            self.dict4runtime["norm_anat"]["Manufacturer"] = (
                self.patient_info.get("MR")
            )

        if (
            self.patient_info.get("Gas") is None
            or self.patient_info["Gas"] == Undefined
        ):
            if self.dict4runtime["norm_func"]["Gas"] != "Undefined":
                self.patient_info["Gas"] == self.dict4runtime["norm_func"][
                    "Gas"
                ]

        else:
            self.dict4runtime["norm_func"]["Gas"] = self.patient_info.get(
                "Gas"
            )

        if (
            self.patient_info.get("GasAdmin") is None
            or self.patient_info["GasAdmin"] == Undefined
        ):
            if self.dict4runtime["norm_func"]["GasAdmin"] != "Undefined":
                self.patient_info["GasAdmin"] == self.dict4runtime[
                    "norm_func"
                ]["GasAdmin"]

        else:
            self.dict4runtime["norm_func"]["GasAdmin"] = self.patient_info.get(
                "GasAdmin"
            )

        # Site
        if self.dict4runtime["norm_anat"]["Institution"] in ("", "Undefined"):
            self.dict4runtime["norm_anat"][
                "Institution"
            ] = "Grenoble University Hospital - CLUNI (default value)"

        # MriScanner
        if self.dict4runtime["norm_anat"]["Manufacturer"] in ("", "Undefined"):
            self.dict4runtime["norm_anat"][
                "Manufacturer"
            ] = "Philips (default value)"
        if self.dict4runtime["norm_anat"]["Manufacturer's Model"] in (
            "",
            "Undefined",
        ):
            self.dict4runtime["norm_anat"][
                "Manufacturer's Model"
            ] = "Achieva dStream (default value)"

        # Reference population file for result stat anal
        config = Config()
        self.dict4runtime["CVR_ref_data"] = os.path.join(
            config.get_resources_path(),
            "reference_population_data",
            "CVR_ref_pop_SEMVIE-16.xlsx",
        )

        # Generate an output name
        if (
            self.norm_anat
            and self.norm_anat not in ["<undefined>", traits.Undefined]
        ) and (
            self.norm_func
            and self.norm_func not in ["<undefined>", traits.Undefined]
        ):
            self.outputs["report"] = os.path.join(
                self.output_directory,
                "{0}_CO2_inhal_CVR_Report_{1}.pdf".format(
                    self.dict4runtime["norm_anat"]["PatientRef"],
                    datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f")[:22],
                ),
            )

        else:
            return self.make_initResult()

        self.tags_inheritance(
            self.norm_anat,
            self.outputs["report"],
        )
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ReportCO2inhalCvr, self).run_process_mia()

        report = Report(
            self.report,
            self.dict4runtime,
            CVR=True,
            norm_anat=self.norm_anat,
            norm_anat_fig_rows=self.norm_anat_fig_rows,
            norm_anat_fig_cols=self.norm_anat_fig_cols,
            norm_anat_inf_slice_start=self.norm_anat_inf_slice_start,
            norm_anat_slices_gap=self.norm_anat_slices_gap,
            norm_anat_cmap=self.norm_anat_cmap,
            norm_anat_vmin=self.norm_anat_vmin,
            norm_anat_vmax=self.norm_anat_vmax,
            norm_func=self.norm_func,
            norm_func_fig_rows=self.norm_func_fig_rows,
            norm_func_fig_cols=self.norm_func_fig_cols,
            norm_func_inf_slice_start=self.norm_func_inf_slice_start,
            norm_func_slices_gap=self.norm_func_slices_gap,
            norm_func_cmap=self.norm_func_cmap,
            norm_func_vmin=self.norm_func_vmin,
            norm_func_vmax=self.norm_func_vmax,
            realignment_parameters=self.realignment_parameters,
            regressor_physio=self.regressor_physio,
            beta_image=self.beta_image,
            beta_cmap=self.beta_cmap,
            beta_vmin=self.beta_vmin,
            beta_vmax=self.beta_vmax,
            spmT_image=self.spmT_image,
            spmT_cmap=self.spmT_cmap,
            spmT_vmin=self.spmT_vmin,
            spmT_vmax=self.spmT_vmax,
            mask_003=self.mask_003,
            output_directory=self.output_directory,
        )

        report.make_report()


class ReportFuncMriqc(ProcessMIA):
    """
    *Generates the report for functional data in MRIQC pipeline*

    Please, see the complete documentation for the `ReportFuncMriqc
    brick in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/ReportFuncMriqc.html>`_

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
        brain_mask_desc = "Brain mask"

        IQMs_file_desc = "A .JSON file containing the IQMs"

        IQMs_plot_desc = (
            "A figure with carpet and outliers/dvars/FD/spikes plot"
        )

        func_desc = "An existing, functional image file"

        func_mean_desc = "An existing, mean functional image file"

        func_fig_rows_desc = (
            "The number of lines for the functional slice planes plot"
        )

        func_fig_cols_desc = (
            "The number of columns for the functional slice planes plot"
        )

        func_inf_slice_start_desc = (
            "The first index displayed in functional slice planes plot"
        )

        func_slices_gap_desc = (
            "Gap between slices in functional slice planes plot"
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
            "Gap between slices in normalised functional slice planes plot"
        )

        stddev_func_desc = "Functional standard deviation image"

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
            "func_mean",
            File(
                copyfile=False,
                output=False,
                optional=True,
                desc=func_mean_desc,
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
            "IQMs_plot",
            File(
                copyfile=False,
                output=False,
                optional=True,
                desc=IQMs_plot_desc,
            ),
        )

        self.add_trait(
            "brain_mask",
            ImageFileSPM(
                copyfile=False,
                output=False,
                optional=False,
                desc=brain_mask_desc,
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

        self.add_trait(
            "stddev_func",
            File(
                copyfile=False,
                output=False,
                optional=True,
                desc=stddev_func_desc,
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

    def list_outputs(self, is_plugged=None, iteration=False):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :param iteration: the state, iterative or not, of the process.
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
            "Institution",
            "Manufacturer",
            "Manufacturer's Model",
            "SoftwareVersions",
            "Age",
        )

        # Site
        if self.dict4runtime["Institution"] in ("", "Undefined"):
            self.dict4runtime["Institution"] = (
                "Grenoble University Hospital - CLUNI (default value)"
            )

        # MriScanner
        if self.dict4runtime["Manufacturer"] in ("", "Undefined"):
            self.dict4runtime["Manufacturer"] = "Philips (default value)"
        if self.dict4runtime["Manufacturer's Model"] in ("", "Undefined"):
            self.dict4runtime["Manufacturer's Model"] = (
                "Achieva dStream (default value)"
            )

        # Generate an output name
        if self.func and self.func not in ["<undefined>", traits.Undefined]:
            self.outputs["report"] = os.path.join(
                self.output_directory,
                "{0}_functional_mriqcReport_{1}.pdf".format(
                    (
                        self.dict4runtime["PatientName"]
                        if self.dict4runtime["PatientName"] != "Undefined"
                        else "Undefined_name_ref"
                    ),
                    datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f")[:22],
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
                "Standard deviation of the distribution of background",
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
                "95% percentile of the distribution of background",
                "{:.2e}",
                None,
            ],
            "summary_bg_p05": [
                "5% percentile of the distribution of background",
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
                "Number of voxels in the distribution of background",
                "{:.2e}",
                None,
            ],
            "summary_fg_mean": [
                "Mean of the distribution of foreground",
                "{:.2e}",
                None,
            ],
            "summary_fg_stdv": [
                "Standard deviation of the distribution of foreground",
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
                "95% percentile of the distribution of foreground",
                "{:.2e}",
                None,
            ],
            "summary_fg_p05": [
                "5% percentile of the distribution of foreground",
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
                "Number of voxels in the distribution of foreground",
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

        self.tags_inheritance(
            self.func,
            self.outputs["report"],
        )

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
            IQMS_plot=self.IQMs_plot,
            stddev=self.stddev_func,
            func_mean=self.func_mean,
            func_fig_rows=self.func_fig_rows,
            func_fig_cols=self.func_fig_cols,
            func_inf_slice_start=self.func_inf_slice_start,
            func_slices_gap=self.func_slices_gap,
            brain_mask=self.brain_mask,
            norm_func=self.norm_func,
            norm_func_fig_rows=self.norm_func_fig_rows,
            norm_func_fig_cols=self.norm_func_fig_cols,
            norm_func_inf_slice_start=self.norm_func_inf_slice_start,
            norm_func_slices_gap=self.norm_func_slices_gap,
        )

        report.make_report()


class ReportGE2REC(ProcessMIA):
    """
    Report for GE2REC pipeline

    Please, see the complete documentation for the `ReportGE2REC brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/ReportGE2REC.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ReportGE2REC, self).__init__()

        # Third party software required for the execution of the brick
        self.requirement = []

        # Inputs description
        norm_anat_desc = (
            "An existing, uncompressed normalised anatomical "
            "image file (valid extensions: .nii)"
        )

        norm_func_gene_desc = (
            "An existing, uncompressed normalised functional image file "
            "(valid extensions: .nii, .nii.gz)"
        )

        norm_func_reco_desc = (
            "An existing, uncompressed normalised functional image file "
            "(valid extensions: .nii, .nii.gz)"
        )

        norm_func_recall_desc = (
            "An existing, uncompressed normalised functional image file "
            "(valid extensions: .nii, .nii.gz)"
        )

        norm_func_mask_desc = "Mask"
        realignment_parameters_gene_desc = (
            "Estimation of translation and rotation parameters when "
            "realigning functional data (a pathlike object or a"
            "string representing a file, or a list "
            "of pathlike objects or strings "
            "representing a file)"
        )
        realignment_parameters_reco_desc = (
            "Estimation of translation and rotation parameters when "
            "realigning functional data (a pathlike object or a"
            "string representing a file, or a list "
            "of pathlike objects or strings "
            "representing a file)"
        )
        realignment_parameters_recall_desc = (
            "Estimation of translation and rotation parameters when "
            "realigning functional data (a pathlike object or a"
            "string representing a file, or a list "
            "of pathlike objects or strings "
            "representing a file)"
        )

        spmT_gene_desc = (
            "Stat images from a t-contrast from GENE task "
            "(valid extensions:.nii)"
        )

        spmT_reco_desc = (
            "Stat images from a t-contrast from GENE task "
            "(valid extensions:.nii)"
        )
        spmT_recall_desc = (
            "Stat images from a t-contrast from GENE task "
            "(valid extensions:.nii)"
        )
        spmT_gene_enco_desc = ()
        correct_response_desc = ()

        spmT_vmin_desc = (
            "Minimum value in the data range covered by the color map"
        )

        spmT_vmax_desc = (
            "Maximum value in the data range covered by the color map"
        )
        li_curves_desc = ""

        patient_info_desc = (
            "Optional dictionary with information about the patient "
            "(e.g. {"
            "'AcquisitionDate': '20240415', 'Pathology': 'ACMD', "
            "'Age': 64, 'Sex': 'M'}"
        )

        # Outputs description
        report_desc = "The generated report (pdf)"

        # Inputs traits
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
            "norm_func_gene",
            ImageFileSPM(
                copyfile=False,
                output=False,
                optional=False,
                desc=norm_func_gene_desc,
            ),
        )

        self.add_trait(
            "norm_func_reco",
            ImageFileSPM(
                copyfile=False,
                output=False,
                optional=False,
                desc=norm_func_reco_desc,
            ),
        )

        self.add_trait(
            "norm_func_recall",
            ImageFileSPM(
                copyfile=False,
                output=False,
                optional=False,
                desc=norm_func_recall_desc,
            ),
        )

        self.add_trait(
            "norm_func_mask",
            ImageFileSPM(
                copyfile=False,
                output=False,
                optional=False,
                desc=norm_func_mask_desc,
            ),
        )

        self.add_trait(
            "realignment_parameters_gene",
            File(
                output=False,
                optional=False,
                desc=realignment_parameters_gene_desc,
            ),
        )

        self.add_trait(
            "realignment_parameters_reco",
            File(
                output=False,
                optional=False,
                desc=realignment_parameters_reco_desc,
            ),
        )

        self.add_trait(
            "realignment_parameters_recall",
            File(
                output=False,
                optional=False,
                desc=realignment_parameters_recall_desc,
            ),
        )

        self.add_trait(
            "spmT_gene",
            File(output=False, optional=False, desc=spmT_gene_desc),
        )

        self.add_trait(
            "spmT_reco",
            File(output=False, optional=False, desc=spmT_reco_desc),
        )
        self.add_trait(
            "spmT_recall",
            File(output=False, optional=False, desc=spmT_recall_desc),
        )
        self.add_trait(
            "spmT_gene_enco",
            File(output=False, optional=True, desc=spmT_gene_enco_desc),
        )

        self.add_trait(
            "spmT_vmin",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=spmT_vmin_desc,
            ),
        )
        self.spmT_vmin = 2.5

        self.add_trait(
            "spmT_vmax",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=spmT_vmax_desc,
            ),
        )
        self.spmT_vmax = 5.0

        self.add_trait(
            "li_curves",
            traits.List(
                File(),
                output=False,
                optional=True,
                desc=li_curves_desc,
            ),
        )

        self.add_trait(
            "correct_response",
            File(
                output=False,
                optional=True,
                desc=correct_response_desc,
            ),
        )

        self.add_trait(
            "patient_info",
            traits.Dict(output=False, optional=True, desc=patient_info_desc),
        )
        self.patient_info = dict(
            Pathology=Undefined,
            Age=Undefined,
            Sex=Undefined,
            AcquisitionDate=Undefined,
            DominantHand=Undefined,
            LateralizationPathology=Undefined,
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

    def list_outputs(self, is_plugged=None, iteration=False):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :param iteration: the state, iterative or not, of the process.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(ReportGE2REC, self).list_outputs()

        file_position = (
            self.norm_anat.find(self.project.getName())
            + len(self.project.getName())
            + 1
        )
        db_file_norm_anat = self.norm_anat[file_position:]
        db_file_norm_func_gene = self.norm_func_gene[file_position:]
        db_file_norm_func_reco = self.norm_func_reco[file_position:]
        db_file_norm_func_recall = self.norm_func_recall[file_position:]

        # As we do not have access to the database at the runtime (see #272),
        # we prepare here the data that the run_process_mia method will need
        # via dict4runtime parameter:
        self.dict4runtime["norm_anat"] = {}
        dict4runtime_update(
            self.dict4runtime["norm_anat"],
            self.project.session,
            db_file_norm_anat,
            "AcquisitionDate",
            "Acquisition nbr",
            "Affine regularization type",
            "Age",
            "Dataset dimensions (Count, X,Y,Z,T...)",
            "EchoTime",
            "FlipAngle",
            "FOV",
            "Grid spacings (X,Y,Z,T,...)",
            "Institution",
            "Manufacturer",
            "Manufacturer's Model",
            "Pathology",
            "PatientName",
            "ProtocolName",
            "RepetitionTime",
            "SequenceName",
            "Sex",
            "SliceThickness",
            "Start/end slice",
            "StudyName",
            "Voxel sizes",
            "DominantHand",
            "LateralizationPathology",
        )

        if (
            self.patient_info.get("AcquisitionDate") is None
            or self.patient_info["AcquisitionDate"] == Undefined
        ):
            if (
                self.dict4runtime["norm_anat"]["AcquisitionDate"]
                != "Undefined"
            ):
                self.patient_info["AcquisitionDate"] = self.dict4runtime[
                    "norm_anat"
                ]["AcquisitionDate"]

        else:
            self.dict4runtime["norm_anat"]["AcquisitionDate"] = (
                self.patient_info.get("AcquisitionDate")
            )
        if (
            self.patient_info.get("Age") is None
            or self.patient_info["Age"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Age"] != "Undefined":
                self.patient_info["Age"] = self.dict4runtime["norm_anat"][
                    "Age"
                ]

        else:
            self.dict4runtime["norm_anat"]["Age"] = self.patient_info.get(
                "Age"
            )

        if (
            self.patient_info.get("Sex") is None
            or self.patient_info["Sex"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Sex"] != "Undefined":
                self.patient_info["Sex"] = self.dict4runtime["norm_anat"][
                    "Sex"
                ]

        else:
            self.dict4runtime["norm_anat"]["Sex"] = self.patient_info.get(
                "Sex"
            )
        if (
            self.patient_info.get("Pathology") is None
            or self.patient_info["Pathology"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Pathology"] != "Undefined":
                self.patient_info["Pathology"] = self.dict4runtime[
                    "norm_anat"
                ]["Pathology"]

        else:
            self.dict4runtime["norm_anat"]["Pathology"] = (
                self.patient_info.get("Pathology")
            )
        if (
            self.patient_info.get("LateralizationPathology") is None
            or self.patient_info["LateralizationPathology"] == Undefined
        ):
            self.patient_info["LateralizationPathology"] = self.dict4runtime[
                "norm_anat"
            ]["LateralizationPathology"]
        else:
            self.dict4runtime["norm_anat"]["LateralizationPathology"] = (
                self.patient_info.get("LateralizationPathology")
            )
        if (
            self.patient_info.get("DominantHand") is None
            or self.patient_info["DominantHand"] == Undefined
        ):
            self.patient_info["DominantHand"] = self.dict4runtime["norm_anat"][
                "DominantHand"
            ]
        else:
            self.dict4runtime["norm_anat"]["DominantHand"] = (
                self.patient_info.get("DominantHand")
            )

        if (
            self.patient_info.get("Age") is None
            or self.patient_info["Age"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Age"] != "Undefined":
                self.patient_info["Age"] = self.dict4runtime["norm_anat"][
                    "Age"
                ]

        else:
            self.dict4runtime["norm_anat"]["Age"] = self.patient_info.get(
                "Age"
            )

        if (
            self.patient_info.get("Sex") is None
            or self.patient_info["Sex"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Sex"] != "Undefined":
                self.patient_info["Sex"] = self.dict4runtime["norm_anat"][
                    "Sex"
                ]
        else:
            self.dict4runtime["norm_anat"]["Sex"] = self.patient_info.get(
                "Sex"
            )
        self.dict4runtime["norm_func_gene"] = {}
        dict4runtime_update(
            self.dict4runtime["norm_func_gene"],
            self.project.session,
            db_file_norm_func_gene,
            "Acquisition nbr",
            "Affine regularization type",
            "Dataset dimensions (Count, X,Y,Z,T...)",
            "EchoTime",
            "FlipAngle",
            "FOV",
            "Gas",
            "GasAdmin",
            "Grid spacings (X,Y,Z,T,...)",
            "ProtocolName",
            "RepetitionTime",
            "SequenceName",
            "SliceThickness",
            "Start/end slice",
            "Voxel sizes",
        )
        self.dict4runtime["norm_func_reco"] = {}
        dict4runtime_update(
            self.dict4runtime["norm_func_reco"],
            self.project.session,
            db_file_norm_func_reco,
            "Acquisition nbr",
            "Affine regularization type",
            "Dataset dimensions (Count, X,Y,Z,T...)",
            "EchoTime",
            "FlipAngle",
            "FOV",
            "Gas",
            "GasAdmin",
            "Grid spacings (X,Y,Z,T,...)",
            "ProtocolName",
            "RepetitionTime",
            "SequenceName",
            "SliceThickness",
            "Start/end slice",
            "Voxel sizes",
        )
        self.dict4runtime["norm_func_recall"] = {}
        dict4runtime_update(
            self.dict4runtime["norm_func_recall"],
            self.project.session,
            db_file_norm_func_recall,
            "Acquisition nbr",
            "Affine regularization type",
            "Dataset dimensions (Count, X,Y,Z,T...)",
            "EchoTime",
            "FlipAngle",
            "FOV",
            "Gas",
            "GasAdmin",
            "Grid spacings (X,Y,Z,T,...)",
            "ProtocolName",
            "RepetitionTime",
            "SequenceName",
            "SliceThickness",
            "Start/end slice",
            "Voxel sizes",
        )

        # Site
        if self.dict4runtime["norm_anat"]["Institution"] in ("", "Undefined"):
            self.dict4runtime["norm_anat"][
                "Institution"
            ] = "Grenoble University Hospital - CLUNI (default value)"

        # MriScanner
        if self.dict4runtime["norm_anat"]["Manufacturer"] in ("", "Undefined"):
            self.dict4runtime["norm_anat"][
                "Manufacturer"
            ] = "Philips (default value)"
        if self.dict4runtime["norm_anat"]["Manufacturer's Model"] in (
            "",
            "Undefined",
        ):
            self.dict4runtime["norm_anat"][
                "Manufacturer's Model"
            ] = "Achieva dStream (default value)"

        # Generate an output name
        if self.norm_anat and self.norm_anat not in [
            "<undefined>",
            traits.Undefined,
        ]:
            date = datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f")[:22]
            patient_name = get_dbFieldValue(
                self.project, self.norm_anat, "PatientName"
            )
            self.outputs["report"] = os.path.join(
                self.output_directory,
                f"{patient_name}_GE2REC_Report_" f"{date}.pdf",
            )
        else:
            return self.make_initResult()

        self.tags_inheritance(
            self.norm_anat,
            self.outputs["report"],
        )
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ReportGE2REC, self).run_process_mia()

        report = Report(
            self.report,
            self.dict4runtime,
            GE2REC=True,
            norm_anat=self.norm_anat,
            norm_func_gene=self.norm_func_gene,
            norm_func_reco=self.norm_func_reco,
            norm_func_recall=self.norm_func_recall,
            norm_func_mask=self.norm_func_mask,
            realignment_parameters_gene=self.realignment_parameters_gene,
            realignment_parameters_reco=self.realignment_parameters_reco,
            realignment_parameters_recall=self.realignment_parameters_recall,
            spmT_gene=self.spmT_gene,
            spmT_reco=self.spmT_reco,
            spmT_recall=self.spmT_recall,
            spmT_gene_enco=self.spmT_gene_enco,
            spmT_vmin=self.spmT_vmin,
            spmT_vmax=self.spmT_vmax,
            li_curves=self.li_curves,
            correct_response=self.correct_response,
            output_directory=self.output_directory,
        )

        report.make_report()


class ReportGroupMriqc(ProcessMIA):
    """
    *Generates a group report for mriqc pipelines*

    Please, see the complete documentation for the `ReportGroupMriqc brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/ReportGroupMriqc.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ReportGroupMriqc, self).__init__()

        # Third party software required for the execution of the brick
        self.requirement = []

        # Inputs description
        modality_desc = "Modality"

        # Outputs description
        report_pdf_desc = "The generated report (pdf)"
        out_tsv_desc = "The tsv file"

        # Inputs traits
        self.add_trait(
            "modality",
            traits.Enum(
                "anat", "bold", output=False, optional=True, desc=modality_desc
            ),
        )
        # Outputs traits
        self.add_trait("out_tsv", File(output=True, desc=out_tsv_desc))
        self.add_trait("report_pdf", File(output=True, desc=report_pdf_desc))

        # Special parameter used as a messenger for the run_process_mia method
        self.add_trait(
            "dict4runtime",
            traits.Dict(output=False, optional=True, userlevel=1),
        )

        self.init_default_traits()

    def list_outputs(self, is_plugged=None, iteration=False):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :param iteration: the state, iterative or not, of the process.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """

        # Using the inheritance to ProcessMIA class, list_outputs method
        super(ReportGroupMriqc, self).list_outputs()

        # Update dict4runtime
        # As we do not have access to the database at the runtime (see #272),
        # we prepare here the data that the run_process_mia method will need
        # via dict4runtime parameter
        # for mriqc group we need to get all subjects involved
        # and info about each subject
        files_name = mriqc_get_all_run(
            self.modality, self.project, self.output_directory
        )

        if files_name is None:
            print(
                "There is no mriqc report for this project."
                "Please used mriqc anat and func pipeline first"
            )
            return
        i = 0
        self.dict4runtime = {}
        for file_name in files_name:
            dict4runtime_file = {"FileName": os.path.basename(file_name)}
            file_position = (
                file_name.find(self.project.getName())
                + len(self.project.getName())
                + 1
            )
            database_filename = file_name[file_position:]
            dict4runtime_update(
                dict4runtime_file,
                self.project.session,
                database_filename,
                "PatientName",
                "StudyName",
                "AcquisitionDate",
                "Sex",
                "Institution",
                "Manufacturer",
                "Manufacturer's Model",
                "SoftwareVersions",
                "Age",
            )
            # Site
            if self.dict4runtime["norm_anat"]["Institution"] in (
                "",
                "Undefined",
            ):
                self.dict4runtime["norm_anat"][
                    "Institution"
                ] = "Grenoble University Hospital - CLUNI (default value)"

            # MriScanner
            if self.dict4runtime["norm_anat"]["Manufacturer"] in (
                "",
                "Undefined",
            ):
                self.dict4runtime["norm_anat"][
                    "Manufacturer"
                ] = "Philips (default value)"
            if self.dict4runtime["norm_anat"]["Manufacturer's Model"] in (
                "",
                "Undefined",
            ):
                self.dict4runtime["norm_anat"][
                    "Manufacturer's Model"
                ] = "Achieva dStream (default value)"

            self.dict4runtime[str(i)] = dict4runtime_file
            i += 1

        # Get report name
        date = datetime.now().strftime("_%Y_%m_%d_%H_%M")
        if self.output_directory:
            self.outputs["report_pdf"] = os.path.join(
                self.output_directory,
                "mriqc_group_report_" + self.modality + date + ".pdf",
            )
            self.outputs["out_tsv"] = os.path.join(
                self.output_directory,
                "mriqc_group_report_" + self.modality + date + ".tsv",
            )
        else:
            print("No output_directory was found, please select a project!\n")
            return

        # FIXME: Do we need tags inheritance ?

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ReportGroupMriqc, self).run_process_mia()
        # Get all info in a tsv file and in a panda dataframe
        df, out_tsv = mriqc_group_iqms_tsv(
            self.modality, self.output_directory
        )

        # Plot boxplot and save them in a '.png' image
        date = datetime.now().strftime("_%Y_%m_%d_%H_%M")
        iqms_df = {
            "EFC": df[["efc"]],
            "FBER": df[["fber"]],
            "FWHM": df[["fwhm_avg", "fwhm_x", "fwhm_y", "fwhm_z"]],
        }
        if self.modality == "anat":
            iqms_df.update(
                {
                    "SNR": df[["snr_csf", "snr_wm", "snr_gm"]],
                    "CJV": df[["cjv"]],
                    "CNR": df[["cnr"]],
                    "WM2MAX": df[["wm2max"]],
                    "SNRD": df[["snrd_csf", "snrd_wm", "snrd_gm"]],
                    "QI": df[["qi_1", "qi_2"]],
                    "ICVS": df[["icvs_csf", "icvs_gm", "icvs_wm"]],
                    "INU": df[["inu_range", "inu_med"]],
                    "RPVE": df[["rpve_csf", "rpve_gm", "rpve_wm"]],
                    "TMP_OVERLAP": df[
                        ["tpm_overlap_csf", "tpm_overlap_gm", "tpm_overlap_wm"]
                    ],
                    "SUMMARY_WM": df[
                        [
                            "summary_wm_mean",
                            "summary_wm_stdv",
                            "summary_wm_median",
                            "summary_wm_mad",
                            "summary_wm_p95",
                            "summary_wm_p05",
                            "summary_wm_k",
                            "summary_wm_n",
                        ]
                    ],
                    "SUMMARY_CSF": df[
                        [
                            "summary_csf_mean",
                            "summary_csf_stdv",
                            "summary_csf_median",
                            "summary_csf_mad",
                            "summary_csf_p95",
                            "summary_csf_p05",
                            "summary_csf_k",
                            "summary_csf_n",
                        ]
                    ],
                    "SUMMARY_GM": df[
                        [
                            "summary_gm_mean",
                            "summary_gm_stdv",
                            "summary_gm_median",
                            "summary_gm_mad",
                            "summary_gm_p95",
                            "summary_gm_p05",
                            "summary_gm_k",
                            "summary_gm_n",
                        ]
                    ],
                    "SUMMARY_BG": df[
                        [
                            "summary_bg_mean",
                            "summary_bg_stdv",
                            "summary_bg_median",
                            "summary_bg_mad",
                            "summary_bg_p95",
                            "summary_bg_p05",
                            "summary_bg_k",
                            "summary_bg_n",
                        ]
                    ],
                }
            )

        if self.modality == "bold":
            iqms_df.update(
                {
                    "SNR": df[["snr"]],
                    "GSR": df[["gsr_x", "gsr_y"]],
                    "DVARS": df[["dvars_std", "dvars_vstd"]],
                    "DVARSN": df[["dvars_nstd"]],
                    "FD": df[["fd_mean"]],
                    "FD_NUM": df[["fd_num"]],
                    "FD_PERC": df[["fd_perc"]],
                    "DUMMY": df[["dummy_trs"]],
                    "GCOR": df[["gcor"]],
                    "TSNR": df[["tsnr"]],
                    "AOR": df[["aor"]],
                    "AQI": df[["aqi"]],
                    "SUMMARY_FG": df[
                        [
                            "summary_fg_mean",
                            "summary_fg_stdv",
                            "summary_fg_p95",
                            "summary_fg_p05",
                            "summary_fg_k",
                        ]
                    ],
                    "SUMMARY_BG": df[
                        [
                            "summary_bg_mean",
                            "summary_bg_stdv",
                            "summary_bg_p95",
                            "summary_bg_p05",
                            "summary_bg_k",
                        ]
                    ],
                }
            )

        mriqc_iqms_group = {"modality": self.modality}

        tpm_dir = tempfile.TemporaryDirectory()

        for iqm in list(iqms_df.keys()):
            df_i = iqms_df[iqm]
            out_image_path = os.path.join(
                tpm_dir.name,
                "group_plot_" + self.modality + "_" + iqm + date + ".png",
            )
            boxplot_path = plot_boxplot_points(df_i, iqm, iqm, out_image_path)
            mriqc_iqms_group[iqm] = boxplot_path

        report = Report(
            self.report_pdf, self.dict4runtime, mriqc_group=mriqc_iqms_group
        )

        report.make_report()

        tpm_dir.cleanup()


class ReportPerfDsc(ProcessMIA):
    """
    *Generates report for Perfusion study using DSC-MRI*

    Please, see the complete documentation for the `ReportPerfDsc brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/ReportPerfDsc.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ReportPerfDsc, self).__init__()

        # Third party software required for the execution of the brick
        self.requirement = []

        # Inputs description
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
            "Gap between slices in normalised anatomical slice planes plot"
        )

        norm_anat_cmap_desc = (
            "Colormap name used for the normalised anatomical plot "
            "(a string, default: 'Greys_r')"
        )

        norm_anat_vmin_desc = (
            "Minimum value in the data range covered by the color map"
        )

        norm_anat_vmax_desc = (
            "Maximum value in the data range covered by the color map"
        )

        norm_func_desc = (
            "An existing, uncompressed normalised functional image file "
            "(valid extensions: .nii, .nii.gz)"
        )

        norm_func_fig_rows_desc = (
            "The number of lines for the normalised functional slice "
            "planes plot"
        )

        norm_func_fig_cols_desc = (
            "The number of columns for the normalised functional slice "
            "planes plot"
        )

        norm_func_inf_slice_start_desc = (
            "The first index displayed in the normalised functional slice "
            "planes plot"
        )

        norm_func_slices_gap_desc = (
            "Gap between slices in the normalised functional slice planes plot"
        )

        norm_func_cmap_desc = (
            "Colormap name used for the normalised functional plot "
            "(a string, default: 'Greys_r')"
        )

        norm_func_vmin_desc = (
            "Minimum value in the data range covered by the color map ("
        )

        norm_func_vmax_desc = (
            "Maximum value in the data range covered by the color map"
        )

        realignment_parameters_desc = (
            "Estimation of translation and rotation parameters when "
            "realigning functional data (a pathlike object or a string "
            "representing a file, or a list of pathlike objects or strings "
            "representing a file)"
        )

        CBV_image_desc = (
            "Cerebral Blood Volume image map (valid extensions:.nii): The "
            "volume of blood within a given volume of tissue (mL/100g)"
        )

        CBV_cmap_desc = (
            "Colormap name used for the CBV image plot (a string, "
            "default: 'rainbow')"
        )

        CBV_vmin_desc = (
            "Minimum value in the data range covered by the color map"
        )

        CBV_vmax_desc = (
            "Maximum value in the data range covered by the color map"
        )

        CBF_image_desc = (
            "Cerebral Blood Flow image map (valid extensions:.nii): The rate "
            "at which blood is delivered to a given volume of brain tissue, "
            "usually expressed in milliliters per 100 grams of brain tissue "
            "per minute (mL/100g/min)"
        )

        CBF_cmap_desc = (
            "Colormap name used for the CBF image plot (a string, "
            "default: 'rainbow')"
        )

        CBF_vmin_desc = (
            "Minimum value in the data range covered by the color map"
        )

        CBF_vmax_desc = (
            "Maximum value in the data range covered by the color map"
        )

        Tmax_image_desc = (
            "Tmax image map (valid extensions:.nii): Time to the maximum of "
            "the residue function. This parameter measures the delay between "
            "the arrival of the contrast agent in the arterial input function "
            "and its peak concentration in the tissue of interest (s)"
        )

        Tmax_cmap_desc = (
            "Colormap name used for the Tmax image plot (a string, "
            "default: 'rainbow')"
        )

        Tmax_vmin_desc = (
            "Minimum value in the data range covered by the color map"
        )

        Tmax_vmax_desc = (
            "Maximum value in the data range covered by the color map"
        )

        MTT_image_desc = (
            "Mean Transit Time image map (valid extensions:.nii): The "
            "average time it takes for blood to pass through a given volume "
            "of brain tissue (s)"
        )

        MTT_cmap_desc = (
            "Colormap name used for the MTT image plot (a string, "
            "default: 'rainbow')"
        )

        MTT_vmin_desc = (
            "Minimum value in the data range covered by the color map"
        )

        MTT_vmax_desc = (
            "Maximum value in the data range covered by the color map"
        )

        aif_file_desc = "The file containing the AIF data (a .json file)"

        patient_info_desc = (
            "Optional dictionary with information about the patient "
            "(e.g. {"
            "'PatientName': 'sub-1', 'Pathology': 'ACMD', "
            "'Age': 64, 'Sex': 'M', 'MR': '3T'}"
        )

        # Outputs description
        report_desc = "The generated report (pdf)"

        # Inputs traits
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

        self.add_trait(
            "norm_anat_cmap",
            traits.String(
                "Greys_r",
                output=False,
                optional=True,
                desc=norm_anat_cmap_desc,
            ),
        )

        self.add_trait(
            "norm_anat_vmin",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=norm_anat_vmin_desc,
            ),
        )
        self.norm_anat_vmin = Undefined

        self.add_trait(
            "norm_anat_vmax",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=norm_anat_vmax_desc,
            ),
        )
        self.norm_anat_vmax = Undefined

        # TODO: We use the normalized functional. It seems that it's the
        #       normalized - smothed functional that's used in Amigo. We'll
        #       have to check and decide what we'll use here finally?
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

        self.add_trait(
            "norm_func_cmap",
            traits.String(
                "Greys_r",
                output=False,
                optional=True,
                desc=norm_func_cmap_desc,
            ),
        )

        self.add_trait(
            "norm_func_vmin",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=norm_func_vmin_desc,
            ),
        )
        self.norm_func_vmin = Undefined

        self.add_trait(
            "norm_func_vmax",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=norm_func_vmax_desc,
            ),
        )
        self.norm_func_vmax = Undefined

        self.add_trait(
            "realignment_parameters",
            File(
                output=False, optional=False, desc=realignment_parameters_desc
            ),
        )

        self.add_trait(
            "CBV_image",
            File(output=False, optional=False, desc=CBV_image_desc),
        )

        self.add_trait(
            "CBV_cmap",
            traits.String(
                "rainbow", output=False, optional=True, desc=CBV_cmap_desc
            ),
        )

        self.add_trait(
            "CBV_vmin",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=CBV_vmin_desc,
            ),
        )
        self.CBV_vmin = Undefined

        self.add_trait(
            "CBV_vmax",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=CBV_vmax_desc,
            ),
        )
        self.CBV_vmax = Undefined

        self.add_trait(
            "CBF_image",
            File(output=False, optional=False, desc=CBF_image_desc),
        )

        self.add_trait(
            "CBF_cmap",
            traits.String(
                "rainbow", output=False, optional=True, desc=CBF_cmap_desc
            ),
        )

        self.add_trait(
            "CBF_vmin",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=CBF_vmin_desc,
            ),
        )
        self.CBF_vmin = Undefined

        self.add_trait(
            "CBF_vmax",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=CBF_vmax_desc,
            ),
        )
        self.CBF_vmax = Undefined

        self.add_trait(
            "Tmax_image",
            File(output=False, optional=False, desc=Tmax_image_desc),
        )

        self.add_trait(
            "Tmax_cmap",
            traits.String(
                "rainbow", output=False, optional=True, desc=Tmax_cmap_desc
            ),
        )

        self.add_trait(
            "Tmax_vmin",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=Tmax_vmin_desc,
            ),
        )
        self.Tmax_vmin = Undefined

        self.add_trait(
            "Tmax_vmax",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=Tmax_vmax_desc,
            ),
        )
        self.Tmax_vmax = Undefined

        self.add_trait(
            "MTT_image",
            File(output=False, optional=False, desc=MTT_image_desc),
        )

        self.add_trait(
            "MTT_cmap",
            traits.String(
                "rainbow", output=False, optional=True, desc=MTT_cmap_desc
            ),
        )

        self.add_trait(
            "MTT_vmin",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=MTT_vmin_desc,
            ),
        )
        self.MTT_vmin = Undefined

        self.add_trait(
            "MTT_vmax",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=MTT_vmax_desc,
            ),
        )
        self.MTT_vmax = Undefined

        self.add_trait(
            "aif_file",
            traits.Either(
                File(),
                Undefined,
                copyfile=False,
                output=False,
                optional=False,
                desc=aif_file_desc,
            ),
        )
        self.aif_file = traits.Undefined

        self.add_trait(
            "patient_info",
            traits.Dict(output=False, optional=True, desc=patient_info_desc),
        )
        self.patient_info = dict(
            PatientRef=Undefined,
            Pathology=Undefined,
            Age=Undefined,
            Sex=Undefined,
            MR=Undefined,
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

    def list_outputs(self, is_plugged=None, iteration=False):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :param iteration: the state, iterative or not, of the process.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(ReportPerfDsc, self).list_outputs()

        if iteration is True:
            self.patient_info = dict(
                PatientRef=Undefined,
                Pathology=Undefined,
                Age=Undefined,
                Sex=Undefined,
                MR=Undefined,
            )

        file_position = (
            self.norm_anat.find(self.project.getName())
            + len(self.project.getName())
            + 1
        )
        db_file_norm_anat = self.norm_anat[file_position:]
        db_file_norm_func = self.norm_func[file_position:]
        db_file_CBV_image = self.CBV_image[file_position:]
        # TODO: Do we need to explicitly take the smoothed func as input,
        #       or can we simply add the prefix 's' to the normalized func?
        folder, fil = os.path.split(db_file_norm_func)
        db_file_smooth_norm_func = os.path.join(folder, "s" + fil)

        # As we do not have access to the database at the runtime (see #272),
        # we prepare here the data that the run_process_mia method will need
        # via dict4runtime parameter:
        self.dict4runtime["norm_anat"] = {}
        dict4runtime_update(
            self.dict4runtime["norm_anat"],
            self.project.session,
            db_file_norm_anat,
            "AcquisitionDuration",
            "AcquisitionDate",
            "AcquisitionNumber",
            "Affine regularization type",
            "Age",
            "Dataset dimensions (Count, X,Y,Z,T...)",
            "EchoTime",
            "FlipAngle",
            "FOV",
            "Grid spacings (X,Y,Z,T,...)",
            "Institution",
            "Manufacturer",
            "Manufacturer's Model",
            "SoftwareVersions",
            "MaxNumOfSlices",
            "Pathology",
            "PatientName",
            "PatientRef",
            "ProtocolName",
            "RepetitionTime",
            "SequenceName",
            "Sex",
            "ScanResolution",
            "SliceGap",
            "SliceThickness",
            "Start/end slice",
            "StudyName",
            "Voxel sizes",
        )
        self.dict4runtime["norm_func"] = {}
        dict4runtime_update(
            self.dict4runtime["norm_func"],
            self.project.session,
            db_file_norm_func,
            "AcquisitionDuration",
            "AcquisitionNumber",
            "Affine regularization type",
            "Dataset dimensions (Count, X,Y,Z,T...)",
            "EchoTime",
            "FlipAngle",
            "FOV",
            "Grid spacings (X,Y,Z,T,...)",
            "MaxNumOfSlices",
            "ProtocolName",
            "RepetitionTime",
            "SequenceName",
            "ScanResolution",
            "SliceGap",
            "SliceThickness",
            "Start/end slice",
            "Voxel sizes",
        )
        self.dict4runtime["smooth_norm_func"] = {}
        dict4runtime_update(
            self.dict4runtime["smooth_norm_func"],
            self.project.session,
            db_file_smooth_norm_func,
            "FWHM (X, Y, Z) for Smooth",
        )
        self.dict4runtime["CBV_image"] = {}
        dict4runtime_update(
            self.dict4runtime["CBV_image"],
            self.project.session,
            db_file_CBV_image,
            "Perf Normalisation Factor",
        )
        # FIXME: the data should be anonymized and we should use PatientRef
        #        instead of PatientName !
        if self.dict4runtime["norm_anat"]["PatientName"] == "Undefined":
            print(
                "\nReportPerfDsc brick:\nThe PatientName tags was not "
                "found in the database for the {} file...\n The "
                "initialization is aborted...".format(db_file_norm_anat)
            )
            return self.make_initResult()

        if (
            self.patient_info.get("PatientRef") is None
            or self.patient_info["PatientRef"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["PatientRef"] == "Undefined":
                self.dict4runtime["norm_anat"]["PatientRef"] = (
                    self.dict4runtime["norm_anat"]["PatientName"]
                )
                self.patient_info["PatientRef"] = self.dict4runtime[
                    "norm_anat"
                ]["PatientName"]

            else:
                self.patient_info["PatientRef"] = self.dict4runtime[
                    "norm_anat"
                ]["PatientRef"]

        else:
            self.dict4runtime["norm_anat"]["PatientRef"] = (
                self.patient_info.get("PatientRef")
            )

        if (
            self.patient_info.get("Pathology") is None
            or self.patient_info["Pathology"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Pathology"] != "Undefined":
                self.patient_info["Pathology"] = self.dict4runtime[
                    "norm_anat"
                ]["Pathology"]

        else:
            self.dict4runtime["norm_anat"]["Pathology"] = (
                self.patient_info.get("Pathology")
            )

        if (
            self.patient_info.get("Age") is None
            or self.patient_info["Age"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Age"] != "Undefined":
                self.patient_info["Age"] = self.dict4runtime["norm_anat"][
                    "Age"
                ]

        else:
            self.dict4runtime["norm_anat"]["Age"] = self.patient_info.get(
                "Age"
            )

        if (
            self.patient_info.get("Sex") is None
            or self.patient_info["Sex"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Sex"] != "Undefined":
                self.patient_info["Sex"] = self.dict4runtime["norm_anat"][
                    "Sex"
                ]

        else:
            self.dict4runtime["norm_anat"]["Sex"] = self.patient_info.get(
                "Sex"
            )

        if (
            self.patient_info.get("MR") is None
            or self.patient_info["MR"] == Undefined
        ):
            if self.dict4runtime["norm_anat"]["Manufacturer"] != "Undefined":
                self.patient_info["MR"] = self.dict4runtime["norm_anat"][
                    "Manufacturer"
                ]

        else:
            self.dict4runtime["norm_anat"]["Manufacturer"] = (
                self.patient_info.get("MR")
            )

        # Site
        if self.dict4runtime["norm_anat"]["Institution"] in ("", "Undefined"):
            self.dict4runtime["norm_anat"][
                "Institution"
            ] = "Grenoble University Hospital - CLUNI (default value)"

        # MriScanner
        if self.dict4runtime["norm_anat"]["Manufacturer"] in ("", "Undefined"):
            self.dict4runtime["norm_anat"][
                "Manufacturer"
            ] = "Philips (default value)"
        if self.dict4runtime["norm_anat"]["Manufacturer's Model"] in (
            "",
            "Undefined",
        ):
            self.dict4runtime["norm_anat"][
                "Manufacturer's Model"
            ] = "Achieva dStream (default value)"
        # Generate an output name
        if (
            (
                self.norm_anat
                and self.norm_anat not in ["<undefined>", traits.Undefined]
            )
            and (
                self.norm_func
                and self.norm_func not in ["<undefined>", traits.Undefined]
            )
            and (
                self.CBV_image
                and self.CBV_image not in ["<undefined>", traits.Undefined]
            )
            and (
                self.CBF_image
                and self.CBF_image not in ["<undefined>", traits.Undefined]
            )
            and (
                self.Tmax_image
                and self.Tmax_image not in ["<undefined>", traits.Undefined]
            )
            and (
                self.MTT_image
                and self.MTT_image not in ["<undefined>", traits.Undefined]
            )
        ):
            self.outputs["report"] = os.path.join(
                self.output_directory,
                "{0}_Perf_DSC_Report_{1}.pdf".format(
                    self.dict4runtime["norm_anat"]["PatientRef"],
                    datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f")[:22],
                ),
            )

        else:
            return self.make_initResult()

        self.tags_inheritance(
            self.norm_anat,
            self.outputs["report"],
        )
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ReportPerfDsc, self).run_process_mia()

        report = Report(
            self.report,
            self.dict4runtime,
            Perf_DSC=True,
            norm_anat=self.norm_anat,
            norm_anat_fig_rows=self.norm_anat_fig_rows,
            norm_anat_fig_cols=self.norm_anat_fig_cols,
            norm_anat_inf_slice_start=self.norm_anat_inf_slice_start,
            norm_anat_slices_gap=self.norm_anat_slices_gap,
            norm_anat_cmap=self.norm_anat_cmap,
            norm_anat_vmin=self.norm_anat_vmin,
            norm_anat_vmax=self.norm_anat_vmax,
            norm_func=self.norm_func,
            norm_func_fig_rows=self.norm_func_fig_rows,
            norm_func_fig_cols=self.norm_func_fig_cols,
            norm_func_inf_slice_start=self.norm_func_inf_slice_start,
            norm_func_slices_gap=self.norm_func_slices_gap,
            norm_func_cmap=self.norm_func_cmap,
            norm_func_vmin=self.norm_func_vmin,
            norm_func_vmax=self.norm_func_vmax,
            realignment_parameters=self.realignment_parameters,
            CBF_image=self.CBF_image,
            CBF_cmap=self.CBF_cmap,
            CBF_vmin=self.CBF_vmin,
            CBF_vmax=self.CBF_vmax,
            CBV_image=self.CBV_image,
            CBV_cmap=self.CBV_cmap,
            CBV_vmin=self.CBV_vmin,
            CBV_vmax=self.CBV_vmax,
            Tmax_image=self.Tmax_image,
            Tmax_cmap=self.Tmax_cmap,
            Tmax_vmin=self.Tmax_vmin,
            Tmax_vmax=self.Tmax_vmax,
            MTT_image=self.MTT_image,
            MTT_cmap=self.MTT_cmap,
            MTT_vmin=self.MTT_vmin,
            MTT_vmax=self.MTT_vmax,
            aif_file=self.aif_file,
            output_directory=self.output_directory,
        )

        report.make_report()
