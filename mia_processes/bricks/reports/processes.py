# -*- coding: utf-8 -*-
"""The report preprocess library of the mia_processes package.

The purpose of this module is to provide bricks and functions to
compute necessary values for reporting.

:Contains:
    :Class:
        - AnatIQMs
        - BoldIQMs
        - BoldIQMsPlot
        - CarpetParcellation
        - ComputeDVARS
        - FramewiseDisplacement
        - Mean_stdDev_calc
        - Result_collector
        - Spikes

    :Function:
        - art_qi1
        - art_qi2
        - cjv
        - cnr
        - efc
        - fber
        - fuzzy_jaccard
        - find_peaks
        - find_spikes
        - gsr
        - image_binary_dilation
        - normalize_mc_params
        - regress_poly
        - rpve
        - snr
        - snr_dietrich
        - summary_stats
        - volume_fraction
        - wm2max
        - _AR_est_YW
        - _flatten_dict
        - _prepare_mask
        - _robust_zscore


"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# Other import
import json
import os
import re
import shutil
import tempfile
from math import sqrt

import nibabel as nb
import numpy as np
import pandas as pd
import scipy.ndimage as nd
from nilearn.signal import clean
from nipy.algorithms.registration import aff2euler, to_matrix44

# nipype import
from nipype.interfaces.base import (
    File,
    InputMultiPath,
    OutputMultiPath,
    Undefined,
    traits,
)
from nipype.interfaces.spm.base import ImageFileSPM
from nitime.algorithms import AR_est_YW
from niworkflows.utils.timeseries import _nifti_timeseries
from niworkflows.viz.plots import fMRIPlot
from numpy.polynomial import Legendre

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA
from scipy.stats import kurtosis  # pylint: disable=E0611
from skimage.transform import resize

# mia_processes import
from mia_processes.utils import (
    checkFileExt,
    get_dbFieldValue,
    set_dbFieldValue,
)

DIETRICH_FACTOR = 0.6551364  # 1.0 / sqrt(2 / (4 - pi))
FSL_FAST_LABELS = {"csf": 1, "gm": 2, "wm": 3, "bg": 0}
RAS_AXIS_ORDER = {"x": 0, "y": 1, "z": 2}
EXT = {"NIFTI_GZ": "nii.gz", "NIFTI": "nii"}


class AnatIQMs(ProcessMIA):
    """
    *Computes the anatomical IQMs*

    Please, see the complete documentation for the `AnatIQMs brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/AnatIQMs.html>`_

    adapted from `mriqc
    <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/anatomical.py#L332>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(AnatIQMs, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_ras_desc = (
            "RAS input image (a pathlike object or string "
            "representing a file)."
        )
        airmask_desc = (
            "Air mask image (a pathlike object or string "
            "representing a file)."
        )
        artmask_desc = (
            "Artifact mask image (a pathlike object or string "
            "representing a file)."
        )
        headmask_desc = (
            "Head mask image (a pathlike object or string "
            "representing a file)."
        )
        rotmask_desc = (
            "Rotation mask image (a pathlike object or string "
            "representing a file)."
        )
        hatmask_desc = (
            "Hat mask image (a pathlike object or string "
            "representing a file)."
        )
        segmentation_desc = (
            "Segmentation mask image (a pathlike object or "
            "string representing a file)."
        )
        in_inu_desc = (
            "Input INU image (a pathlike object or string "
            "representing a file)."
        )
        in_noinu_desc = (
            "Input no-INU image (a pathlike object or string "
            "representing a file)."
        )
        pvms_desc = (
            "PVMS image (a pathlike object or string " "representing a file)."
        )
        mni_tpms_desc = (
            "MNI TPMS file (a pathlike "
            "object or string representing a file)."
        )
        in_fwhm_desc = "FWHM (a float)."

        # Outputs description
        out_file_desc = "a json file containing IQMs"

        # Inputs traits
        self.add_trait(
            "in_ras", File(output=False, optional=False, desc=in_ras_desc)
        )
        self.add_trait(
            "airmask", File(output=False, optional=True, desc=airmask_desc)
        )
        self.add_trait(
            "artmask", File(output=False, optional=True, desc=artmask_desc)
        )
        self.add_trait(
            "headmask", File(output=False, optional=True, desc=headmask_desc)
        )
        self.add_trait(
            "rotmask", File(output=False, optional=True, desc=rotmask_desc)
        )
        self.add_trait(
            "hatmask", File(output=False, optional=True, desc=hatmask_desc)
        )
        self.add_trait(
            "segmentation",
            File(output=False, optional=True, desc=segmentation_desc),
        )
        self.add_trait(
            "in_inu", File(output=False, optional=True, desc=in_inu_desc)
        )
        self.add_trait(
            "in_noinu", File(output=False, optional=True, desc=in_noinu_desc)
        )
        self.add_trait(
            "pvms",
            traits.List(File(), output=False, optional=True, desc=pvms_desc),
        )
        self.add_trait(
            "mni_tpms",
            traits.List(
                File(), output=False, optional=True, desc=mni_tpms_desc
            ),
        )
        self.add_trait(
            "in_fwhm", File(output=False, optional=True, desc=in_fwhm_desc)
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

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
        super(AnatIQMs, self).list_outputs()

        if self.in_ras:
            valid_ext, in_ext, fileName = checkFileExt(self.in_ras, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized ...!")
                return

            report_file = os.path.join(
                self.output_directory, (fileName + "_anat_qc.json")
            )

            if fileName:
                self.outputs["out_file"] = report_file

            else:
                print(
                    "- There was no output file deducted during "
                    "initialisation. Please check the input parameters...!"
                )

            # tags inheritance (optional)
            if self.outputs:
                self.inheritance_dict[self.outputs["out_file"]] = self.in_ras

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(AnatIQMs, self).run_process_mia()

        results_dict = {}

        imdata = nb.load(self.in_ras).get_fdata()

        # Try to load files
        has_in_noinu = True
        if self.in_noinu:
            try:
                imnii = nb.load(self.in_noinu)
            except (
                nb.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                has_in_noinu = False
                print("\nError with in_noinu file: ", e)
            else:
                # Load image corrected for INU
                inudata = np.nan_to_num(imnii.get_fdata())
                inudata[inudata < 0] = 0
        else:
            has_in_noinu = False

        has_segmentation = True
        if self.segmentation:
            try:
                # Load binary segmentation from FSL FAST
                segnii = nb.load(self.segmentation)
            except (
                nb.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                has_segmentation = False
                print("\nError with segmentation file: ", e)
                pass
            else:
                segdata = segnii.get_fdata().astype(np.uint8)
        else:
            has_segmentation = False

        has_airmask = True
        if self.airmask:
            try:
                # Load binary segmentation from FSL FAST
                airnii = nb.load(self.airmask)
            except (
                nb.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                has_airmask = False
                print("\nError with airmask file: ", e)
                pass
            else:
                airdata = airnii.get_fdata()
        else:
            has_airmask = False

        has_artmask = True
        if self.artmask:
            try:
                # Load binary segmentation from FSL FAST
                artnii = nb.load(self.artmask)
            except (
                nb.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                has_artmask = False
                print("\nError with artmask file: ", e)
                pass
            else:
                artdata = artnii.get_fdata().astype(np.uint8)
        else:
            has_artmask = False

        has_headmask = True
        if self.headmask:
            try:
                # Load binary segmentation from FSL FAST
                headnii = nb.load(self.headmask)
            except (
                nb.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                has_headmask = False
                print("\nError with headmask file: ", e)
                pass
            else:
                headdata = headnii.get_fdata().astype(np.uint8)
        else:
            has_headmask = False

        has_rotmask = True
        if self.rotmask:
            try:
                # Load binary segmentation from FSL FAST
                rotnii = nb.load(self.rotmask)
            except (
                nb.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                has_rotmask = False
                print("\nError with rotmask file: ", e)
                pass
            else:
                rotdata = rotnii.get_fdata().astype(np.uint8)
        else:
            has_rotmask = False

        has_pvms = True
        if self.rotmask:
            try:
                # Load Partial Volume Maps (pvms) from FSL FAST
                pvmniis = []
                for fname in self.pvms:
                    pvmniis.append(nb.load(fname))
            except (
                nb.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                has_pvms = False
                print("\nError with pvms files: ", e)
                pass
            else:
                # Load Partial Volume Maps (pvms) from FSL FAST
                pvmdata = []
                for pvmnii in pvmniis:
                    pvmdata.append(pvmnii.get_fdata(dtype="float32"))
                    if np.sum(pvmdata[-1] > 1e-4) < 10:
                        raise RuntimeError(
                            "Detected less than 10 voxels "
                            "belonging to one tissue prob. "
                            "map. MRIQC failed to process this "
                            "dataset."
                        )
        else:
            has_pvms = False

        has_mni_tpms = True
        if self.mni_tpms:
            try:
                mni_tpmsniis = []
                for fname in self.mni_tpms:
                    mni_tpmsniis.append(nb.load(fname))
            except (
                nb.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                has_mni_tpms = False
                print("\nError with mni_tpms files: ", e)
                pass
        else:
            has_mni_tpms = False

        has_stats = False
        if has_in_noinu and has_pvms and has_airmask:
            erode = np.all(
                np.array(imnii.header.get_zooms()[:3], dtype=np.float32) < 1.9
            )

            # Summary stats
            stats = summary_stats(inudata, pvmdata, airdata, erode=erode)
            has_stats = True
            results_dict["summary"] = stats

            # SNR
            snrvals = []
            results_dict["snr"] = {}
            for tlabel in ["csf", "wm", "gm"]:
                snrvals.append(
                    snr(
                        stats[tlabel]["median"],
                        stats[tlabel]["stdv"],
                        int(stats[tlabel]["n"]),
                    )
                )
                results_dict["snr"][tlabel] = snrvals[-1]
            results_dict["snr"]["total"] = float(np.mean(snrvals))

            snrvals = []
            results_dict["snrd"] = {
                tlabel: snr_dietrich(
                    stats[tlabel]["median"],
                    mad_air=stats["bg"]["mad"],
                    sigma_air=stats["bg"]["stdv"],
                )
                for tlabel in ["csf", "wm", "gm"]
            }
            results_dict["snrd"]["total"] = float(
                np.mean([val for _, val in list(results_dict["snrd"].items())])
            )

            # CNR
            results_dict["cnr"] = cnr(
                stats["wm"]["median"],
                stats["gm"]["median"],
                sqrt(sum(stats[k]["stdv"] ** 2 for k in ["bg", "gm", "wm"])),
            )

        if has_in_noinu and has_headmask and has_rotmask:
            # FBER
            results_dict["fber"] = fber(inudata, headdata, rotdata)

        if has_in_noinu and has_rotmask:
            # EFC
            results_dict["efc"] = efc(inudata, rotdata)

        if has_stats and has_in_noinu:
            # M2WM
            results_dict["wm2max"] = wm2max(inudata, stats["wm"]["median"])

        if has_airmask and has_artmask:
            # Artifacts
            results_dict["qi_1"] = art_qi1(airdata, artdata)

        if self.hatmask != Undefined:
            # Artifacts QI2
            results_dict["qi_2"], results_dict["histogram_qi2"] = art_qi2(
                imdata, nb.load(self.hatmask).get_fdata()
            )

        if has_stats:
            # CJV
            results_dict["cjv"] = cjv(
                # mu_wm, mu_gm, sigma_wm, sigma_gm
                stats["wm"]["median"],
                stats["gm"]["median"],
                stats["wm"]["mad"],
                stats["gm"]["mad"],
            )

        # FWHM
        if self.in_fwhm and has_in_noinu:
            try:
                f = open(self.in_fwhm)
                lines = f.readlines()
                str_fwhm = re.findall(r"[-+]?\d*\.*\d+", lines[0])
                fwhm = []
                for item in str_fwhm:
                    fwhm.append(float(item))
            except (FileNotFoundError, TypeError) as e:
                print("\nError with fwhm file: ", e)
                pass
            else:
                fwhm = np.array(fwhm[:3]) / np.array(
                    imnii.header.get_zooms()[:3]
                )
                results_dict["fwhm"] = {
                    "x": float(fwhm[0]),
                    "y": float(fwhm[1]),
                    "z": float(fwhm[2]),
                    "avg": float(np.average(fwhm)),
                }

        if has_pvms:
            # ICVs
            results_dict["icvs"] = volume_fraction(pvmdata)

        if has_pvms and has_segmentation:
            # RPVE
            results_dict["rpve"] = rpve(pvmdata, segdata)

        if has_in_noinu:
            # Image specs
            results_dict["size"] = {
                "x": int(inudata.shape[0]),
                "y": int(inudata.shape[1]),
                "z": int(inudata.shape[2]),
            }
            results_dict["spacing"] = {
                i: float(v)
                for i, v in zip(["x", "y", "z"], imnii.header.get_zooms()[:3])
            }

            try:
                results_dict["size"]["t"] = int(inudata.shape[3])
            except IndexError:
                pass

            try:
                results_dict["spacing"]["tr"] = float(
                    imnii.header.get_zooms()[3]
                )
            except IndexError:
                pass

        if has_segmentation:
            # Bias
            bias = nb.load(self.in_inu).get_fdata()[segdata > 0]
            results_dict["inu"] = {
                "range": float(
                    np.abs(
                        np.percentile(bias, 95.0) - np.percentile(bias, 5.0)
                    )
                ),
                "med": float(np.median(bias)),
            }  # pylint: disable=E1101

        if has_mni_tpms and has_pvms:
            mni_tpms = [tpm.get_fdata() for tpm in mni_tpmsniis]
            in_tpms = [pvm.get_fdata() for pvm in pvmniis]
            overlap = fuzzy_jaccard(in_tpms, mni_tpms)
            results_dict["tpm_overlap"] = {
                "csf": overlap[0],
                "gm": overlap[1],
                "wm": overlap[2],
            }

        # Flatten the dictionary
        flat_results_dict = _flatten_dict(results_dict)

        # Save results as json
        _, file_name = os.path.split(self.in_ras)
        file_name_no_ext, file_extension = os.path.splitext(file_name)
        if file_extension == ".gz":
            (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                file_name_no_ext
            )
            if file_extension_2 == ".nii":
                file_name_no_ext = file_name_no_ext_2

        report_file = os.path.join(
            self.output_directory, (file_name_no_ext + "_anat_qc.json")
        )

        with open(report_file, "w") as fp:
            json.dump(flat_results_dict, fp)


class BoldIQMs(ProcessMIA):
    """
    *Computes the functional IQMs*

    Please, see the complete documentation for the `BoldIQMs brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/BoldIQMs.html>`_

    adapted from `mriqc
    <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/functional.py#L243>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(BoldIQMs, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_epi_desc = (
            "EPI input image (a pathlike object or string "
            "representing a file)."
        )
        in_hmc_desc = (
            "Motion corrected input image (a pathlike object or"
            "string representing a file)."
        )
        in_tsnr_desc = (
            "tSNR input volume (a pathlike object or string "
            "representing a file)."
        )
        in_mask_desc = (
            "Input mask image (a pathlike object or string "
            "representing a file)."
        )
        in_fd_thresh_desc = "Motion threshold for FD computation (a float)"
        in_outliers_file_desc = (
            "outliers file (a pathlike object or string "
            "representing a file)."
        )
        in_QI_file_desc = (
            "Quality index file (a pathlike object or string "
            "representing a file)."
        )
        in_fwhm_file_desc = (
            "FWHM file (a pathlike object or string " "representing a file)."
        )
        in_dvars_file_desc = (
            "DVARS file (a pathlike object or string " "representing a file)."
        )
        in_fd_file_desc = (
            "FD file (a pathlike object or string " "representing a file)."
        )
        in_spikes_file_desc = (
            "Spikes file (a pathlike object or string " "representing a file)."
        )
        in_gcor_desc = "Global correlation value (a float)"
        in_dummy_TRs_desc = "Number of dummy scans (an int)"

        # Outputs description
        out_file_desc = "a json file containing IQMs"

        # Inputs traits
        self.add_trait(
            "in_epi", File(output=False, optional=False, desc=in_epi_desc)
        )
        self.add_trait(
            "in_hmc", File(output=False, optional=True, desc=in_hmc_desc)
        )
        self.add_trait(
            "in_tsnr", File(output=False, optional=True, desc=in_tsnr_desc)
        )
        self.add_trait(
            "in_mask", File(output=False, optional=True, desc=in_mask_desc)
        )
        self.add_trait(
            "in_outliers_file",
            File(output=False, optional=True, desc=in_outliers_file_desc),
        )
        self.add_trait(
            "in_QI_file",
            File(output=False, optional=True, desc=in_QI_file_desc),
        )
        self.add_trait(
            "in_fwhm_file",
            File(output=False, optional=True, desc=in_fwhm_file_desc),
        )
        self.add_trait(
            "in_dvars_file",
            File(output=False, optional=True, desc=in_dvars_file_desc),
        )
        self.add_trait(
            "in_fd_file",
            File(output=False, optional=True, desc=in_fd_file_desc),
        )
        self.add_trait(
            "in_fd_thresh",
            traits.Float(
                0.2, output=False, optional=True, desc=in_fd_thresh_desc
            ),
        )
        self.add_trait(
            "in_spikes_file",
            File(output=False, optional=True, desc=in_spikes_file_desc),
        )
        self.add_trait(
            "in_gcor",
            traits.Float(output=False, optional=True, desc=in_gcor_desc),
        )
        self.add_trait(
            "in_dummy_TRs",
            traits.Int(output=False, optional=True, desc=in_dummy_TRs_desc),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

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
        super(BoldIQMs, self).list_outputs()

        if self.in_epi:
            valid_ext, in_ext, fileName = checkFileExt(self.in_epi, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized ...!")
                return

            report_file = os.path.join(
                self.output_directory, (fileName + "_bold_qc.json")
            )

            if fileName:
                self.outputs["out_file"] = report_file

            else:
                print(
                    "- There was no output file deducted during "
                    "initialisation. Please check the input parameters...!"
                )

            # tags inheritance (optional)
            if self.outputs:
                self.inheritance_dict[self.outputs["out_file"]] = self.in_epi

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(BoldIQMs, self).run_process_mia()

        # Get the mean EPI data and get it ready
        epinii = nb.load(self.in_epi)
        epidata = np.nan_to_num(epinii.get_fdata())
        epidata = epidata.astype(np.float32)
        epidata[epidata < 0] = 0

        # Get EPI data (with mc done) and get it ready
        has_in_hmc = True
        if self.in_hmc:
            try:
                hmcnii = nb.load(self.in_hmc)
            except (
                nb.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                has_in_hmc = False
                print("\nError with in_hmc file: ", e)
            else:
                hmcdata = np.nan_to_num(hmcnii.get_fdata())
                hmcdata = hmcdata.astype(np.float32)
                hmcdata[hmcdata < 0] = 0
        else:
            has_in_hmc = False

        # Get EPI data (with mc done) and get it ready
        has_in_mask = True
        if self.in_mask:
            try:
                msknii = nb.load(self.in_mask)
            except (
                nb.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                has_in_mask = False
                print("\nError with in_mask file: ", e)
            else:
                mskdata = np.nan_to_num(msknii.get_fdata())
                mskdata[mskdata < 0] = 0
                mskdata[mskdata > 0] = 1
                mskdata = mskdata.astype(np.uint8)
        else:
            has_in_mask = False

        # Get tsnr data and get it ready
        has_in_tsnr = True
        if self.in_tsnr:
            try:
                tsnr_nii = nb.load(self.in_tsnr)
            except (
                nb.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                has_in_tsnr = False
                print("\nError with in_mask file: ", e)
            else:
                # nibabel: get_data() is deprecated in favour of get_fdata().
                # We only keep here get_data() to achieve reproducibility
                # with native mriqc results (which use get_data() in V22.06)
                tsnr_data = tsnr_nii.get_data()
        else:
            has_in_tsnr = False

        results_dict = {}

        if has_in_mask:
            # Summary stats
            stats = summary_stats(epidata, mskdata, erode=True)
            results_dict["summary"] = stats

            # SNR
            results_dict["snr"] = snr(
                stats["fg"]["median"], stats["fg"]["stdv"], stats["fg"]["n"]
            )
            # FBER
            results_dict["fber"] = fber(epidata, mskdata)

        # EFC
        results_dict["efc"] = efc(epidata)

        if has_in_mask:
            results_dict["gsr"] = {}
            epidir = ["x", "y"]
            for axis in epidir:
                results_dict["gsr"][axis] = gsr(
                    epidata, mskdata, direction=axis
                )

        # aor
        if self.in_outliers_file:
            try:
                outliers = np.loadtxt(self.in_outliers_file, skiprows=0)
            except (FileNotFoundError, TypeError) as e:
                print("\nError with aor file: ", e)
                pass
            else:
                results_dict["vec_outliers"] = list(outliers)
                results_dict["aor"] = outliers.mean()

        # aqi
        if self.in_QI_file:
            try:
                with open(self.in_QI_file, "r") as fin:
                    lines = fin.readlines()
            except (FileNotFoundError, TypeError) as e:
                print("\nError with aqi file: ", e)
                pass
            else:
                results_dict["aqi"] = np.mean(
                    [
                        float(line.strip())
                        for line in lines
                        if not line.startswith("++")
                    ]
                )

        # DVARS
        if self.in_dvars_file:
            try:
                dvars = np.loadtxt(
                    self.in_dvars_file, skiprows=1, usecols=list(range(3))
                )
            except (FileNotFoundError, TypeError) as e:
                print("\nError with dvars file: ", e)
                pass
            else:
                results_dict["vec_dvars"] = list(dvars[:, 1])
                dvars_avg = dvars.mean(axis=0)
                dvars_col = ["std", "nstd", "vstd"]
                results_dict["dvars"] = {
                    key: float(val) for key, val in zip(dvars_col, dvars_avg)
                }

        # tSNR
        if has_in_tsnr:
            results_dict["tsnr"] = float(np.median(tsnr_data[mskdata > 0]))

        # FD
        if self.in_fd_file:
            try:
                fd_data = np.loadtxt(self.in_fd_file, skiprows=1)
            except (FileNotFoundError, TypeError) as e:
                print("\nError with fd file: ", e)
                pass
            else:
                results_dict["vec_fd"] = list(fd_data)
                if self.in_fd_thresh:
                    num_fd = np.float((fd_data > self.in_fd_thresh).sum())
                    results_dict["fd"] = {
                        "mean": float(fd_data.mean()),
                        "num": int(num_fd),
                        "perc": float(num_fd * 100 / (len(fd_data) + 1)),
                    }
                else:
                    results_dict["fd"] = {"mean": float(fd_data.mean())}

        # FWHM
        if self.in_fwhm_file and has_in_hmc:
            try:
                f = open(self.in_fwhm_file)
                lines = f.readlines()
                str_fwhm = re.findall(r"[-+]?\d*\.*\d+", lines[0])
                fwhm = []
                for item in str_fwhm:
                    fwhm.append(float(item))
            except (FileNotFoundError, TypeError) as e:
                print("\nError with fwhm file: ", e)
                pass
            else:
                fwhm = np.array(fwhm[:3]) / np.array(
                    hmcnii.header.get_zooms()[:3]
                )
                results_dict["fwhm"] = {
                    "x": float(fwhm[0]),
                    "y": float(fwhm[1]),
                    "z": float(fwhm[2]),
                    "avg": float(np.average(fwhm)),
                }

        if has_in_hmc:
            # Image specs
            results_dict["size"] = {
                "x": int(hmcdata.shape[0]),
                "y": int(hmcdata.shape[1]),
                "z": int(hmcdata.shape[2]),
            }
            results_dict["spacing"] = {
                i: float(v)
                for i, v in zip(["x", "y", "z"], hmcnii.header.get_zooms()[:3])
            }

            try:
                results_dict["size"]["t"] = int(hmcdata.shape[3])
            except IndexError:
                pass

            try:
                results_dict["spacing"]["tr"] = float(
                    hmcnii.header.get_zooms()[3]
                )
            except IndexError:
                pass

        # GCOR
        try:
            gcor = float(self.in_gcor)
        except ValueError:
            print("\ngcor value error")
        else:
            results_dict["gcor"] = gcor

        # Dummy TRs
        try:
            dummy_trs = int(self.in_dummy_TRs)
        except ValueError:
            print("\ndummy_trs value error")
        else:
            results_dict["dummy_trs"] = dummy_trs

        # Spikes
        if self.in_spikes_file:
            try:
                spikes_data = np.loadtxt(self.in_spikes_file, skiprows=0)
            except (FileNotFoundError, TypeError) as e:
                print("\nError with fd file: ", e)
                pass
            else:
                for i in range(0, len(spikes_data[:, 0])):
                    spike_name = "vec_spikes_" + str(i)
                    results_dict[spike_name] = list(spikes_data[i, :])

        # Flatten the dictionary
        flat_results_dict = _flatten_dict(results_dict)

        # Save results as json
        _, file_name = os.path.split(self.in_epi)
        file_name_no_ext, file_extension = os.path.splitext(file_name)
        if file_extension == ".gz":
            (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                file_name_no_ext
            )
            if file_extension_2 == ".nii":
                file_name_no_ext = file_name_no_ext_2

        report_file = os.path.join(
            self.output_directory, (file_name_no_ext + "_bold_qc.json")
        )

        with open(report_file, "w") as fp:
            json.dump(flat_results_dict, fp)


class BoldIQMsPlot(ProcessMIA):
    """
    Plot a figure showing the slice-wise signal intensity at the extremes for
    the identification of spikes,
    the outliers metric, the DVARS, the FD and the carpetplot.
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(BoldIQMsPlot, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_func_desc = "Bold image"
        in_outliers_file_desc = "Outliers file."
        in_dvars_file_desc = "DVARS file."
        in_fd_file_desc = "FD file."
        fd_thresh_desc = "Motion threshold for FD computation."
        in_spikes_file_desc = " Spikes file."
        carpet_seg_desc = "Carpet segmentation."
        drop_trs_desc = "Number of dummy scans drops."
        tr_desc = "Repetition time"

        # Outputs description
        out_file_desc = (
            "A figure with carpet and outliers/dvars/FD/spikes plot"
        )

        # Inputs traits
        self.add_trait(
            "in_func", File(output=False, optional=False, desc=in_func_desc)
        )

        self.add_trait(
            "in_outliers_file",
            File(output=False, optional=True, desc=in_outliers_file_desc),
        )

        self.add_trait(
            "in_dvars_file",
            File(output=False, optional=True, desc=in_dvars_file_desc),
        )

        self.add_trait(
            "in_fd_file",
            File(output=False, optional=True, desc=in_fd_file_desc),
        )

        self.add_trait(
            "in_spikes_file",
            File(output=False, optional=True, desc=in_spikes_file_desc),
        )

        self.add_trait(
            "carpet_seg",
            File(output=False, optional=True, desc=carpet_seg_desc),
        )

        self.add_trait(
            "drop_trs",
            traits.Int(0, output=False, optional=True, desc=drop_trs_desc),
        )

        self.add_trait(
            "fd_thresh",
            traits.Float(
                0.2, output=False, optional=True, desc=fd_thresh_desc
            ),
        )

        self.add_trait(
            "tr",
            traits.Either(
                Undefined,
                traits.Float(),
                output=False,
                optional=True,
                desc=tr_desc,
            ),
        )
        self.tr = Undefined

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=True, desc=out_file_desc)
        )

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key 'notInDb' of the
        self.outputs dictionary. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(BoldIQMsPlot, self).list_outputs()

        if self.tr is Undefined:
            # Get TR in the database
            tr = get_dbFieldValue(self.project, self.in_func, "RepetitionTime")
            if tr is None:
                print("Please add TR in database for functional file")
                return self.make_initResult()
            if tr:
                self.tr = tr[0]

        # Outputs definition and tags inheritance (optional)
        if self.in_func:
            valid_ext, in_ext, fileName = checkFileExt(self.in_func, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized ...!")
                return

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_fmriplot.png"
                )

            else:
                print("No output_directory was found...!\n")
                return

            self.inheritance_dict[self.outputs["out_file"]] = self.in_func

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(BoldIQMsPlot, self).run_process_mia()

        dataframe = pd.DataFrame(
            {
                "outliers": np.loadtxt(
                    self.in_outliers_file, usecols=[0]
                ).tolist(),
                # Pick non-standardize dvars (col 1)
                # First timepoint is NaN (difference)
                "DVARS": [np.nan]
                + np.loadtxt(
                    self.in_dvars_file, skiprows=1, usecols=[1]
                ).tolist(),
                # First timepoint is zero (reference volume)
                "FD": [0.0]
                + np.loadtxt(
                    self.in_fd_file, skiprows=1, usecols=[0]
                ).tolist(),
            }
        )

        input_data = nb.load(self.in_func)
        seg_file = self.carpet_seg
        dataset, segments = _nifti_timeseries(input_data, seg_file)

        fig = fMRIPlot(
            dataset,
            segments=segments,
            spikes_files=([self.in_spikes_file]),
            tr=(self.tr) / 1000,
            confounds=dataframe,
            units={"outliers": "%", "FD": "mm"},
            vlines={"FD": [self.fd_thresh]},
            nskip=self.drop_trs,
        ).plot()
        fig.savefig(self.out_file, bbox_inches="tight")


class CarpetParcellation(ProcessMIA):
    """
    *Dilate the brain mask, remove it from itself then generate the union
    of the obtained crown mask and the EPI parcellation*

    Please, see the complete documentation for the `CarpetParcellation brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/CarpetParcellation.html>`_

    Adapted from
    `niworkflows binary dilatation
    <https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/morphology.py#L46>`_,
    `niworkflows binary subtraction
    <https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/morphology.py#L79>`_,
    `mriqc
    <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/functional.py#L1022>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(CarpetParcellation, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        segmentation_desc = (
            "EPI segmentation (a pathlike object or string "
            "representing a file)."
        )
        brainmask_desc = (
            "Brain mask (a pathlike object or string " "representing a file)."
        )

        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "segmentation filename of the output file (a string)."
        )

        # Outputs description
        out_file_desc = (
            "The output file (a pathlike object or a "
            "string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "segmentation",
            File(output=False, optional=False, desc=segmentation_desc),
        )

        self.add_trait(
            "brainmask",
            File(output=False, optional=False, desc=brainmask_desc),
        )

        self.add_trait(
            "out_prefix",
            traits.String(
                "cseg_", output=False, optional=True, desc=out_prefix_desc
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=True, desc=out_file_desc)
        )

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key 'notInDb' of the
        self.outputs dictionary. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(CarpetParcellation, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.segmentation and self.brainmask:
            if self.out_prefix == Undefined:
                self.out_prefix = "cseg_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "cseg" ...'
                )

            if self.output_directory:
                # brainmask name used for carpetplot base name
                # in order to get an unique path for each subject
                valid_ext, in_ext, fileName = checkFileExt(self.brainmask, EXT)

                if not valid_ext:
                    print("\nThe input image format is not recognized ...!")
                    return

                self.outputs["out_file"] = os.path.join(
                    self.output_directory,
                    self.out_prefix + fileName + "." + in_ext,
                )

            else:
                print("No output_directory was found...!\n")
                return

            self.inheritance_dict[self.outputs["out_file"]] = self.brainmask

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(CarpetParcellation, self).run_process_mia()

        # Binary dilation
        brainmask_img = nb.load(self.brainmask)
        brainmaskdata = np.bool_(brainmask_img.dataobj)

        # Obtain dilated brain mask
        dilated = image_binary_dilation(
            brainmaskdata,
            radius=2,
        )

        # Binary subtraction
        subtraction = np.bool_(dilated)
        subtraction[np.bool_(brainmaskdata)] = False

        # Carpet parcellation
        img = nb.load(self.segmentation)

        lut = np.zeros((256,), dtype="uint8")
        lut[100:201] = 1  # Ctx GM
        lut[30:99] = 2  # dGM
        lut[1:11] = 3  # WM+CSF
        lut[255] = 4  # Cerebellum
        # Apply lookup table
        seg = lut[np.asanyarray(img.dataobj, dtype="uint16")]
        seg[np.asanyarray(subtraction, dtype=int) > 0] = 5

        # Out file name
        _, file_name = os.path.split(self.brainmask)

        file_out = os.path.join(
            self.output_directory, (self.out_prefix + file_name)
        )

        outimg = img.__class__(seg.astype("uint8"), img.affine, img.header)
        outimg.set_data_dtype("uint8")
        outimg.to_filename(file_out)


class ComputeDVARS(ProcessMIA):
    """
    *Computes the DVARS*

    Please, see the complete documentation for the `ComputeDVARS brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/ComputeDVARS.html>`_

    adapted from `nipype
    <https://github.com/nipy/nipype/blob/f662acfce8def4717e0c3414618f3a5de5913b31/nipype/algorithms/confounds.py#L100>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ComputeDVARS, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = (
            "Bold image after HMC (a pathlike object or string "
            "representing a file)."
        )
        in_mask_desc = (
            "Brain mask (a pathlike object or string " "representing a file)."
        )
        remove_zero_variance_desc = (
            "Remove voxels with zero variance" "(a bool)."
        )
        intensity_normalization_desc = (
            "Divide value in each voxel at each"
            "timepoint by the median calculated"
            "across all voxels and timepoints"
            "within the mask (if specified) and"
            "then multiply by the value specified"
            "by this parameter."
            "By using the default (1000)"
            "output DVARS will be expressed in "
            "x10 % BOLD units compatible"
            "with Power et al.2012."
            "Set this to 0 to disable intensity"
            "normalization altogether. (a float)"
        )
        variance_tol_desc = (
            'Maximum variance to consider "close to" zero'
            "for the purposes of removal"
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filename of the output file (a string)."
        )

        # Outputs description
        out_file_desc = (
            "The output file (a pathlike object or a "
            "string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "in_mask", File(output=False, optional=False, desc=in_mask_desc)
        )

        self.add_trait(
            "remove_zero_variance",
            traits.Bool(
                True,
                optional=True,
                output=False,
                desc=remove_zero_variance_desc,
            ),
        )

        self.add_trait(
            "intensity_normalization",
            traits.Float(
                1000.0,
                optional=True,
                output=False,
                desc=intensity_normalization_desc,
            ),
        )

        self.add_trait(
            "variance_tol",
            traits.Float(
                0.0000001000,
                optional=True,
                output=False,
                desc=variance_tol_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            traits.String(
                "dvars_", output=False, optional=True, desc=out_prefix_desc
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key 'notInDb' of the
        self.outputs dictionary. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(ComputeDVARS, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file and self.in_mask:
            if self.out_prefix == Undefined:
                self.out_prefix = "dvars_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "dvars" ...'
                )

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is not recognized ...!")
                    return

                self.outputs["out_file"] = os.path.join(
                    self.output_directory, self.out_prefix + fileName + ".out"
                )

            else:
                print("No output_directory was found...!\n")
                return

            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ComputeDVARS, self).run_process_mia()

        # Out file name
        _, file_name = os.path.split(self.in_file)
        file_name_no_ext, file_extension = os.path.splitext(file_name)
        if file_extension == ".gz":
            (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                file_name_no_ext
            )
            if file_extension_2 == ".nii":
                file_name_no_ext = file_name_no_ext_2

        file_out = os.path.join(
            self.output_directory,
            (self.out_prefix + file_name_no_ext + ".out"),
        )

        # Load data
        func = np.float32(nb.load(self.in_file).dataobj)
        mask = np.bool_(nb.load(self.in_mask).dataobj)

        if len(func.shape) != 4:
            raise RuntimeError("Input fMRI dataset should be 4-dimensional")

        mfunc = func[mask]

        if self.intensity_normalization != 0:
            mfunc = (mfunc / np.median(mfunc)) * self.intensity_normalization

        # Robust standard deviation (we are using "lower" interpolation
        # because this is what FSL is doing
        func_sd = (
            np.percentile(mfunc, 75, axis=1, interpolation="lower")
            - np.percentile(mfunc, 25, axis=1, interpolation="lower")
        ) / 1.349

        if self.remove_zero_variance:
            zero_variance_voxels = func_sd > self.variance_tol
            mfunc = mfunc[zero_variance_voxels, :]
            func_sd = func_sd[zero_variance_voxels]

        # Compute (non-robust) estimate of lag-1 autocorrelation
        ar1 = np.apply_along_axis(
            _AR_est_YW,
            1,
            regress_poly(0, mfunc, remove_mean=True)[0].astype(np.float32),
            1,
        )[:, 0]

        # Compute (predicted) standard deviation of temporal
        # difference time series
        diff_sdhat = np.squeeze(np.sqrt(((1 - ar1) * 2).tolist())) * func_sd
        diff_sd_mean = diff_sdhat.mean()

        # Compute temporal difference time series
        func_diff = np.diff(mfunc, axis=1)

        # DVARS (no standardization)
        dvars_nstd = np.sqrt(np.square(func_diff).mean(axis=0))

        # standardization
        dvars_stdz = dvars_nstd / diff_sd_mean

        try:
            # voxelwise standardization
            diff_vx_stdz = np.square(
                func_diff / np.array([diff_sdhat] * func_diff.shape[-1]).T
            )
            dvars_vx_stdz = np.sqrt(diff_vx_stdz.mean(axis=0))
        except Exception:
            print("\nError calculating vx-wise std DVARS...!")
            np.savetxt(
                file_out,
                np.vstack((dvars_stdz, dvars_nstd)).T,
                fmt=b"%0.8f",
                delimiter=b"\t",
                header="std DVARS\tnon-std DVARS",
                comments="",
            )
        else:
            np.savetxt(
                file_out,
                np.vstack((dvars_stdz, dvars_nstd, dvars_vx_stdz)).T,
                fmt=b"%0.8f",
                delimiter=b"\t",
                header="std DVARS\tnon-std DVARS\tvx-wise std DVARS",
                comments="",
            )


class FramewiseDisplacement(ProcessMIA):
    """
    *Calculate the FD (framewise displacement) as in [Power2012]*

    This implementation reproduces the calculation in fsl_motion_outliers

    `[Power2012] <https://doi.org/10.1016/j.neuroimage.2011.10.018>`_ Power
    et al., Spurious but systematic correlations in functional connectivity
    MRI networks arise from subject motion, NeuroImage 59(3).

    Please, see the complete documentation for the `FramewiseDisplacement
    brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/FramewiseDisplacement.html>`_

    adapted from `nipype
    <https://github.com/nipy/nipype/blob/f662acfce8def4717e0c3414618f3a5de5913b31/nipype/algorithms/confounds.py#L298>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(FramewiseDisplacement, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = (
            "Motion parameters file (a pathlike object or string "
            "representing a file)."
        )
        parameter_source_desc = (
            "Source of movement parameters"
            "(a string which is FSL "
            "or AFNI or SPM or FSFAST or NIPY"
        )
        radius_desc = "Radius in mm to calculate angular FDs (a float)."
        normalize_desc = "Calculate FD in mm/s (a bool)."
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filename of the output file (a string)."
        )

        # Outputs description
        out_file_desc = (
            "The output file (a pathlike object or a "
            "string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "parameter_source",
            traits.Enum(
                "FSL",
                "AFNI",
                "SPM",
                "FSFAST",
                "NIPY",
                output=False,
                optional=True,
                desc=parameter_source_desc,
            ),
        )

        self.add_trait(
            "radius",
            traits.Float(50.0, optional=True, output=False, desc=radius_desc),
        )

        self.add_trait(
            "normalize",
            traits.Bool(
                False, optional=True, output=False, desc=normalize_desc
            ),
        )

        self.add_trait(
            "out_prefix",
            traits.String(
                "fd_", output=False, optional=True, desc=out_prefix_desc
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key 'notInDb' of the
        self.outputs dictionary. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(FramewiseDisplacement, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.out_prefix == Undefined:
                self.out_prefix = "fd_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "fd" ...'
                )

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(
                    self.in_file, {"TXT": "txt"}
                )

                if not valid_ext:
                    print("\nThe input image format is not recognized ...!")
                    return

                self.outputs["out_file"] = os.path.join(
                    self.output_directory, self.out_prefix + fileName + ".out"
                )

            else:
                print("No output_directory was found...!\n")
                return

            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(FramewiseDisplacement, self).run_process_mia()

        mpars = np.loadtxt(self.in_file)  # mpars is N_t x 6
        mpars = np.apply_along_axis(
            func1d=normalize_mc_params,
            axis=1,
            arr=mpars,
            source=self.parameter_source,
        )
        diff = mpars[:-1, :6] - mpars[1:, :6]
        diff[:, 3:6] *= self.radius
        fd_res = np.abs(diff).sum(axis=1)

        # Image save
        _, file_name = os.path.split(self.in_file)
        file_name_no_ext, file_extension = os.path.splitext(file_name)
        if file_extension == ".gz":
            (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                file_name_no_ext
            )
            if file_extension_2 == ".nii":
                file_name_no_ext = file_name_no_ext_2

        file_out = os.path.join(
            self.output_directory,
            (self.out_prefix + file_name_no_ext + ".out"),
        )

        np.savetxt(
            file_out, fd_res, header="FramewiseDisplacement", comments=""
        )


class Mean_stdDev_calc(ProcessMIA):
    """
    *Makes the mean and standard deviation of parametric_maps*

    - The parametric_maps are first resized, if necessary, to the size of the
      rois_files. Next, the parametric_maps and the rois_files are convolved.
      Finally, the mean and standard deviation are calculated for the
      corresponding ROIs.
    - The “PatientName_data/ROI_data/ROI_analysis” directory is created to
      receive the results. If this directory exists at runtime, it is
      overwritten.
    - To work correctly, the database entry for the first element of
      parametric_maps must have the PatientName tag filled in.

    Please, see the complete documentation for the `Mean_stdDev_calc
    brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/Mean_stdDev_calc.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Mean_stdDev_calc, self).__init__()

        # Inputs description
        parametric_maps_desc = "A list of files (existing, uncompressed file)"
        rois_files_desc = (
            "A list of regions of interest applied to the parametric maps to "
            "calculate their mean and standard deviation"
        )
        contrast_type_desc = "The Contrast used (a string, ex. BOLD)"
        prefix_to_delete_desc = (
            "The prefix to delete from the " "deduced ROI name (a string)"
        )

        # Outputs description
        mean_out_files_desc = (
            "A list of .txt files with the calculated average for each "
            "ROI after convolution with parametric_maps_desc"
        )
        std_out_files_desc = (
            "A list of .txt files with the standard deviation for each ROI "
            "after convolution with parametric_maps_desc"
        )

        # Inputs traits
        self.add_trait(
            "parametric_maps",
            traits.List(
                traits.File(exists=True),
                output=False,
                optional=False,
                desc=parametric_maps_desc,
            ),
        )
        self.parametric_maps = Undefined

        self.add_trait(
            "rois_files",
            InputMultiPath(
                ImageFileSPM(),
                output=False,
                optional=False,
                desc=rois_files_desc,
            ),
        )

        self.add_trait(
            "contrast_type",
            traits.String(
                "BOLD", output=False, optional=True, desc=contrast_type_desc
            ),
        )

        self.add_trait(
            "prefix_to_delete",
            traits.String(
                "conv", output=False, optional=True, desc=prefix_to_delete_desc
            ),
        )

        # Output traits
        self.add_trait(
            "mean_out_files",
            OutputMultiPath(File(), output=True, desc=mean_out_files_desc),
        )

        self.add_trait(
            "std_out_files",
            OutputMultiPath(File(), output=True, desc=std_out_files_desc),
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
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key 'notInDb' of the
        self.outputs dictionary. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Mean_stdDev_calc, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.parametric_maps != Undefined and self.rois_files != Undefined:
            # FIXME: We retrieve the name of the patient from the first element
            #        of parametric_maps. This is only fine if all the elements
            #        of parametric_maps correspond to the same patient.
            patient_name = get_dbFieldValue(
                self.project, self.parametric_maps[0], "PatientName"
            )

            if patient_name is None:
                print(
                    "\nMean_stdDev_cal brick:\nThe PatientName tag is not "
                    "filled in the database for the {} file ...\n The "
                    "initialization is "
                    "aborted...".format(self.parametric_maps[0])
                )
                return self.make_initResult()

            self.dict4runtime["patient_name"] = patient_name
            analysis_dir = os.path.join(
                self.output_directory,
                patient_name + "_data",
                "ROI_data",
                "ROI_analysis",
            )
            mean_out_files = []
            std_out_files = []

            if self.contrast_type.isspace() or not self.contrast_type:
                self.contrast_type = "UnknownContrast"

            if (
                (not self.prefix_to_delete)
                or (self.prefix_to_delete.isspace())
                or (self.prefix_to_delete in [Undefined, "<undefined>"])
            ):
                self.prefix_to_delete = " "

            for parametric_map in self.parametric_maps:
                param_file_name = os.path.basename(parametric_map)

                if "_" not in param_file_name:
                    parameter = os.path.splitext(param_file_name)

                else:
                    parameter = param_file_name[: param_file_name.index("_")]

                for roi in self.rois_files:
                    roi_name, _ = os.path.splitext(os.path.basename(roi))

                    if self.prefix_to_delete != " " and roi_name.startswith(
                        self.prefix_to_delete
                    ):
                        # fmt: off
                        roi_name = roi_name[len(self.prefix_to_delete):]
                        # fmt: on

                    mean_out_file = os.path.join(
                        analysis_dir,
                        roi_name
                        + "_mean_"
                        + parameter
                        + "_"
                        + self.contrast_type
                        + ".txt",
                    )
                    # FIXME: In the latest version of mia, indexing of the
                    #        database with particular tags defined in the
                    #        processes is done only at the end of the
                    #        initialisation of the whole pipeline. So we
                    #        cannot use the value of these tags in other
                    #        processes of the pipeline at the time of
                    #        initialisation (see populse_mia #290). Until
                    #        better we use a quick and dirty hack with the
                    #        set_dbFieldValue() function !
                    tag_to_add = dict()
                    tag_to_add["name"] = "PatientName"
                    tag_to_add["field_type"] = "string"
                    tag_to_add["description"] = ""
                    tag_to_add["visibility"] = True
                    tag_to_add["origin"] = "user"
                    tag_to_add["unit"] = None
                    tag_to_add["default_value"] = None
                    tag_to_add["value"] = patient_name
                    set_dbFieldValue(self.project, mean_out_file, tag_to_add)
                    mean_out_files.append(mean_out_file)
                    self.inheritance_dict[mean_out_file] = parametric_map

                    std_out_file = os.path.join(
                        analysis_dir,
                        roi_name
                        + "_std_"
                        + parameter
                        + "_"
                        + self.contrast_type
                        + ".txt",
                    )
                    tag_to_add = dict()
                    tag_to_add["name"] = "PatientName"
                    tag_to_add["field_type"] = "string"
                    tag_to_add["description"] = ""
                    tag_to_add["visibility"] = True
                    tag_to_add["origin"] = "user"
                    tag_to_add["unit"] = None
                    tag_to_add["default_value"] = None
                    tag_to_add["value"] = patient_name
                    set_dbFieldValue(self.project, std_out_file, tag_to_add)
                    std_out_files.append(std_out_file)
                    self.inheritance_dict[std_out_file] = parametric_map

            if mean_out_files != []:
                self.outputs["mean_out_files"] = mean_out_files

            if std_out_files != []:
                self.outputs["std_out_files"] = std_out_files

        # FIXME: What are we doing about the tags inheritance?

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        # No need the next line (we don't use self.process and SPM)
        # super(Mean_stdDev_calc, self).run_process_mia()

        pat_name_dir = os.path.join(
            self.output_directory, self.dict4runtime["patient_name"] + "_data"
        )

        if not os.path.exists(pat_name_dir):
            os.mkdir(pat_name_dir)

        roi_data_dir = os.path.join(pat_name_dir, "ROI_data")

        if not os.path.exists(roi_data_dir):
            os.mkdir(roi_data_dir)

        analysis_dir = os.path.join(roi_data_dir, "ROI_analysis")

        tmp = "None"

        if os.path.isdir(analysis_dir):
            tmp = tempfile.mktemp(dir=os.path.dirname(roi_data_dir))
            os.mkdir(tmp)
            shutil.move(analysis_dir, os.path.join(tmp, "ROI_analysis"))
            print(
                '\nMean_stdDev_calc brick:\nA "{}" folder already exists, '
                "it will be overwritten by this new "
                "calculation...".format(analysis_dir)
            )
        os.mkdir(analysis_dir)

        if os.path.isdir(tmp):
            shutil.rmtree(tmp)

        for parametric_map in self.parametric_maps:
            # Reading parametric map
            map_img = nb.load(parametric_map)
            map_data = map_img.get_fdata()
            # Setting the NaN to 0 in parametric map
            map_data = np.nan_to_num(map_data)
            # Deduction of the parameter from the name of the parametric map
            param_file_name = os.path.basename(parametric_map)

            if "_" not in param_file_name:
                parameter = os.path.splitext(param_file_name)

            else:
                parameter = param_file_name[: param_file_name.index("_")]

            for roi_file in self.rois_files:
                roi_img = nb.load(roi_file)
                roi_data = roi_img.get_fdata()
                # Setting the NaN to 0 in ROI images
                roi_data = np.nan_to_num(roi_data)
                # Deduction of the ROI name from the roi_file
                roi_name, _ = os.path.splitext(os.path.basename(roi_file))

                if self.prefix_to_delete != " " and roi_name.startswith(
                    self.prefix_to_delete
                ):
                    # fmt: off
                    roi_name = roi_name[len(self.prefix_to_delete):]
                    # fmt: on

                # Making sure that the ROIs and parametric images are at the
                # same size
                if roi_data.shape[:3] != map_data.shape[:3]:
                    map_data_max = max(map_data.max(), -map_data.min())
                    final_map_data = (
                        resize(map_data / map_data_max, roi_data.shape[:3])
                        * map_data_max
                    )

                # Convolution of the parametric map with the ROI images
                roi_thresh = (roi_data > 0).astype(float)
                result = final_map_data * roi_thresh

                # Calculating mean and standard deviation
                if np.size(result[result.nonzero()]) == 0:
                    print(
                        "\nMean_stdDev_cal brick:\nWarning: No result found "
                        "after convolution of the {0} ROI and the {1} "
                        "data".format(roi_file, parametric_map)
                    )
                    mean_result = 0
                    std_result = 0

                else:
                    mean_result = result[result.nonzero()].mean()
                    std_result = result[result.nonzero()].std()

                for calculation in ["mean", "std"]:
                    out_file = os.path.join(
                        analysis_dir,
                        roi_name
                        + "_"
                        + calculation
                        + "_"
                        + parameter
                        + "_"
                        + self.contrast_type
                        + ".txt",
                    )

                    with open(out_file, "w") as f:
                        if calculation == "mean":
                            f.write("%.3f" % mean_result)

                        else:
                            f.write("%.3f" % std_result)


class Result_collector(ProcessMIA):
    """
    *Save a file.xlm with the data collection for a patient*

    - To work correctly, the database entry for the first element of
      parameter_files must have the "PatientName" tag filled in.
    - The “PatientName_data/results_aggregation” directory is created to
      receive the results. If this directory exists at runtime, new results
      can overwrite old results with the same name.
    - To work correctly, the name of each file in parameter_files must be
      exactly like this: roi_hemi_calcul_param_contrast.txt, where
        - roi: region of interest (ex. ACA)
        - hemi: hemisphere (ex. L)
        - calcul: type of calcul (ex. mean)
        - param: the parameter studied (ex. spmT)
        - contrast: the type of contrast/effect used (ex. BOLD)

    Please, see the complete documentation for the `Result_collector
    brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/Result_collector.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Result_collector, self).__init__()

        # Inputs description
        parameter_files_desc = (
            "A list of .txt files. Each file contains one (and "
            "only one) value. The name of each file is used to "
            "define this value. The name must exactly be like "
            "this: roi_hemi_calcul_param_contrast"
        )
        laterality_index_desc = "Calculates the laterality indexes (a boolean)"
        patient_info_desc = (
            "A dictionary whose keys/values correspond to "
            "information about the patient "
            "(e.g. {"
            "'PatientName': 'ablair','Pathology': 'ACMD', "
            "'Age': 64, 'Sex': 'M', 'MR': '3T', "
            "'Gas': 'BACTAL', 'GasAdmin': 'MASK'}"
        )

        # Outputs description
        out_files_desc = (
            "A list of .xml files containing a summary of the input parameters"
        )

        # Inputs traits
        self.add_trait(
            "parameter_files",
            traits.List(
                traits.File(exists=True),
                output=False,
                optional=False,
                desc=parameter_files_desc,
            ),
        )
        self.parameter_files = Undefined

        self.add_trait(
            "laterality_index",
            traits.Bool(
                default_value=False,
                output=False,
                optional=True,
                desc=laterality_index_desc,
            ),
        )

        self.add_trait(
            "patient_info",
            traits.Dict(output=False, optional=True, desc=patient_info_desc),
        )
        self.patient_info = dict(
            PatientName=Undefined,
            Pathology=Undefined,
            Age=Undefined,
            Sex=Undefined,
            MR=Undefined,
            Gas=Undefined,
            GasAdmin=Undefined,
        )

        # Outputs traits
        self.add_trait(
            "out_files",
            traits.List(traits.File(), output=True, desc=out_files_desc),
        )
        self.out_files = Undefined

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
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key 'notInDb' of the
        self.outputs dictionary. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Result_collector, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.parameter_files != Undefined:
            # FIXME 1: We retrieve the name of the patient from the first
            #          element of parameter_files. This is only fine if all the
            #          elements of parameter_files correspond to the same
            #          patient.
            # FIXME 2: The data should be anonymized and we should use
            #          PatientRef instead of PatientName !
            patient_name = get_dbFieldValue(
                self.project, self.parameter_files[0], "PatientName"
            )
            if patient_name is None:
                print(
                    "\nResult_collector brick:\nThe PatientName tag is not "
                    "filled in the database for the {} file ...\n The "
                    "initialization is "
                    "aborted...".format(self.parameter_files[0])
                )

                return self.make_initResult()

            self.dict4runtime["patient_name"] = patient_name

            aggreg_results_dir = os.path.join(
                self.output_directory,
                patient_name + "_data",
                "results_aggregation",
            )
            out_files = set()
            res = dict()

            for data in self.parameter_files:
                # FIXME 1: We need to protect the following command line. Case
                #          where there are not 5 objects returned.
                # FIXME 2: Ideally, we should also make sure that they are well
                #          roi, hemi, etc ... This is difficult to achieve for
                #          this last point ...
                # roi: region of interest (ex. ACA)
                # hemi: hemisphere (ex. L)
                # calcul: type of calcul (ex. mean)
                # param: the parameter studied (ex. spmT)
                # contrast: The type of contrast/effect used (ex. BOLD)

                try:
                    roi, hemi, calcul, param, contrast = os.path.splitext(
                        os.path.basename(data)
                    )[0].split("_")

                except Exception as e:
                    print(
                        "\nResult_collector brick: initialization stopped "
                        "due to the following issue:\n"
                        " {}\n".format(e)
                    )
                    return self.make_initResult()

                out_files.add(
                    os.path.join(
                        aggreg_results_dir,
                        "{0}_{1}_{2}.xls".format(contrast, calcul, param),
                    )
                )

                if contrast not in res:
                    res[contrast] = {}

                if param not in res[contrast]:
                    res[contrast][param] = {}

                if calcul not in res[contrast][param]:
                    res[contrast][param][calcul] = {}

                if (
                    self.laterality_index is True
                    and "IL_" + calcul not in res[contrast][param]
                ):
                    res[contrast][param]["IL_" + calcul] = {}

                if roi not in res[contrast][param][calcul]:
                    res[contrast][param][calcul][roi] = {}

                if (
                    self.laterality_index is True
                    and roi not in res[contrast][param]["IL_" + calcul]
                ):
                    res[contrast][param]["IL_" + calcul][roi] = {}

                if hemi in res[contrast][param][calcul][roi]:
                    print(
                        "\nResult_collector brick:\nThe data for "
                        "{0}-{1}-{2}-{3}-{4} in {5} already exists "
                        "in the final result. Overwriting with the new "
                        "data ...\n".format(
                            contrast, param, calcul, roi, hemi, data
                        )
                    )

                res[contrast][param][calcul][roi][hemi] = True

                if (
                    self.laterality_index is True
                    and "L" in res[contrast][param][calcul][roi]
                    and "R" in res[contrast][param][calcul][roi]
                ):
                    out_files.add(
                        os.path.join(
                            aggreg_results_dir,
                            "{0}_{1}_{2}.xls".format(
                                contrast, "IL_" + calcul, param
                            ),
                        )
                    )

            if out_files:
                self.outputs["out_files"] = list(out_files)

            # FIXME: the data should be anonymized and we should use PatientRef
            #        instead of PatientName !

            if (
                self.patient_info.get("PatientName") is None
                or self.patient_info["PatientName"] == Undefined
            ):
                patient_ref = get_dbFieldValue(
                    self.project, self.parameter_files[0], "PatientRef"
                )

                if patient_ref is None:
                    self.patient_info["PatientName"] = patient_name

                else:
                    self.patient_info["PatientName"] = patient_ref

            if (
                self.patient_info.get("Pathology") is None
                or self.patient_info["Pathology"] == Undefined
            ):
                pathology = get_dbFieldValue(
                    self.project, self.parameter_files[0], "Pathology"
                )

                if pathology is None:
                    self.patient_info["Pathology"] = "Undefined"

                else:
                    self.patient_info["Pathology"] = pathology

            if (
                self.patient_info.get("Age") is None
                or self.patient_info["Age"] == Undefined
            ):
                age = get_dbFieldValue(
                    self.project, self.parameter_files[0], "Age"
                )

                if age is None:
                    self.patient_info["Age"] = "Undefined"

                else:
                    self.patient_info["Age"] = age

            if (
                self.patient_info.get("Sex") is None
                or self.patient_info["Sex"] == Undefined
            ):
                sex = get_dbFieldValue(
                    self.project, self.parameter_files[0], "Sex"
                )

                if sex is None:
                    self.patient_info["Sex"] = "Undefined"

                else:
                    self.patient_info["Sex"] = sex

            if (
                self.patient_info.get("MR") is None
                or self.patient_info["MR"] == Undefined
            ):
                mr = get_dbFieldValue(
                    self.project, self.parameter_files[0], "MR"
                )

                if mr is None:
                    self.patient_info["MR"] = "Undefined"

                else:
                    self.patient_info["MR"] = mr

            if (
                self.patient_info.get("Gas") is None
                or self.patient_info["Gas"] == Undefined
            ):
                gas = get_dbFieldValue(
                    self.project, self.parameter_files[0], "Gas"
                )

                if gas is None:
                    self.patient_info["Gas"] = "Undefined"

                else:
                    self.patient_info["Gas"] = gas

            if (
                self.patient_info.get("GasAdmin") is None
                or self.patient_info["GasAdmin"] == Undefined
            ):
                gas_admin = get_dbFieldValue(
                    self.project, self.parameter_files[0], "GasAdmin"
                )

                if gas_admin is None:
                    self.patient_info["GasAdmin"] = "Undefined"

                else:
                    self.patient_info["GasAdmin"] = gas_admin

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        # No need the next line (we don't use self.process and SPM)
        # super(Result_collector, self).run_process_mia()

        pat_name_dir = os.path.join(
            self.output_directory, self.dict4runtime["patient_name"] + "_data"
        )

        if not os.path.exists(pat_name_dir):
            os.mkdir(pat_name_dir)

        aggreg_results_dir = os.path.join(pat_name_dir, "results_aggregation")

        if not os.path.exists(aggreg_results_dir):
            os.mkdir(aggreg_results_dir)

        else:
            print(
                "Result_collector brick warning: The {} directory exists "
                "before the start of the calculation. Its contents can be "
                "overwritten by the result of the current "
                "calculation ...!".format(aggreg_results_dir)
            )

        res = dict()

        for data in self.parameter_files:
            # roi: region of interest (ex. ACA)
            # hemi: hemisphere (ex. L)
            # calcul: type of calcul (ex. mean)
            # param: the parameter studied (ex. spmT)
            # contrast: The type of contrast/effect used (ex. BOLD)
            roi, hemi, calcul, param, contrast = os.path.splitext(
                os.path.basename(data)
            )[0].split("_")

            if contrast not in res:
                res[contrast] = {}

            if param not in res[contrast]:
                res[contrast][param] = {}

            if calcul not in res[contrast][param]:
                res[contrast][param][calcul] = {}

            if (
                self.laterality_index is True
                and "IL_" + calcul not in res[contrast][param]
            ):
                res[contrast][param]["IL_" + calcul] = {}

            if roi not in res[contrast][param][calcul]:
                res[contrast][param][calcul][roi] = {}

            if (
                self.laterality_index is True
                and roi not in res[contrast][param]["IL_" + calcul]
            ):
                res[contrast][param]["IL_" + calcul][roi] = {}

            try:
                with open(data, "r") as f_read:
                    data_val = float(f_read.read())

            except Exception:
                print(
                    "\nResult_collector brick:\nNo result "
                    "found in {} ...\n".format(data)
                )
                data_val = "Undefined"

            if hemi in res[contrast][param][calcul][roi]:
                print(
                    "\nResult_collector brick:\nThe data for {0}-{1}-{2}-{3}"
                    "-{4} in {5} already exists in the final result. "
                    "Overwriting with the new data ...\n".format(
                        contrast, param, calcul, roi, hemi, data
                    )
                )

            res[contrast][param][calcul][roi][hemi] = data_val

            if (
                self.laterality_index is True
                and "L" in res[contrast][param][calcul][roi]
                and "R" in res[contrast][param][calcul][roi]
            ):
                left = res[contrast][param][calcul][roi]["L"]
                right = res[contrast][param][calcul][roi]["R"]
                res[contrast][param]["IL_" + calcul][roi] = (left - right) / (
                    left + right
                )

        for contrast in res:
            for param in res[contrast]:
                for calcul in res[contrast][param]:
                    out_file = os.path.join(
                        aggreg_results_dir,
                        "{0}_{1}_{2}.xls".format(contrast, calcul, param),
                    )

                    with open(out_file, "w") as f:
                        f.write("{0}\t".format("subjects"))
                        f.write("{0}\t".format("patho"))
                        f.write("{0}\t".format("age"))
                        f.write("{0}\t".format("sex"))
                        f.write("{0}\t".format("MR"))
                        f.write("{0}\t".format("Gaz"))
                        f.write("{0}\t".format("Admin"))

                        if not calcul.startswith("IL_"):
                            for roi in res[contrast][param][calcul]:
                                for hemi in res[contrast][param][calcul][roi]:
                                    f.write(
                                        "{0}_{1}_{2}\t".format(
                                            param, roi, hemi
                                        )
                                    )

                        else:
                            for roi in res[contrast][param][calcul]:
                                f.write("{0}_{1}\t".format(param, roi))

                        f.write(
                            "\n{0}\t".format(self.patient_info["PatientName"])
                        )
                        f.write("{0}\t".format(self.patient_info["Pathology"]))
                        f.write("{0}\t".format(self.patient_info["Age"]))
                        f.write("{0}\t".format(self.patient_info["Sex"]))
                        f.write("{0}\t".format(self.patient_info["MR"]))
                        f.write("{0}\t".format(self.patient_info["Gas"]))
                        f.write("{0}\t".format(self.patient_info["GasAdmin"]))

                        if not calcul.startswith("IL_"):
                            for roi in res[contrast][param][calcul]:
                                for hemi in res[contrast][param][calcul][roi]:
                                    f.write(
                                        "{0}\t".format(
                                            res[contrast][param][calcul][roi][
                                                hemi
                                            ]
                                        )
                                    )

                        else:
                            for roi in res[contrast][param][calcul]:
                                f.write(
                                    "{0}\t".format(
                                        res[contrast][param][calcul][roi]
                                    )
                                )


class Spikes(ProcessMIA):
    """
    *Computes the number of spikes*

    Please, see the complete documentation for the `Spikes brick in the
    populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/reports/Spikes.html>`_

    adapted from `mriqc spikes_mask function
    <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/functional.py#L939>`_
    and
    `mriqc Spikes class
    <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/functional.py#L223>`_


    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Spikes, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = (
            "A bold file (a pathlike object or string " "representing a file)."
        )
        no_zscore_desc = "Do not zscore (a boolean)"
        detrend_desc = "Do detrend (a boolean)."
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the output image file "
            "(a string)."
        )
        spike_thresh_desc = (
            "z-score to call one timepoint of one axial"
            " slice a spike (a float)"
        )
        skip_frames_desc = (
            "number of frames to skip in the beginning "
            "of the time series (an int)"
        )

        # Outputs description
        out_file_desc = (
            "The output file (a pathlike object or a "
            "string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )
        self.add_trait(
            "no_zscore",
            traits.Bool(
                True, output=False, optional=True, desc=no_zscore_desc
            ),
        )
        self.add_trait(
            "detrend",
            traits.Bool(False, output=False, optional=True, desc=detrend_desc),
        )
        self.add_trait(
            "spike_thresh",
            traits.Float(
                6.0, output=False, optional=True, desc=spike_thresh_desc
            ),
        )
        self.add_trait(
            "skip_frames",
            traits.Int(0, output=False, optional=True, desc=skip_frames_desc),
        )
        self.add_trait(
            "out_prefix",
            traits.String(
                "spikes_", output=False, optional=True, desc=out_prefix_desc
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key 'notInDb' of the
        self.outputs dictionary. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Spikes, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.out_prefix == Undefined:
                self.out_prefix = "spikes_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "spikes" ...'
                )

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is not recognized ...!")
                    return

                self.outputs["out_file"] = os.path.join(
                    self.output_directory, self.out_prefix + fileName + ".out"
                )

            else:
                print("No output_directory was found...!\n")
                return

            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Spikes, self).run_process_mia()

        func_nii = nb.load(self.in_file)
        func_data = func_nii.get_fdata()
        func_shape = func_data.shape

        orientation = nb.aff2axcodes(func_nii.affine)

        new_mask_3d = np.zeros(func_nii.shape[:3]) == 1

        if orientation[0] in ["L", "R"]:
            new_mask_3d[0:2, :, :] = True
            new_mask_3d[-3:-1, :, :] = True
        else:
            new_mask_3d[:, 0:2, :] = True
            new_mask_3d[:, -3:-1, :] = True

        ntsteps = func_shape[-1]
        tr = func_nii.header.get_zooms()[-1]
        nskip = self.skip_frames

        if self.detrend:
            data = func_data.reshape(-1, ntsteps)
            clean_data = clean(data[:, nskip:].T, t_r=tr, standardize=False).T
            new_shape = (
                func_shape[0],
                func_shape[1],
                func_shape[2],
                clean_data.shape[-1],
            )
            func_data = np.zeros(func_shape)
            func_data[..., nskip:] = clean_data.reshape(new_shape)

        mask_data = new_mask_3d.astype(np.uint8)
        mask_data[..., :nskip] = 0
        mask_data = np.stack([mask_data] * ntsteps, axis=-1)

        mask_data[..., : self.skip_frames] = 1
        brain = np.ma.array(func_data, mask=(mask_data == 1))

        if self.no_zscore:
            ts_z = find_peaks(brain)
        else:
            _, ts_z = find_spikes(brain, self.spike_thresh)

        # File save
        _, file_name = os.path.split(self.in_file)
        file_name_no_ext, file_extension = os.path.splitext(file_name)
        if file_extension == ".gz":
            (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                file_name_no_ext
            )
            if file_extension_2 == ".nii":
                file_name_no_ext = file_name_no_ext_2

        file_out = os.path.join(
            self.output_directory,
            (self.out_prefix + file_name_no_ext + ".out"),
        )

        np.savetxt(file_out, ts_z)


def art_qi1(airmask, artmask):
    r"""
    Detect artifacts in the image using the method described in [Mortamet2009].
    Caculates :math:`\text{QI}_1`, as the proportion of voxels with intensity
    corrupted by artifacts normalized by the number of voxels in the
    background:
    .. math ::
        \text{QI}_1 = \frac{1}{N} \sum\limits_{x\in X_\text{art}} 1
    Lower values are better.
    :param numpy.ndarray airmask: input air mask, without artifacts
    :param numpy.ndarray artmask: input artifacts mask
    """

    if airmask.sum() < 1:
        return -1.0

    # Count the number of voxels that remain after the opening operation.
    # These are artifacts.
    return float(artmask.sum() / (airmask.sum() + artmask.sum()))


def art_qi2(img, airmask, min_voxels=int(1e3), max_voxels=int(3e5)):
    r"""
    Calculates :math:`\text{QI}_2`, based on the goodness-of-fit of a centered
    :math:`\chi^2` distribution onto the intensity distribution of
    non-artifactual background (within the "hat" mask):
    .. math ::
        \chi^2_n = \frac{2}{(\sigma \sqrt{2})^{2n} \,
                            (n - 1)!}x^{2n - 1}\, e^{-\frac{x}{2}}
    where :math:`n` is the number of coil elements.
    :param numpy.ndarray img: input data
    :param numpy.ndarray airmask: input air mask without artifacts
    """

    from scipy.stats import chi2
    from sklearn.neighbors import KernelDensity

    # S. Ogawa was born
    np.random.seed(1191935)

    data = img[airmask > 0]
    data = data[data > 0]

    if len(data) < min_voxels:
        return 0.0

    modelx = (
        data
        if len(data) < max_voxels
        else np.random.choice(data, size=max_voxels)
    )

    x_grid = np.linspace(0.0, np.percentile(data, 99), 1000)

    # Estimate data pdf with KDE on a random subsample
    kde_skl = KernelDensity(
        bandwidth=0.05 * np.percentile(data, 98), kernel="gaussian"
    ).fit(modelx[:, np.newaxis])
    kde = np.exp(kde_skl.score_samples(x_grid[:, np.newaxis]))

    # Find cutoff
    kdethi = np.argmax(kde[::-1] > kde.max() * 0.5)

    # Fit X^2
    param = chi2.fit(modelx[modelx < np.percentile(data, 95)], 32)
    chi_pdf = chi2.pdf(x_grid, *param[:-2], loc=param[-2], scale=param[-1])

    # Compute goodness-of-fit (gof)
    gof = float(np.abs(kde[-kdethi:] - chi_pdf[-kdethi:]).mean())

    hist_dict = {
        "x_grid": x_grid.tolist(),
        "ref_pdf": kde.tolist(),
        "fit_pdf": chi_pdf.tolist(),
        "ref_data": modelx.tolist(),
        "cutoff_idx": float(kdethi),
    }
    return gof, hist_dict


def cjv(mu_wm, mu_gm, sigma_wm, sigma_gm):
    r"""
    Calculate the :abbr:`CJV (coefficient of joint variation)`, a measure
    related to :abbr:`SNR (Signal-to-Noise Ratio)` and
    :abbr:`CNR (Contrast-to-Noise Ratio)` that is presented as a proxy for
    the :abbr:`INU (intensity non-uniformity)` artifact [Ganzetti2016]_.
    Lower is better.
    .. math::
        \text{CJV} = \frac{\sigma_\text{WM} +
                           \sigma_\text{GM}}{|\mu_\text{WM} -
                           \mu_\text{GM}|}.
    :param float mu_wm: mean of signal within white-matter mask.
    :param float mu_gm: mean of signal within gray-matter mask.
    :param float sigma_wm: standard deviation of signal within
                           white-matter mask.
    :param float sigma_gm: standard deviation of signal within
                           gray-matter mask.
    :return: the computed CJV
    """
    return float((sigma_wm + sigma_gm) / abs(mu_wm - mu_gm))


def cnr(mu_wm, mu_gm, sigma_air):
    r"""
    Calculate the :abbr:`CNR (Contrast-to-Noise Ratio)` [Magnota2006]_.
    Higher values are better.
    .. math::
        \text{CNR} = \frac{|\mu_\text{GM} -
                           \mu_\text{WM} |}{\sqrt{\sigma_B^2 +
                           \sigma_\text{WM}^2 +
                           \sigma_\text{GM}^2}},
    where :math:`\sigma_B` is the standard deviation of the noise distribution
    withinmthe air (background) mask.

    :param float mu_wm: mean of signal within white-matter mask.
    :param float mu_gm: mean of signal within gray-matter mask.
    :param float sigma_air: standard deviation of the air surrounding the
                            head ("hat" mask).

    :return: the computed CNR
    """
    return float(abs(mu_wm - mu_gm) / sigma_air)


def efc(img, framemask=None):
    r"""
    Calculate the :abbr:`EFC (Entropy Focus Criterion)` [Atkinson1997]_.
    Uses the Shannon entropy of voxel intensities as an indication of ghosting
    and blurring induced by head motion. A range of low values is better,
    with EFC = 0 for all the energy concentrated in one pixel.
    .. math::
        \text{E} = - \sum_{j=1}^N \frac{x_j}{x_\text{max}}
        \ln \left[\frac{x_j}{x_\text{max}}\right]
    with :math:`x_\text{max} = \sqrt{\sum_{j=1}^N x^2_j}`.
    The original equation is normalized by the maximum entropy, so that the
    :abbr:`EFC (Entropy Focus Criterion)` can be compared across images with
    different dimensions:
    .. math::
        \text{EFC} =
        \left( \frac{N}{\sqrt{N}} \, \log{\sqrt{N}^{-1}} \right) \text{E}
    :param numpy.ndarray img: input data
    :param numpy.ndarray framemask: a mask of empty voxels inserted after
                                    a rotation of data
    """

    if framemask is None:
        framemask = np.zeros_like(img, dtype=np.uint8)

    n_vox = np.sum(1 - framemask)
    # Calculate the maximum value of the EFC (which occurs any time all
    # voxels have the same value)
    efc_max = (
        1.0 * n_vox * (1.0 / np.sqrt(n_vox)) * np.log(1.0 / np.sqrt(n_vox))
    )

    # Calculate the total image energy
    b_max = np.sqrt((img[framemask == 0] ** 2).sum())

    # Calculate EFC (add 1e-16 to the image data to keep log happy)
    return float(
        (1.0 / efc_max)
        * np.sum(
            (img[framemask == 0] / b_max)
            * np.log((img[framemask == 0] + 1e-16) / b_max)
        )
    )


def fber(img, headmask, rotmask=None):
    r"""
    Calculate the
    :abbr:`FBER (Foreground-Background Energy Ratio)` [Shehzad2015]_,
    defined as the mean energy of image values within the head relative
    to outside the head. Higher values are better.
    .. math::
        \text{FBER} = \frac{E[|F|^2]}{E[|B|^2]}
    :param numpy.ndarray img: input data
    :param numpy.ndarray headmask: a mask of the head (including skull,
                                   skin, etc.)
    :param numpy.ndarray rotmask: a mask of empty voxels inserted after
                                  a rotation of data
    """

    fg_mu = np.median(np.abs(img[headmask > 0]) ** 2)

    airmask = np.ones_like(headmask, dtype=np.uint8)
    airmask[headmask > 0] = 0
    if rotmask is not None:
        airmask[rotmask > 0] = 0
    bg_mu = np.median(np.abs(img[airmask == 1]) ** 2)
    if bg_mu < 1.0e-3:
        return -1.0
    return float(fg_mu / bg_mu)


def find_peaks(data):
    """
    :param data: a numpy ndarray
    :returns: a list
    """
    t_z = [
        data[:, :, i, :].mean(axis=0).mean(axis=0)
        for i in range(data.shape[2])
    ]
    return t_z


def find_spikes(data, spike_thresh):
    """blabla"""
    data -= np.median(np.median(np.median(data, axis=0), axis=0), axis=0)
    slice_mean = np.median(np.median(data, axis=0), axis=0)
    t_z = _robust_zscore(slice_mean)
    spikes = np.abs(t_z) > spike_thresh
    spike_inds = np.transpose(spikes.nonzero())
    # mask out the spikes and recompute z-scores using variance uncontaminated
    # with spikes. This will catch smaller spikes that may have been swamped
    # by big ones.
    data.mask[:, :, spike_inds[:, 0], spike_inds[:, 1]] = True
    slice_mean2 = np.median(np.median(data, axis=0), axis=0)
    t_z = _robust_zscore(slice_mean2)

    spikes = np.logical_or(spikes, np.abs(t_z) > spike_thresh)
    spike_inds = [tuple(i) for i in np.transpose(spikes.nonzero())]
    return spike_inds, t_z


def fuzzy_jaccard(in_tpms, in_mni_tpms):
    """blabla"""
    overlaps = []
    for tpm, mni_tpm in zip(in_tpms, in_mni_tpms):
        tpm = tpm.reshape(-1)
        mni_tpm = mni_tpm.reshape(-1)

        num = np.min([tpm, mni_tpm], axis=0).sum()
        den = np.max([tpm, mni_tpm], axis=0).sum()
        overlaps.append(float(num / den))
    return overlaps


def gsr(epi_data, mask, direction="y", ref_file=None, out_file=None):
    """
    Compute the :abbr:`GSR (ghost to signal ratio)` [Giannelli2010]_.
    The procedure is as follows:
      #. Create a Nyquist ghost mask by circle-shifting the original
         mask by :math:`N/2`.
      #. Rotate by :math:`N/2`
      #. Remove the intersection with the original mask
      #. Generate a non-ghost background
      #. Calculate the :abbr:`GSR (ghost to signal ratio)`
    .. warning ::
      This should be used with EPI images for which the phase
      encoding direction is known.
    :param str epi_file: path to epi file
    :param str mask_file: path to brain mask
    :param str direction: the direction of phase encoding (x, y, all)
    :return: the computed gsr
    """
    direction = direction.lower()
    if direction[-1] not in ["x", "y", "all"]:
        raise Exception(
            "Unknown direction {}, should be one of x, -x, y, -y, all".format(
                direction
            )
        )

    if direction == "all":
        result = []
        for newdir in ["x", "y"]:
            ofile = None
            if out_file is not None:
                fname, ext = os.path.splitext(ofile)
                if ext == ".gz":
                    fname, ext2 = os.path.splitext(fname)
                    ext = ext2 + ext
                ofile = "{0}_{1}{2}".format(fname, newdir, ext)
            result += [
                gsr(epi_data, mask, newdir, ref_file=ref_file, out_file=ofile)
            ]
        return result

    # Roll data of mask through the appropriate axis
    axis = RAS_AXIS_ORDER[direction]
    n2_mask = np.roll(mask, mask.shape[axis] // 2, axis=axis)

    # Step 3: remove from n2_mask pixels inside the brain
    n2_mask = n2_mask * (1 - mask)

    # Step 4: non-ghost background region is labeled as 2
    n2_mask = n2_mask + 2 * (1 - n2_mask - mask)

    # Step 5: signal is the entire foreground image
    ghost = np.mean(epi_data[n2_mask == 1]) - np.mean(epi_data[n2_mask == 2])
    signal = np.median(epi_data[n2_mask == 0])
    return float(ghost / signal)


def image_binary_dilation(in_mask, radius=2):
    """
    Dilate the input binary mask.
    Parameters
    ----------
    in_mask: :obj:`numpy.ndarray`
        A 3D binary array.
    radius: :obj:`int`, optional
        The radius of the ball-shaped footprint for dilation of the mask.
    """
    from scipy import ndimage as ndi
    from skimage.morphology import ball

    return ndi.binary_dilation(in_mask.astype(bool), ball(radius)).astype(int)


def normalize_mc_params(params, source):
    """
    Normalize a single row of motion parameters to the SPM format.
    SPM saves motion parameters as:
        x   Right-Left          (mm)
        y   Anterior-Posterior  (mm)
        z   Superior-Inferior   (mm)
        rx  Pitch               (rad)
        ry  Roll                (rad)
        rz  Yaw                 (rad)
    """
    if source.upper() == "FSL":
        params = params[[3, 4, 5, 0, 1, 2]]
    elif source.upper() in ("AFNI", "FSFAST"):
        params = params[np.asarray([4, 5, 3, 1, 2, 0]) + (len(params) > 6)]
        params[3:] = params[3:] * np.pi / 180.0
    elif source.upper() == "NIPY":
        matrix = to_matrix44(params)
        params = np.zeros(6)
        params[:3] = matrix[:3, 3]
        params[-1:2:-1] = aff2euler(matrix)

    return params


def regress_poly(
    degree, data, remove_mean=True, axis=-1, failure_mode="error"
):
    """
    Returns data with degree polynomial regressed out.
    :param bool remove_mean: whether or not demean data (i.e. degree 0),
    :param int axis: numpy array axes along which regression is performed
    """

    datashape = data.shape
    timepoints = datashape[axis]
    if datashape[0] == 0 and failure_mode != "error":
        return data, np.array([])

    # Rearrange all voxel-wise time-series in rows
    data = data.reshape((-1, timepoints))

    # Generate design matrix
    X = np.ones((timepoints, 1))  # quick way to calc degree 0
    for i in range(degree):
        polynomial_func = Legendre.basis(i + 1)
        value_array = np.linspace(-1, 1, timepoints)
        X = np.hstack((X, polynomial_func(value_array)[:, np.newaxis]))

    non_constant_regressors = X[:, :-1] if X.shape[1] > 1 else np.array([])

    # Calculate coefficients
    betas = np.linalg.pinv(X).dot(data.T)

    # Estimation
    if remove_mean:
        datahat = X.dot(betas).T
    else:  # disregard the first layer of X, which is degree 0
        datahat = X[:, 1:].dot(betas[1:, ...]).T
    regressed_data = data - datahat

    # Back to original shape
    return regressed_data.reshape(datashape), non_constant_regressors


def rpve(pvms, seg):
    """
    * | Computes the :abbr:`rPVe (residual partial voluming error)`
      | of each tissue class.

    .. math ::
        \\text{rPVE}^k = \\frac{1}{N} \\left[ \\sum\\limits_{p^k_i \
        \\in [0.5, P_{98}]} p^k_i + \\sum\\limits_{p^k_i \
        \\in [P_{2}, 0.5)} 1 - p^k_i \\right]
    """

    pvfs = {}
    for k, lid in list(FSL_FAST_LABELS.items()):
        if lid == 0:
            continue
        pvmap = pvms[lid - 1]
        pvmap[pvmap < 0.0] = 0.0
        pvmap[pvmap >= 1.0] = 1.0
        totalvol = np.sum(pvmap > 0.0)
        upth = np.percentile(pvmap[pvmap > 0], 98)
        loth = np.percentile(pvmap[pvmap > 0], 2)
        pvmap[pvmap < loth] = 0
        pvmap[pvmap > upth] = 0
        pvfs[k] = (
            pvmap[pvmap > 0.5].sum() + (1.0 - pvmap[pvmap <= 0.5]).sum()
        ) / totalvol
    return {k: float(v) for k, v in list(pvfs.items())}


def snr(mu_fg, sigma_fg, n):
    r"""
    Calculate the :abbr:`SNR (Signal-to-Noise Ratio)`.
    The estimation may be provided with only one foreground region in
    which the noise is computed as follows:
    .. math::
        \text{SNR} = \frac{\mu_F}{\sigma_F\sqrt{n/(n-1)}},
    where :math:`\mu_F` is the mean intensity of the foreground and
    :math:`\sigma_F` is the standard deviation of the same region.
    :param float mu_fg: mean of foreground.
    :param float sigma_fg: standard deviation of foreground.
    :param int n: number of voxels in foreground mask.
    :return: the computed SNR
    """
    return float(mu_fg / (sigma_fg * sqrt(n / (n - 1))))


def snr_dietrich(mu_fg, mad_air=0.0, sigma_air=1.0):
    r"""
    Calculate the :abbr:`SNR (Signal-to-Noise Ratio)`.

    This must be an air mask around the head, and it should not contain
    artifacts. The computation is done following the eq. A.12 of
    [Dietrich2007]_, which includes a correction factor in the estimation
    of the standard deviation of air and its Rayleigh distribution:

    .. math::
        \text{SNR} =
        \frac{\mu_F}{\sqrt{\frac{2}{4-\pi}}\,\sigma_\text{air}}.

    :param float mu_fg: mean of foreground.
    :param float sigma_air: standard deviation of the air surrounding
                            the head ("hat" mask).

    :return: the computed SNR for the foreground segmentation

    """
    if mad_air > 1.0:
        return float(DIETRICH_FACTOR * mu_fg / mad_air)

    return (
        float(DIETRICH_FACTOR * mu_fg / sigma_air)
        if sigma_air > 1e-3
        else -1.0
    )


def summary_stats(img, pvms, airmask=None, erode=True):
    r"""
    Estimates the mean, the standard deviation, the 95\%
    and the 5\% percentiles of each tissue distribution.
    .. warning ::
        Sometimes (with datasets that have been partially processed), the air
        mask will be empty. In those cases, the background stats will be zero
        for the mean, median, percentiles and kurtosis, the sum of voxels in
        the other remaining labels for ``n``, and finally the MAD and the
        :math:`\sigma` will be calculated as:
        .. math ::
            \sigma_\text{BG} = \sqrt{\sum \sigma_\text{i}^2}
    """
    from statsmodels.robust.scale import mad

    # Check type of input masks
    dims = np.squeeze(np.array(pvms)).ndim
    if dims == 4:
        # If pvms is from FSL FAST, create the bg mask
        stats_pvms = [np.zeros_like(img)] + pvms
    elif dims == 3:
        stats_pvms = [np.ones_like(pvms) - pvms, pvms]
    else:
        raise RuntimeError(
            "Incorrect image dimensions ({0:d})".format(np.array(pvms).ndim)
        )

    if airmask is not None:
        stats_pvms[0] = airmask

    labels = list(FSL_FAST_LABELS.items())
    if len(stats_pvms) == 2:
        labels = list(zip(["bg", "fg"], list(range(2))))

    output = {}
    for k, lid in labels:
        mask = np.zeros_like(img, dtype=np.uint8)
        mask[stats_pvms[lid] > 0.85] = 1

        if erode:
            struc = nd.generate_binary_structure(3, 2)
            mask = nd.binary_erosion(mask, structure=struc).astype(np.uint8)

        nvox = float(mask.sum())
        if nvox < 1e3:
            print(
                'calculating summary stats of label "%s" in a very small '
                "mask (%d voxels)",
                k,
                int(nvox),
            )
            if k == "bg":
                output["bg"] = {
                    "mean": 0.0,
                    "median": 0.0,
                    "p95": 0.0,
                    "p05": 0.0,
                    "k": 0.0,
                    "stdv": 0.0,
                    "mad": 0.0,
                    "n": nvox,
                }
                continue

        output[k] = {
            "mean": float(img[mask == 1].mean()),
            "stdv": float(img[mask == 1].std()),
            "median": float(np.median(img[mask == 1])),
            "mad": float(mad(img[mask == 1])),
            "p95": float(np.percentile(img[mask == 1], 95)),
            "p05": float(np.percentile(img[mask == 1], 5)),
            "k": float(kurtosis(img[mask == 1])),
            "n": nvox,
        }

    return output


def volume_fraction(pvms):
    r"""
    Computes the :abbr:`ICV (intracranial volume)` fractions corresponding to
    the (partial volume maps).
    .. math ::
        \text{ICV}^k =
        \frac{\sum_i p^k_i}{\sum\limits_{x \in X_\text{brain}} 1}
    :param list pvms: list of :code:`numpy.ndarray` of partial volume maps.
    """
    tissue_vfs = {}
    total = 0
    for k, lid in list(FSL_FAST_LABELS.items()):
        if lid == 0:
            continue
        tissue_vfs[k] = pvms[lid - 1].sum()
        total += tissue_vfs[k]

    for k in list(tissue_vfs.keys()):
        tissue_vfs[k] /= total
    return {k: float(v) for k, v in list(tissue_vfs.items())}


def wm2max(img, mu_wm):
    r"""
    Calculate the :abbr:`WM2MAX (white-matter-to-max ratio)`,
    defined as the maximum intensity found in the volume w.r.t. the
    mean value of the white matter tissue. Values close to 1.0 are
    better:
    .. math ::
        \text{WM2MAX} = \frac{\mu_\text{WM}}{P_{99.95}(X)}
    """
    return float(mu_wm / np.percentile(img.reshape(-1), 99.95))


def _AR_est_YW(x, order, rxx=None):
    """Retrieve AR coefficients while dropping the sig_sq return value"""
    return AR_est_YW(x, order, rxx=rxx)[0]


def _flatten_dict(indict):
    """Reduce a dictionary of dictionaries to a dictionary

    Does not exceed two levels in depth.

    :param indict: a dictionary
    :returns: the flattened dictionary
    """
    out_qc = {}
    for k, value in list(indict.items()):
        if not isinstance(value, dict):
            out_qc[k] = value
        else:
            for subk, subval in list(value.items()):
                if not isinstance(subval, dict):
                    out_qc["_".join([k, subk])] = subval
                else:
                    for ssubk, ssubval in list(subval.items()):
                        out_qc["_".join([k, subk, ssubk])] = ssubval
    return out_qc


def _prepare_mask(mask, label, erode=True):
    """blabla


    :param mask: blabla
    :param label: blabla
    :param erode: blabla
    :returns: blabla
    """
    fgmask = mask.copy()

    if np.issubdtype(fgmask.dtype, np.integer):
        if isinstance(label, (str, bytes)):
            label = FSL_FAST_LABELS[label]

        fgmask[fgmask != label] = 0
        fgmask[fgmask == label] = 1
    else:
        fgmask[fgmask > 0.95] = 1.0
        fgmask[fgmask < 1.0] = 0

    if erode:
        # Create a structural element to be used in an opening operation.
        struc = nd.generate_binary_structure(3, 2)
        # Perform an opening operation on the background data.
        fgmask = nd.binary_opening(fgmask, structure=struc).astype(np.uint8)

    return fgmask


def _robust_zscore(data):
    """Blabla

    :return: blabla
    """
    return (data - np.atleast_2d(np.median(data, axis=1)).T) / np.atleast_2d(
        data.std(axis=1)
    ).T
