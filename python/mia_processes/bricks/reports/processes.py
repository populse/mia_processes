"""The report preprocess library of the mia_processes package.

The purpose of this module is to provide bricks and functions to
compute necessary values for reporting.

:Contains:
    :Class:
        - AnatIQMs

    :Function:
        - art_qi1
        - art_qi2
        - cjv
        - cnr
        - efc
        - fber
        - fuzzy_jaccard
        - rpve
        - snr
        - snr_dietrich
        - summary_stats
        - volume_fraction
        - wm2max
"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# nibabel import
import nibabel as nb

# nipype import
from nipype.interfaces.base import (OutputMultiPath, InputMultiPath, File,
                                    traits, TraitListObject, Undefined,
                                    DictStrStr, Str)

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# soma-base imports
from soma.qt_gui.qt_backend.Qt import QMessageBox

# Other import
import os
import numpy as np
from math import pi, sqrt
import scipy.ndimage as nd
from scipy.stats import kurtosis  # pylint: disable=E0611
import json
import re

DIETRICH_FACTOR = 1.0 / sqrt(2 / (4 - pi))
FSL_FAST_LABELS = {"csf": 1, "gm": 2, "wm": 3, "bg": 0}


class AnatIQMs(ProcessMIA):
    """
        * Computes the anatomical IQMs *

        Please, see the complete documentation for the `AnatIQMs' brick in the populse.mia_processes website
        https://populse.github.io/mia_processes/documentation/bricks/preprocess/other/AnatIQMs.html

        """

    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(AnatIQMs, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_ras_desc = ('RAS input image (a pathlike object or string '
                       'representing a file).')
        airmask_desc = ('Air mask image (a pathlike object or string '
                        'representing a file).')
        artmask_desc = ('Artifact mask image (a pathlike object or string '
                        'representing a file).')
        headmask_desc = ('Head mask image (a pathlike object or string '
                        'representing a file).')
        rotmask_desc = ('Rotation mask image (a pathlike object or string '
                         'representing a file).')
        hatmask_desc = ('Hat mask image (a pathlike object or string '
                        'representing a file).')
        segmentation_desc = ('Segmentation mask image (a pathlike object or '
                             'string representing a file).')
        in_inu_desc = ('Input INU image (a pathlike object or string '
                       'representing a file).')
        in_noinu_desc = ('Input no-INU image (a pathlike object or string '
                         'representing a file).')
        pvms_desc = ('PVMS image (a pathlike object or string '
                     'representing a file).')
        mni_tpms_desc = ('MNI TPMS file (a pathlike '
                         'object or string representing a file).')
        in_fwhm_desc = ('FWHM (a float).')

        # Outputs description
        out_file_desc = 'a json file containing IQMs'

        # Inputs traits
        self.add_trait("in_ras",
                       File(output=False,
                            optional=False,
                            desc=in_ras_desc))
        self.add_trait("airmask",
                       File(output=False,
                            optional=True,
                            desc=airmask_desc))
        self.add_trait("artmask",
                       File(output=False,
                            optional=True,
                            desc=artmask_desc))
        self.add_trait("headmask",
                       File(output=False,
                            optional=True,
                            desc=headmask_desc))
        self.add_trait("rotmask",
                       File(output=False,
                            optional=True,
                            desc=rotmask_desc))
        self.add_trait("hatmask",
                       File(output=False,
                            optional=True,
                            desc=hatmask_desc))
        self.add_trait("segmentation",
                       File(output=False,
                            optional=True,
                            desc=segmentation_desc))
        self.add_trait("in_inu",
                       File(output=False,
                            optional=True,
                            desc=in_inu_desc))
        self.add_trait("in_noinu",
                       File(output=False,
                            optional=True,
                            desc=in_noinu_desc))
        self.add_trait("pvms",
                       traits.List(File(),
                                   output=False,
                                   optional=True,
                                   desc=pvms_desc))
        self.add_trait("mni_tpms",
                       traits.List(File(),
                                   output=False,
                                   optional=True,
                                   desc=mni_tpms_desc))
        self.add_trait("in_fwhm",
                       File(output=False,
                            optional=True,
                            desc=in_fwhm_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

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

            file_name = self.in_ras

            path, file_name = os.path.split(file_name)

            (file_name_no_ext,
             file_extension) = os.path.splitext(file_name)
            if file_extension == '.gz':
                (file_name_no_ext_2,
                 file_extension_2) = os.path.splitext(file_name_no_ext)
                if file_extension_2 == '.nii':
                    file_name_no_ext = file_name_no_ext_2

            report_file = os.path.join(self.output_directory,
                                       (file_name_no_ext +
                                        '_anat_qc.json'))

            if file_name:
                self.outputs['out_file'] = report_file

            else:
                print('- There was no output file deducted during '
                      'initialisation. Please check the input parameters...!')

            # tags inheritance (optional)
            if self.outputs:
                self.inheritance_dict[self.outputs[
                    'out_file']] = self.in_ras

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(AnatIQMs, self).run_process_mia()

        results_dict = {}

        imdata = nb.load(self.in_ras).get_data()

        # Try to load files
        has_in_noinu = True
        if self.in_noinu:
            try:
                imnii = nb.load(self.in_noinu)
            except (nb.filebasedimages.ImageFileError,
                    FileNotFoundError, TypeError) as e:
                has_in_noinu = False
                print("\nError with in_noinu file: ", e)
            else:
                # Load image corrected for INU
                inudata = np.nan_to_num(imnii.get_data())
                inudata[inudata < 0] = 0
        else:
            has_in_noinu = False

        has_segmentation = True
        if self.segmentation:
            try:
                # Load binary segmentation from FSL FAST
                segnii = nb.load(self.segmentation)
            except (nb.filebasedimages.ImageFileError,
                    FileNotFoundError, TypeError) as e:
                has_segmentation = False
                print("\nError with segmentation file: ", e)
                pass
            else:
                segdata = segnii.get_data().astype(np.uint8)
        else:
            has_segmentation = False

        has_airmask = True
        if self.airmask:
            try:
                # Load binary segmentation from FSL FAST
                airnii = nb.load(self.airmask)
            except (nb.filebasedimages.ImageFileError,
                    FileNotFoundError, TypeError) as e:
                has_airmask = False
                print("\nError with airmask file: ", e)
                pass
            else:
                airdata = airnii.get_data().astype(np.uint8)
        else:
            has_airmask = False

        has_artmask = True
        if self.artmask:
            try:
                # Load binary segmentation from FSL FAST
                artnii = nb.load(self.artmask)
            except (nb.filebasedimages.ImageFileError,
                    FileNotFoundError, TypeError) as e:
                has_artmask = False
                print("\nError with artmask file: ", e)
                pass
            else:
                artdata = artnii.get_data().astype(np.uint8)
        else:
            has_artmask = False

        has_headmask = True
        if self.headmask:
            try:
                # Load binary segmentation from FSL FAST
                headnii = nb.load(self.headmask)
            except (nb.filebasedimages.ImageFileError,
                    FileNotFoundError, TypeError) as e:
                has_headmask = False
                print("\nError with headmask file: ", e)
                pass
            else:
                headdata = headnii.get_data().astype(np.uint8)
        else:
            has_headmask = False

        has_rotmask = True
        if self.rotmask:
            try:
                # Load binary segmentation from FSL FAST
                rotnii = nb.load(self.rotmask)
            except (nb.filebasedimages.ImageFileError,
                    FileNotFoundError, TypeError) as e:
                has_rotmask = False
                print("\nError with rotmask file: ", e)
                pass
            else:
                rotdata = rotnii.get_data().astype(np.uint8)
        else:
            has_rotmask = False

        has_pvms = True
        if self.rotmask:
            try:
                # Load Partial Volume Maps (pvms) from FSL FAST
                pvmniis = []
                for fname in self.pvms:
                    pvmniis.append(nb.load(fname))
            except (nb.filebasedimages.ImageFileError,
                    FileNotFoundError, TypeError) as e:
                has_pvms = False
                print("\nError with pvms files: ", e)
                pass
            else:
                # Load Partial Volume Maps (pvms) from FSL FAST
                pvmdata = []
                for pvmnii in pvmniis:
                    pvmdata.append(pvmnii.get_data().astype(np.float32))
        else:
            has_pvms = False

        has_mni_tpms = True
        if self.mni_tpms:
            try:
                mni_tpmsniis = []
                for fname in self.mni_tpms:
                    mni_tpmsniis.append(nb.load(fname))
            except (nb.filebasedimages.ImageFileError,
                    FileNotFoundError, TypeError) as e:
                has_mni_tpms = False
                print("\nError with mni_tpms files: ", e)
                pass
        else:
            has_mni_tpms = False

        has_stats = False
        if has_in_noinu and has_pvms and has_airmask:
            erode = np.all(np.array(imnii.header.get_zooms()[:3], dtype=np.float32) < 1.9)

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
                tlabel: snr_dietrich(stats[tlabel]["median"], stats["bg"]["mad"])
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

        if has_airmask:
            # Artifacts QI2
            results_dict["qi_2"] = art_qi2(imdata, airdata)

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
        try:
            f = open(self.in_fwhm)
            lines = f.readlines()
            str_fwhm = re.findall(r"[-+]?\d*\.*\d+", lines[1])
            fwhm = []
            for item in str_fwhm:
                fwhm.append(float(item))
        except (FileNotFoundError, TypeError):
            print("\nError with fwhm file: ", e)
            fwhm = [0, 0, 0]
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
                i: float(v) for i, v in zip(["x", "y", "z"], imnii.header.get_zooms()[:3])
            }

            try:
                results_dict["size"]["t"] = int(inudata.shape[3])
            except IndexError:
                pass

            try:
                results_dict["spacing"]["tr"] = float(imnii.header.get_zooms()[3])
            except IndexError:
                pass

        if has_segmentation:
            # Bias
            bias = nb.load(self.in_inu).get_data()[segdata > 0]
            results_dict["inu"] = {
                "range": float(
                    np.abs(np.percentile(bias, 95.0) - np.percentile(bias, 5.0))
                ),
                "med": float(np.median(bias)),
            }  # pylint: disable=E1101

        if has_mni_tpms and has_pvms:
            mni_tpms = [tpm.get_data() for tpm in mni_tpmsniis]
            in_tpms = [pvm.get_data() for pvm in pvmniis]
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
        if file_extension == '.gz':
            (file_name_no_ext_2,
             file_extension_2) = os.path.splitext(file_name_no_ext)
            if file_extension_2 == '.nii':
                file_name_no_ext = file_name_no_ext_2

        report_file = os.path.join(self.output_directory,
                                   (file_name_no_ext +
                                    '_anat_qc.json'))

        with open(report_file, 'w') as fp:
            json.dump(flat_results_dict, fp)


class FWHMx(ProcessMIA):
    """
    * Computes FWHMs for all sub-bricks in the input dataset, each one separately *

    Please, see the complete documentation for the `FWHMx' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/afni/FWHMx.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(FWHMx, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_desc = ('Input image (a pathlike object or string '
                        'representing a file).')
        mask_file_desc = ('Mask image (a pathlike object or string '
                        'representing a file).')
        combine_desc = 'Combine the final measurements along each axis (a bool).'
        detrend_desc = 'detrend to the specified order (a bool or an int).'
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the output image file '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The output file (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("mask_file",
                       File(Undefined,
                            output=False,
                            optional=True,
                            desc=mask_file_desc))

        self.add_trait("combine",
                       traits.Bool(True,
                                   optional=True,
                                   output=False,
                                   desc=combine_desc))

        self.add_trait("detrend",
                       traits.Either(traits.Bool,
                                     traits.Int,
                                     default=True,
                                     output=False,
                                     desc=detrend_desc))

        self.add_trait("out_prefix",
                       traits.String('fwhm_',
                                     output=False,
                                     optional=True,
                                     desc=out_prefix_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.afni.FWHMx')

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
        super(FWHMx, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:

            if self.out_prefix == Undefined:
                self.out_prefix = 'fwhm_'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "fwhm" ...')

            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)
                    if trail == 'gz':
                        (fileName_2,
                         trail_2) = os.path.splitext(fileName)
                        if trail_2 == 'nii':
                            trail = 'nii.gz'

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:
                    self.outputs['out_file'] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName + '.out')

            else:
                print('No output_directory was found...!\n')
                return

            self.inheritance_dict[self.outputs[
                'out_file']] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(FWHMx, self).run_process_mia()

        self.process.in_file = self.in_file
        self.process.mask = self.mask_file
        self.process.combine = self.combine
        self.process.detrend = self.detrend
        self.process.out_file = self.out_file

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


class GCOR(ProcessMIA):
    """
    * Computes the average correlation between every voxel and ever other voxel, over any give mask *

    Please, see the complete documentation for the `GCOR' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/afni/GCOR.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(GCOR, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_desc = ('Input image (a pathlike object or string '
                        'representing a file).')
        mask_file_desc = ('Mask image (a pathlike object or string '
                        'representing a file).')

        # Outputs description
        out_desc = 'Global correlation value (a float).'

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("mask_file",
                       File(Undefined,
                            output=False,
                            optional=True,
                            desc=mask_file_desc))

        # Outputs traits
        self.add_trait("out",
                       traits.Float(output=True,
                                    desc=out_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.afni.GCOR')

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
        super(GCOR, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:

            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)
                    if trail == 'gz':
                        (fileName_2,
                         trail_2) = os.path.splitext(fileName)
                        if trail_2 == 'nii':
                            trail = 'nii.gz'

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:
                    self.outputs['out'] = 0
            else:
                print('No output_directory was found...!\n')
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(GCOR, self).run_process_mia()

        self.process.in_file = self.in_file
        self.process.mask = self.mask_file

        return self.process.run(configuration_dict={})


class OutlierCount(ProcessMIA):
    """
    * Computes FWHMs for all sub-bricks in the input dataset, each one separately *

    Please, see the complete documentation for the `FWHMx' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/afni/FWHMx.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(OutlierCount, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_desc = ('Input image (a pathlike object or string '
                        'representing a file).')
        mask_file_desc = ('Mask image (a pathlike object or string '
                        'representing a file).')
        fraction_desc = 'Combine the final measurements along each axis (a bool).'
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the output image file '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The output text file (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("mask_file",
                       File(Undefined,
                            output=False,
                            optional=True,
                            desc=mask_file_desc))

        self.add_trait("fraction",
                       traits.Bool(True,
                                   optional=True,
                                   output=False,
                                   desc=fraction_desc))

        self.add_trait("out_prefix",
                       traits.String('outliers_',
                                     output=False,
                                     optional=True,
                                     desc=out_prefix_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.afni.OutlierCount')

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
        super(OutlierCount, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:

            if self.out_prefix == Undefined:
                self.out_prefix = 'outliers_'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "outliers" ...')

            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)
                    if trail == 'gz':
                        (fileName_2,
                         trail_2) = os.path.splitext(fileName)
                        if trail_2 == 'nii':
                            trail = 'nii.gz'

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:
                    self.outputs['out_file'] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName + '.out')

            else:
                print('No output_directory was found...!\n')
                return

            self.inheritance_dict[self.outputs[
                'out_file']] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(OutlierCount, self).run_process_mia()

        self.process.in_file = self.in_file
        if self.mask_file:
            self.process.autoclip = Undefined
            self.process.automask = Undefined
            self.process.mask = self.mask_file

        self.process.fraction = self.fraction
        self.process.out_file = self.out_file

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


class QualityIndex(ProcessMIA):
    """
    * Computes a quality index for each sub-brick in a 3D+time dataset. The output is a 1D time
    series with the index for each sub-brick *

    Please, see the complete documentation for the `QualityIndex' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/afni/QualityIndex.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(QualityIndex, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_desc = ('A bold file (a pathlike object or string '
                        'representing a file).')
        automask_desc = 'Clip off small voxels (a boolean)'
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the output image file '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The output file (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("automask",
                       traits.Bool(True,
                                   output=False,
                                   optional=False,
                                   desc=automask_desc))

        self.add_trait("out_prefix",
                       traits.String('QI_',
                                     output=False,
                                     optional=True,
                                     desc=out_prefix_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.afni.QualityIndex')

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
        super(QualityIndex, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:

            if self.out_prefix == Undefined:
                self.out_prefix = 'QI_'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "QI" ...')

            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)
                    if trail == 'gz':
                        (fileName_2,
                         trail_2) = os.path.splitext(fileName)
                        if trail_2 == 'nii':
                            trail = 'nii.gz'

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:
                    self.outputs['out_file'] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName + '.out')

            else:
                print('No output_directory was found...!\n')
                return

            self.inheritance_dict[self.outputs[
                'out_file']] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(QualityIndex, self).run_process_mia()

        self.process.in_file = self.in_file
        self.process.automask = self.automask
        self.process.out_file = self.out_file

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


def art_qi1(airmask, artmask):
    r"""
    Detect artifacts in the image using the method described in [Mortamet2009]_.
    Caculates :math:`\text{QI}_1`, as the proportion of voxels with intensity
    corrupted by artifacts normalized by the number of voxels in the background:
    .. math ::
        \text{QI}_1 = \frac{1}{N} \sum\limits_{x\in X_\text{art}} 1
    Lower values are better.
    :param numpy.ndarray airmask: input air mask, without artifacts
    :param numpy.ndarray artmask: input artifacts mask
    """

    # Count the number of voxels that remain after the opening operation.
    # These are artifacts.
    return float(artmask.sum() / (airmask.sum() + artmask.sum()))


def art_qi2(img, airmask, min_voxels=int(1e3), max_voxels=int(3e5)):
    r"""
    Calculates :math:`\text{QI}_2`, based on the goodness-of-fit of a centered
    :math:`\chi^2` distribution onto the intensity distribution of
    non-artifactual background (within the "hat" mask):
    .. math ::
        \chi^2_n = \frac{2}{(\sigma \sqrt{2})^{2n} \, (n - 1)!}x^{2n - 1}\, e^{-\frac{x}{2}}
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

    modelx = data if len(data) < max_voxels else np.random.choice(data, size=max_voxels)

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

    return gof


def cjv(mu_wm, mu_gm, sigma_wm, sigma_gm):
    r"""
    Calculate the :abbr:`CJV (coefficient of joint variation)`, a measure
    related to :abbr:`SNR (Signal-to-Noise Ratio)` and
    :abbr:`CNR (Contrast-to-Noise Ratio)` that is presented as a proxy for
    the :abbr:`INU (intensity non-uniformity)` artifact [Ganzetti2016]_.
    Lower is better.
    .. math::
        \text{CJV} = \frac{\sigma_\text{WM} + \sigma_\text{GM}}{|\mu_\text{WM} - \mu_\text{GM}|}.
    :param float mu_wm: mean of signal within white-matter mask.
    :param float mu_gm: mean of signal within gray-matter mask.
    :param float sigma_wm: standard deviation of signal within white-matter mask.
    :param float sigma_gm: standard deviation of signal within gray-matter mask.
    :return: the computed CJV
    """
    return float((sigma_wm + sigma_gm) / abs(mu_wm - mu_gm))


def cnr(mu_wm, mu_gm, sigma_air):
    r"""
    Calculate the :abbr:`CNR (Contrast-to-Noise Ratio)` [Magnota2006]_.
    Higher values are better.
    .. math::
        \text{CNR} = \frac{|\mu_\text{GM} - \mu_\text{WM} |}{\sqrt{\sigma_B^2 +
        \sigma_\text{WM}^2 + \sigma_\text{GM}^2}},
    where :math:`\sigma_B` is the standard deviation of the noise distribution within
    the air (background) mask.
    :param float mu_wm: mean of signal within white-matter mask.
    :param float mu_gm: mean of signal within gray-matter mask.
    :param float sigma_air: standard deviation of the air surrounding the head ("hat" mask).
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
        \text{EFC} = \left( \frac{N}{\sqrt{N}} \, \log{\sqrt{N}^{-1}} \right) \text{E}
    :param numpy.ndarray img: input data
    :param numpy.ndarray framemask: a mask of empty voxels inserted after a rotation of
      data
    """

    if framemask is None:
        framemask = np.zeros_like(img, dtype=np.uint8)

    n_vox = np.sum(1 - framemask)
    # Calculate the maximum value of the EFC (which occurs any time all
    # voxels have the same value)
    efc_max = 1.0 * n_vox * (1.0 / np.sqrt(n_vox)) * np.log(1.0 / np.sqrt(n_vox))

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
    Calculate the :abbr:`FBER (Foreground-Background Energy Ratio)` [Shehzad2015]_,
    defined as the mean energy of image values within the head relative
    to outside the head. Higher values are better.
    .. math::
        \text{FBER} = \frac{E[|F|^2]}{E[|B|^2]}
    :param numpy.ndarray img: input data
    :param numpy.ndarray headmask: a mask of the head (including skull, skin, etc.)
    :param numpy.ndarray rotmask: a mask of empty voxels inserted after a rotation of
      data
    """

    fg_mu = np.median(np.abs(img[headmask > 0]) ** 2)

    airmask = np.ones_like(headmask, dtype=np.uint8)
    airmask[headmask > 0] = 0
    if rotmask is not None:
        airmask[rotmask > 0] = 0
    bg_mu = np.median(np.abs(img[airmask == 1]) ** 2)
    if bg_mu < 1.0e-3:
        return 0
    return float(fg_mu / bg_mu)


def fuzzy_jaccard(in_tpms, in_mni_tpms):
    overlaps = []
    for tpm, mni_tpm in zip(in_tpms, in_mni_tpms):
        tpm = tpm.reshape(-1)
        mni_tpm = mni_tpm.reshape(-1)

        num = np.min([tpm, mni_tpm], axis=0).sum()
        den = np.max([tpm, mni_tpm], axis=0).sum()
        overlaps.append(float(num / den))
    return overlaps


def rpve(pvms, seg):
    """
    Computes the :abbr:`rPVe (residual partial voluming error)`
    of each tissue class.
    .. math ::
        \\text{rPVE}^k = \\frac{1}{N} \\left[ \\sum\\limits_{p^k_i \
\\in [0.5, P_{98}]} p^k_i + \\sum\\limits_{p^k_i \\in [P_{2}, 0.5)} 1 - p^k_i \\right]
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


def snr_dietrich(mu_fg, sigma_air):
    r"""
    Calculate the :abbr:`SNR (Signal-to-Noise Ratio)`.
    This must be an air mask around the head, and it should not contain artifacts.
    The computation is done following the eq. A.12 of [Dietrich2007]_, which
    includes a correction factor in the estimation of the standard deviation of
    air and its Rayleigh distribution:
    .. math::
        \text{SNR} = \frac{\mu_F}{\sqrt{\frac{2}{4-\pi}}\,\sigma_\text{air}}.
    :param float mu_fg: mean of foreground.
    :param float sigma_air: standard deviation of the air surrounding the head ("hat" mask).
    :return: the computed SNR for the foreground segmentation
    """
    if sigma_air < 1.0:
        from .. import config

        config.loggers.interface.warning(
            f"SNRd - background sigma is too small ({sigma_air})"
        )
        sigma_air += 1.0

    return float(DIETRICH_FACTOR * mu_fg / sigma_air)


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

    if "bg" not in output:
        output["bg"] = {
            "mean": 0.0,
            "median": 0.0,
            "p95": 0.0,
            "p05": 0.0,
            "k": 0.0,
            "stdv": sqrt(sum(val["stdv"] ** 2 for _, val in list(output.items()))),
            "mad": sqrt(sum(val["mad"] ** 2 for _, val in list(output.items()))),
            "n": sum(val["n"] for _, val in list(output.items())),
        }

    if "bg" in output and output["bg"]["mad"] == 0.0 and output["bg"]["stdv"] > 1.0:
        print(
            "estimated MAD in the background was too small (MAD=%f)",
            output["bg"]["mad"],
        )
        output["bg"]["mad"] = output["bg"]["stdv"] / DIETRICH_FACTOR
    return output


def volume_fraction(pvms):
    r"""
    Computes the :abbr:`ICV (intracranial volume)` fractions
    corresponding to the (partial volume maps).
    .. math ::
        \text{ICV}^k = \frac{\sum_i p^k_i}{\sum\limits_{x \in X_\text{brain}} 1}
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


def _flatten_dict(indict):
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
