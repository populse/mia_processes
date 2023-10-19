# -*- coding: utf-8 -*-

"""The dipy preprocess library of the mia_processes package.

The purpose of this module is to customise the main dipy preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - Denoise
        - ComputeDKI

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

import os

import dipy.reconst.dki as dki
import nibabel as nib

# Other import
from dipy.core.gradients import gradient_table
from dipy.io.gradients import read_bvals_bvecs
from dipy.io.image import load_nifti

# nipype imports
from nipype.interfaces.base import File

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA
from traits.api import Either, Enum, Float, Int, String, Undefined

from mia_processes.utils import checkFileExt

EXT = {"NIFTI_GZ": "nii.gz", "NIFTI": "nii"}


class Denoise(ProcessMIA):
    """
    *Non-local means for denoising 3D images*

    Please, see the complete documentation for the `Denoise brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/dipy/Denoise.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Denoise, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype"]

        # Mandatory inputs description
        in_file_desc = (
            "A file to denoise (a pathlike object or string "
            "representing a file)."
        )
        # Optional inputs with default value description
        block_radius_desc = "Block_radius (an integer). Default is 5."
        noise_model_desc = (
            "Noise distribution model (‘rician’ or"
            " ‘gaussian’). Default is rician"
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the smoothed image file(s) "
            "(a string)."
        )
        patch_radius_desc = "Patch radius (an integer). Default is 1."
        # Optional inputs description
        in_mask_desc = (
            "Brain mask (a pathlike object or string representing a file)."
        )
        noise_mask_desc = (
            "Mask in which the mean signal will"
            " be computed. (a pathlike object or"
            "string representing a file)."
        )
        signal_mask_desc = (
            "Mask in which the standard deviation of noise"
            "will be computed. (a pathlike object or"
            "string representing a file)."
        )
        snr_desc = "Set manually Signal to noise ratio (a float)"
        # Outputs description
        out_file_desc = (
            "The denoised file (a pathlike object or a "
            "string representing a file)."
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits
        self.add_trait(
            "block_radius",
            Int(5, output=False, optional=True, desc=block_radius_desc),
        )

        self.add_trait(
            "noise_model",
            Enum(
                "rician",
                "gaussian",
                output=False,
                optional=True,
                desc=noise_model_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String(
                "denoise_", output=False, optional=True, desc=out_prefix_desc
            ),
        )

        self.add_trait(
            "patch_radius",
            Int(1, output=False, optional=True, desc=patch_radius_desc),
        )

        # Optional inputs value traits
        self.add_trait(
            "in_mask", File(output=False, optional=True, desc=in_mask_desc)
        )

        self.add_trait(
            "noise_mask",
            File(output=False, optional=True, desc=noise_mask_desc),
        )

        self.add_trait(
            "signal_mask",
            File(output=False, optional=True, desc=signal_mask_desc),
        )

        self.add_trait(
            "snr",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=snr_desc,
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.dipy.Denoise")

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
        super(Denoise, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if not self.out_prefix:
                self.out_prefix = "denoise_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "denoise" ...'
                )

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is not recognized...!")
                    return
                else:
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName + "." + in_ext,
                    )
                    self.tags_inheritance(
                        in_file=self.in_file,
                        out_file=self.outputs["out_file"],
                    )

            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Denoise, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.block_radius = self.block_radius
        self.process.noise_model = self.noise_model
        self.process._out_file = self.out_file
        self.process.patch_radius = self.patch_radius
        if self.in_mask:
            self.process.in_mask = self.in_mask
        if self.noise_mask:
            self.process.noise_mask = self.noise_mask
        if self.signal_mask:
            self.process.signal_mask = self.signal_mask
        if self.snr:
            self.process.snr = self.snr

        return self.process.run(configuration_dict={})


class ComputeDKI(ProcessMIA):
    """
    *Compute DKI*

    Please, see the complete documentation for the `ComputeDKI brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/dipy/ComputeDKI.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ComputeDKI, self).__init__()

        # Third party softwares required for the execution of the brick
        # self.requirement = ["nipype"]

        # Mandatory inputs description
        in_dwi_desc = (
            "Diffusion image(a pathlike object or string "
            "representing a file)."
        )
        dwi_bvec_desc = (
            "Bvec file (a pathlike object or string representing a file)."
        )
        dwi_bval_desc = (
            "Bval file (a pathlike object or string representing a file)."
        )
        in_mask_desc = (
            "Brain mask (a pathlike object or string representing a file)."
        )

        # Outputs description
        out_FA_desc = (
            "FA file (a pathlike object or a string representing a file)."
        )
        out_MD_desc = (
            "MD file (a pathlike object or a string representing a file)."
        )
        out_AD_desc = (
            "AD file (a pathlike object or a string representing a file)."
        )
        out_RD_desc = (
            "RD file (a pathlike object or a string representing a file)."
        )
        out_MK_desc = (
            "MK file (a pathlike object or a string representing a file)."
        )
        out_AK_desc = (
            "AK file (a pathlike object or a string representing a file)."
        )
        out_RK_desc = (
            "AK file (a pathlike object or a string representing a file)."
        )
        out_kFA_desc = (
            "kFA file (a pathlike object or a string representing a file)."
        )
        out_mKT_desc = (
            "mKT file (a pathlike object or a string representing a file)."
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_dwi", File(output=False, optional=False, desc=in_dwi_desc)
        )
        self.add_trait(
            "dwi_bvec", File(output=False, optional=False, desc=dwi_bvec_desc)
        )
        self.add_trait(
            "dwi_bval", File(output=False, optional=False, desc=dwi_bval_desc)
        )
        self.add_trait(
            "in_mask", File(output=False, optional=True, desc=in_mask_desc)
        )

        # Outputs traits
        self.add_trait("out_FA", File(output=True, desc=out_FA_desc))
        self.add_trait("out_MD", File(output=True, desc=out_MD_desc))
        self.add_trait("out_AD", File(output=True, desc=out_AD_desc))
        self.add_trait("out_RD", File(output=True, desc=out_RD_desc))
        self.add_trait("out_MK", File(output=True, desc=out_MK_desc))
        self.add_trait("out_AK", File(output=True, desc=out_AK_desc))
        self.add_trait("out_RK", File(output=True, desc=out_RK_desc))
        self.add_trait("out_kFA", File(output=True, desc=out_kFA_desc))
        self.add_trait("out_mKT", File(output=True, desc=out_mKT_desc))

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
        super(ComputeDKI, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_dwi:
            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_dwi, EXT)

                if not valid_ext:
                    print("\nThe input image format is not recognized...!")
                    return
                else:
                    self.outputs["out_FA"] = os.path.join(
                        self.output_directory,
                        fileName + "_dki_FA." + in_ext,
                    )
                    self.outputs["out_MD"] = os.path.join(
                        self.output_directory,
                        fileName + "_dki_MD." + in_ext,
                    )
                    self.outputs["out_AD"] = os.path.join(
                        self.output_directory,
                        fileName + "_dki_AD." + in_ext,
                    )
                    self.outputs["out_RD"] = os.path.join(
                        self.output_directory,
                        fileName + "_dki_RD." + in_ext,
                    )
                    self.outputs["out_MK"] = os.path.join(
                        self.output_directory,
                        fileName + "_dki_MK." + in_ext,
                    )
                    self.outputs["out_AK"] = os.path.join(
                        self.output_directory,
                        fileName + "_dki_AK." + in_ext,
                    )
                    self.outputs["out_RK"] = os.path.join(
                        self.output_directory,
                        fileName + "_dki_RK." + in_ext,
                    )
                    self.outputs["out_kFA"] = os.path.join(
                        self.output_directory,
                        fileName + "_dki_kFA." + in_ext,
                    )
                    self.outputs["out_mKT"] = os.path.join(
                        self.output_directory,
                        fileName + "_dki_mKT." + in_ext,
                    )
                    for out in [
                        "out_FA",
                        "out_MD",
                        "out_AD",
                        "out_RD",
                        "out_MK",
                        "out_AK",
                        "out_RK",
                        "out_kFA",
                        "out_mKT",
                    ]:
                        self.tags_inheritance(
                            in_file=self.in_dwi,
                            out_file=self.outputs[out],
                        )
            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ComputeDKI, self).run_process_mia()

        data, affine = load_nifti(self.in_dwi)
        bvals, bvecs = read_bvals_bvecs(self.dwi_bval, self.dwi_bvec)
        gtab = gradient_table(bvals, bvecs)
        dkimodel = dki.DiffusionKurtosisModel(gtab)
        mask, affine_mask = load_nifti(self.in_mask)
        dkifit = dkimodel.fit(data, mask=mask)
        dki_metrics = {
            "FA": dkifit.fa,
            "MD": dkifit.md,
            "AD": dkifit.ad,
            "RD": dkifit.rd,
            "MK": dkifit.mk(0, 3),
            "AK": dkifit.ak(0, 3),
            "RK": dkifit.rk(0, 3),
            "kFA": dkifit.kfa,
            "mKT": dkifit.mkt(0, 3),
        }

        paths = {
            "FA": self.out_FA,
            "MD": self.out_MD,
            "AD": self.out_AD,
            "RD": self.out_RD,
            "MK": self.out_MK,
            "AK": self.out_AK,
            "RK": self.out_RK,
            "kFA": self.out_kFA,
            "mKT": self.out_mKT,
        }

        for metric in list(dki_metrics.keys()):
            img = nib.Nifti1Image(dki_metrics[metric], affine)
            path = paths[metric]
            nib.save(img, path)
