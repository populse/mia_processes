# -*- coding: utf-8 -*-

"""The mrtrix (mrtrix3) preprocess library of the mia_processes package.

The purpose of this module is to customise the main mrtrix preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - BrainMask
        - ConstrainedSphericalDeconvolution
        - DWIBiasCorrect
        - DWIDenoise
        - DWIPreproc
        - FitTensor
        - Generate5tt
        - MRDeGibbs
        - MRTransform -to do
        - ResponseSDDhollander
        - TensorMetrics - to do
        - Tractography - to do

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

from nipype.interfaces.base import File

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA
from traits.api import (
    Bool,
    Either,
    Enum,
    Float,
    Int,
    List,
    String,
    Tuple,
    Undefined,
)

from mia_processes.utils import checkFileExt

EXT = {"NIFTI_GZ": "nii.gz", "NIFTI": "nii"}


class BrainMask(ProcessMIA):
    """
    *Generates a whole brain mask from a DWI image.
    (dwi2mask command)*

    Please, see the complete documentation for the `BrainMask brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/BrainMask.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(BrainMask, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "Input DWI image (a pathlike object"
            "string representing an existing file)"
        )
        # Output descriptions
        out_file_desc = (
            "Output brain mask (a pathlike object or "
            "string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.BrainMask")

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
        super(BrainMask, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_brainmask." + in_ext
                )

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(BrainMask, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.out_file = self.out_file

        return self.process.run(configuration_dict={})


class ConstrainedSphericalDeconvolution(ProcessMIA):
    """
    *Estimate fibre orientation distributions from diffusion data using
    spherical deconvolution (dwi2fod command)*

    Please, see the complete documentation for the
    `ConstrainedSphericalDeconvolution brick in the
    populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/ConstrainedSphericalDeconvolution.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ConstrainedSphericalDeconvolution, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "Input DWI image (a pathlike object"
            "string representing an existing file)"
        )
        # Optional inputs description
        algorithm_desc = "FOD alogorithm (csd or msmt_csd)"
        in_directions_desc = (
            "Specify the directions over which to apply the non-negativity "
            "constraint (by default, the built-in 300 direction set is used)."
            "(a pathlike object or string representing an existing file)"
        )
        in_mask_desc = (
            "Only perform computation within the specified binary brain mask "
            "image.(a pathlike object or string representing an existing file)"
        )
        max_sh_desc = (
            "Maximum harmonic degree of response function (a single "
            "value for single-shell response and a list for multi-shell "
            "response. (a list of items which are an integer) "
        )
        get_predicted_signal_desc = (
            "Get a file with the predicted signal from the FOD "
            "estimates (only for msmt_csd algorithm)"
            "(a boolean) "
        )
        shell_desc = (
            "Specify one or more b-values to use during processing"
            "(a list of items which are a float)"
        )
        # Output descriptions
        predicted_signal_file_desc = (
            "Out predicted signal "
            "(a pathlike object or string representing a file)"
        )
        csf_odf_desc = (
            "Output CSF ODF."
            "(a pathlike object or string representing a file)"
        )
        csf_txt_desc = (
            "CSF response text file."
            "(a pathlike object or string representing a file)"
        )
        gm_odf_desc = (
            "Output GM ODF."
            "(a pathlike object or string representing a file)"
        )
        gm_txt_desc = (
            "GM response text file."
            "(a pathlike object or string representing a file)"
        )
        wm_odf_desc = (
            "Output WM ODF."
            "(a pathlike object or string representing a file)"
        )
        wm_txt_desc = (
            "WM response text file."
            "(a pathlike object or string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optionnal inouts traits
        self.add_trait(
            "algorithm",
            Enum(
                "csd",
                "msmt_csd",
                output=False,
                optional=True,
                desc=algorithm_desc,
            ),
        )

        self.add_trait(
            "in_directions",
            File(output=False, optional=True, desc=in_directions_desc),
        )

        self.add_trait(
            "in_mask", File(output=False, optional=True, desc=in_mask_desc)
        )

        self.add_trait(
            "max_sh",
            List(
                Int(),
                output=False,
                optional=True,
                desc=max_sh_desc,
            ),
        )

        self.add_trait(
            "shell",
            List(
                Float(),
                output=False,
                optional=True,
                desc=shell_desc,
            ),
        )

        self.add_trait(
            "get_predicted_signal",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_predicted_signal_desc,
            ),
        )

        # Ouputs traits
        self.add_trait(
            "csf_odf", File(output=True, optional=True, desc=csf_odf_desc)
        )

        self.add_trait(
            "csf_txt", File(output=True, optional=True, desc=csf_txt_desc)
        )

        self.add_trait(
            "gm_odf", File(output=True, optional=True, desc=gm_odf_desc)
        )

        self.add_trait(
            "gm_txt", File(output=True, optional=True, desc=gm_txt_desc)
        )

        self.add_trait(
            "wm_odf", File(output=True, optional=False, desc=wm_odf_desc)
        )

        self.add_trait(
            "wm_txt", File(output=True, optional=True, desc=wm_txt_desc)
        )

        self.add_trait(
            "predicted_signal_file",
            File(output=True, optional=False, desc=predicted_signal_file_desc),
        )

        self.init_default_traits()

        self.init_process(
            "nipype.interfaces.mrtrix3.ConstrainedSphericalDeconvolution"
        )

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
        super(ConstrainedSphericalDeconvolution, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["csf_odf"] = os.path.join(
                    self.output_directory, fileName + "_csf_odf." + in_ext
                )
                self.outputs["csf_txt"] = os.path.join(
                    self.output_directory, fileName + "_csf_response.txt"
                )
                self.outputs["gm_odf"] = os.path.join(
                    self.output_directory, fileName + "_gm_odf." + in_ext
                )
                self.outputs["gm_txt"] = os.path.join(
                    self.output_directory, fileName + "_gm_response.txt"
                )
                self.outputs["wm_odf"] = os.path.join(
                    self.output_directory, fileName + "_wm_odf." + in_ext
                )
                self.outputs["wm_txt"] = os.path.join(
                    self.output_directory, fileName + "_wm_response.txt"
                )

                if self.get_predicted_signal:
                    self.outputs["predicted_signal_file"] = os.path.join(
                        self.output_directory,
                        fileName + "_odf_predicted_signal." + in_ext,
                    )

        if self.outputs:
            self.inheritance_dict[self.outputs["csf_odf"]] = self.in_file
            self.inheritance_dict[self.outputs["csf_txt"]] = self.in_file
            self.inheritance_dict[self.outputs["gm_odf"]] = self.in_file
            self.inheritance_dict[self.outputs["gm_txt"]] = self.in_file
            self.inheritance_dict[self.outputs["wm_odf"]] = self.in_file
            self.inheritance_dict[self.outputs["wm_txt"]] = self.in_file
            if self.get_predicted_signal:
                self.inheritance_dict[
                    self.outputs["predicted_signal_file"]
                ] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ConstrainedSphericalDeconvolution, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.algorithm = self.algorithm
        self.process.wm_odf = self.wm_odf
        self.process.wm_txt = self.wm_txt
        self.process.gm_odf = self.gm_odf
        self.process.gm_txt = self.gm_txt
        self.process.csf_odf = self.csf_odf
        self.process.csf_txt = self.csf_txt

        if self.get_predicted_signal and self.algorithm == "msmt_csd":
            self.process.predicted_signal = self.predicted_signal_file

        if self.in_directions:
            self.process.in_dirs = self.in_directions

        if self.in_mask:
            self.process.mask_file = self.in_mask

        if self.shell:
            self.process.shell = self.shell

        if self.max_sh:
            self.process.max_sh = self.max_sh

        return self.process.run(configuration_dict={})


class DWIBiasCorrect(ProcessMIA):
    """
    *Perform B1 field inhomogeneity correction for a DWI volume series.
    (dwibias correct command)*

    Please, see the complete documentation for the `DWIBiasCorrect brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/DWIBiasCorrect.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(DWIBiasCorrect, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "Input DWI image (a pathlike object"
            "string representing an existing file)"
        )
        # Optionnal inputs description
        use_fsl_desc = (
            "Use FSL FAST to estimate the inhomogeneity field. (a boolean)"
        )
        use_ants_desc = (
            "Use ANTS N4 to estimate the inhomogeneity field. (a boolean)"
        )
        bias_desc = "Get bias field (a boolean)"
        in_mask_desc = (
            "Input mask image for bias field estimation (a pathlike object or "
            "string representing a file) "
        )
        # Output descriptions
        out_file_desc = (
            "The output denoised DWI image (a pathlike object or "
            "string representing a file)"
        )
        bias_field_map_desc = (
            "The output bias field map (a pathlike object or "
            "string representing a file)"
        )
        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "use_fsl",
            Bool(
                False,
                output=False,
                optional=True,
                desc=use_fsl_desc,
            ),
        )

        self.add_trait(
            "use_ants",
            Bool(
                True,
                output=False,
                optional=True,
                desc=use_ants_desc,
            ),
        )

        self.add_trait(
            "bias",
            Bool(
                True,
                output=False,
                optional=True,
                desc=bias_desc,
            ),
        )

        self.add_trait(
            "in_mask",
            Either(
                Undefined,
                File(),
                default=Undefined,
                output=False,
                optional=True,
                desc=in_mask_desc,
            ),
        )

        self.add_trait(
            "bias_field_map",
            File(output=True, optional=True, desc=bias_field_map_desc),
        )

        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.DWIBiasCorrect")

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
        super(DWIBiasCorrect, self).list_outputs()

        if self.use_ants and self.use_fsl:
            print("use_ants and use_fsl are mutually exclusif")
            return self.make_initResult()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_unbias." + in_ext
                )
                if self.bias:
                    self.outputs["bias_field_map"] = os.path.join(
                        self.output_directory,
                        fileName + "_bias_field_map." + in_ext,
                    )

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

            if self.bias:
                self.inheritance_dict[
                    self.outputs["bias_field_map"]
                ] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(DWIBiasCorrect, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.out_file = self.out_file
        if self.use_fsl:
            self.process.use_fsl = self.use_fsl
        if self.use_ants:
            self.process.use_ants = self.use_ants
        if self.in_mask:
            self.process.mask = self.in_mask
        if self.bias:
            self.process.bias = self.bias_field_map

        return self.process.run(configuration_dict={})


class DWIDenoise(ProcessMIA):
    """
    *Denoise DWI data and estimate the noise level based on the optimal
    threshold for PCA. (dwi denoise command)*

    Please, see the complete documentation for the `DWIDenoise brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/DWIDenoise.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(DWIDenoise, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "Input DWI image (a pathlike object"
            "string representing an existing file)"
        )
        # Optionnal inputs description
        extend_desc = (
            "Set the window size of the denoising filter. "
            " (a tuple of the form: (an integer, an integer, an integer), "
            "default = 5,5,5) "
        )
        mask_desc = (
            "Only process voxels within the specified binary brain mask "
            "image (a pathlike object or string representing an existing file)"
        )
        noise_desc = "Get noise map (a boolean)"
        # Outputs traits description
        noise_map_desc = (
            "The output noise map (a pathlike object or "
            "string representing a file)"
        )
        out_file_desc = (
            "The output denoised DWI image (a pathlike object or "
            "string representing a file)"
        )
        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "extend",
            Tuple(
                Int,
                Int,
                Int,
                default=(5, 5, 5),
                output=False,
                optional=True,
                desc=extend_desc,
            ),
        )

        self.add_trait(
            "noise",
            Bool(
                True,
                output=False,
                optional=True,
                desc=noise_desc,
            ),
        )

        self.add_trait(
            "mask",
            Either(
                Undefined,
                File(),
                default=Undefined,
                output=False,
                optional=True,
                desc=mask_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "noise_map", File(output=True, optional=True, desc=noise_map_desc)
        )

        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.DWIDenoise")

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
        super(DWIDenoise, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()
                return

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_denoised." + in_ext
                )
                if self.noise:
                    self.outputs["noise_map"] = os.path.join(
                        self.output_directory,
                        fileName + "_noise_map." + in_ext,
                    )

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

            if self.noise:
                self.inheritance_dict[self.outputs["noise_map"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(DWIDenoise, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.out_file = self.out_file
        self.process.extend = self.extend
        if self.mask:
            self.process.mask = self.mask
        if self.noise:
            self.process.noise = self.noise_map

        return self.process.run(configuration_dict={})


class DWIPreproc(ProcessMIA):
    """
    *Perform diffusion image pre-processing using FSL’s eddy tool;
    including inhomogeneity distortion correction using FSL’s topup
    tool if possible (dwifslpreproc correct command)*

    Please, see the complete documentation for the `DWIPreproc brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/DWIPreproc.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(DWIPreproc, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix", "fsl"]

        # Mandatory inputs description
        in_file_desc = (
            "Input DWI image (a pathlike object"
            "string representing an existing file)"
        )
        rpe_otions_desc = (
            "Specify acquisition phase-encoding design. “none” for no "
            "reversed phase-encoding image, “all” for all DWIs have opposing "
            "phase-encoding acquisition, “pair” for using a pair of b0 "
            "volumes for inhomogeneity field estimation only, and “header” "
            "for phase-encoding information can be found in the image "
            "header(s). (none or pair or all or header) "
        )
        # Optionnal inputs description
        align_seepi_desc = (
            "Achieve alignment between the SE-EPI images used for "
            "inhomogeneity field estimation, and the DWIs (a boolean)"
        )
        eddy_mask_desc = (
            "Provide a processing mask to use for eddy, instead of having "
            "dwifslpreproc generate one internally using dwi2mask "
            "(a pathlike object or string representing an existing file)"
        )
        eddy_option_desc = (
            "Additional command-line options to the eddy command (a string)"
        )
        eddy_slspec_desc = (
            "Provide a file containing slice groupings for eddy’s "
            "slice-to-volume registration.(a pathlike object or "
            "string representing a file) "
        )
        eddyqc_all_desc = (
            "Copy ALL outputs generated by eddy (including images), and the "
            "output of eddy_qc (if installed), into an output directory "
            "(a boolean)"
        )
        eddyqc_all_desc = (
            "Copy the various text-based statistical outputs generated "
            "by eddy, and the output of eddy_qc (if installed), "
            "into an output directory. (a boolean)"
        )
        epi_corr_desc = (
            "Provide an additional image series consisting of spin-echo EPI "
            "images which is to be used exclusively by topup for estimating "
            "the  inhomogeneity field (i.e. it will not form part of the "
            "output image series) (a pathlike object or string representing "
            "an existing file)"
        )
        pe_dir_desc = (
            "Specify the phase encoding direction of the input series, "
            "can be a signed axis number (e.g. -0, 1, +2), an axis designator "
            "(e.g. RL, PA, IS), or NIfTI axis codes (e.g. i-, j, k) (a string)"
        )
        ro_time_desc = (
            "Total readout time of input series (in seconds) (a float)"
        )
        topup_option_desc = (
            "Additional command-line options to the topup command (a string)"
        )
        # Output descriptions
        out_file_desc = (
            "The output preprocessed DWI image (a pathlike object or "
            "string representing a file)"
        )
        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "rpe_otions",
            Enum(
                "none",
                "pair",
                "all",
                "header",
                output=False,
                optional=True,
                desc=rpe_otions_desc,
            ),
        )

        # Optionnal inputs traits
        self.add_trait(
            "epi_corr", File(output=False, optional=True, desc=epi_corr_desc)
        )

        self.add_trait(
            "pe_dir",
            Enum(
                "ap",
                "pa",
                "lr",
                "rl",
                "is",
                "si",
                output=False,
                optional=True,
                desc=pe_dir_desc,
            ),
        )

        self.add_trait(
            "ro_time", Float(output=False, optional=True, desc=ro_time_desc)
        )

        self.add_trait(
            "align_seepi",
            Bool(
                False,
                output=False,
                optional=True,
                desc=align_seepi_desc,
            ),
        )

        self.add_trait(
            "eddy_mask",
            Either(
                Undefined,
                File(),
                default=Undefined,
                output=False,
                optional=True,
                desc=eddy_mask_desc,
            ),
        )

        self.add_trait(
            "eddy_option",
            String(
                "",
                output=False,
                optional=True,
                desc=eddy_option_desc,
            ),
        )

        self.add_trait(
            "eddy_slspec",
            Either(
                Undefined,
                File(),
                default=Undefined,
                output=False,
                optional=True,
                desc=eddy_slspec_desc,
            ),
        )

        self.add_trait(
            "eddyqc_all",
            Bool(
                False,
                output=False,
                optional=True,
                desc=eddyqc_all_desc,
            ),
        )

        self.add_trait(
            "eddyqc_text",
            Bool(
                False,
                output=False,
                optional=True,
                desc=eddyqc_all_desc,
            ),
        )

        self.add_trait(
            "topup_option",
            String(
                "",
                output=False,
                optional=True,
                desc=topup_option_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.DWIPreproc")

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
        super(DWIPreproc, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory,
                    fileName + "_dwifslpreproc." + in_ext,
                )

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(DWIBiasCorrect, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.out_file = self.out_file
        self.process.rpe_options = self.rpe_options
        self.process.align_seepi = self.align_seepi
        if self.epi_corr:
            self.process.in_epi = self.epi_corr
        self.process.pe_dir = self.pe_dir
        if self.ro_time:
            self.process.ro_time = self.ro_time
        if self.eddy_mask:
            self.process.eddy_mask = self.eddy_mask
        if self.eddy_options:
            self.process.eddy_options = self.eddy_options
        if self.eddy_slspec:
            self.process.eddy_slspec = self.eddy_slspec
        if self.eddyqc_all:
            self.process.eddyqc_all = os.path.dirname(self.out_file)
        if self.eddyqc_text:
            self.process.eddyqc_text = os.path.dirname(self.out_file)
        if self.topup_options:
            self.process.topup_options = self.topup_options

        return self.process.run(configuration_dict={})


class FitTensor(ProcessMIA):
    """
    *Convert diffusion-weighted images to tensor images.
    (dwi2tensor command)*

    Please, see the complete documentation for the `FitTensor brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/FitTensor.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(FitTensor, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "Input DWI image (a pathlike object"
            "string representing an existing file)"
        )
        # Optionnal inputs traits
        in_mask_desc = (
            "Only perform computation within the specified binary brain "
            "mask image (a pathlike object or string representing an "
            "existing file)"
        )
        get_predicted_signal_desc = (
            "Get a file with the predicted signal from the FOD "
            "estimates (only for msmt_csd algorithm)"
            "(a boolean) "
        )
        # Output descriptions
        out_file_desc = (
            "Output tensor image (a pathlike object or "
            "string representing a file)"
        )
        predicted_signal_file_desc = (
            "Out predicted signal "
            "(a pathlike object or string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optionnal inputs traits
        self.add_trait(
            "in_mask", File(output=False, optional=True, desc=in_mask_desc)
        )
        self.add_trait(
            "get_predicted_signal",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_predicted_signal_desc,
            ),
        )

        # TODO: add option dkt for DKI

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )
        self.add_trait(
            "predicted_signal_file",
            File(output=True, optional=True, desc=predicted_signal_file_desc),
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.FitTensor")

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
        super(FitTensor, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_dti." + in_ext
                )
                if self.get_predicted_signal:
                    self.outputs["predicted_signal_file"] = os.path.join(
                        self.output_directory,
                        fileName + "_dti_predicted_signal." + in_ext,
                    )

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file
            if self.get_predicted_signal:
                self.inheritance_dict[
                    self.outputs["predicted_signal_file"]
                ] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(FitTensor, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.out_file = self.out_file
        if self.get_predicted_signal:
            self.process.predicted_signal = self.predicted_signal_file

        return self.process.run(configuration_dict={})


class Generate5tt(ProcessMIA):
    """
    *Generate a 5TT image suitable for ACT using the selected algorithm.
    (5ttgen command)*

    Please, see the complete documentation for the `Generate5tt brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/Generate5tt.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Generate5tt, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix", "fsl", "freesurfer"]

        # Mandatory inputs description
        in_file_desc = (
            "Input image (a pathlike object"
            "string representing an existing file)"
        )
        algorithm_desc = "Tissue segmentation algorithm (fsl or freesurfer)"
        # Output descriptions
        out_file_desc = (
            "Output image (a pathlike object or " "string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs
        self.add_trait(
            "algorithm",
            Enum(
                "fsl",
                "freesurfer",
                output=False,
                optional=True,
                desc=algorithm_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.BrainMask")

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
        super(Generate5tt, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_5tt." + in_ext
                )

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Generate5tt, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.algorithm = self.algorithm
        self.process.out_file = self.out_file

        return self.process.run(configuration_dict={})


class MRDeGibbs(ProcessMIA):
    """
    *Remove Gibbs ringing artifacts. (mrdegibbs command)*

    Please, see the complete documentation for the `MRDeGibbs brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/DWIDenoise.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(MRDeGibbs, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "Input DWI image (a pathlike object"
            "string representing an existing file)"
        )
        # Optionnal inputs description
        axes_desc = (
            "Indicate the plane in which the data was acquire. "
            "(a list of items which are an integer, "
            "axial = [0, 1] ; coronal = [0, 2]; sagittal = [1, 2])"
        )
        maxW_desc = (
            "Right border of window used for total variation (TV) "
            "computation (an integer, default is 3)"
        )
        minW_desc = (
            "Left border of window used for total variation (TV) "
            "computation (an integer, default is 1)"
        )
        nshifts_desc = (
            "Discretization of subpixel spacing (an integer, default is 20)"
        )
        # Outputs description
        out_file_desc = (
            "The output unringed DWI image (a pathlike object or "
            "string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "axes",
            List(
                Int(),
                default=[0, 1],
                output=False,
                optional=True,
                desc=axes_desc,
            ),
        )

        self.add_trait(
            "maxW",
            Int(
                3,
                output=False,
                optional=True,
                desc=maxW_desc,
            ),
        )

        self.add_trait(
            "minW",
            Int(
                1,
                output=False,
                optional=True,
                desc=minW_desc,
            ),
        )

        self.add_trait(
            "nshifts",
            Int(
                20,
                output=False,
                optional=True,
                desc=nshifts_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.DWIDenoise")

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
        super(MRDeGibbs, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()
                return

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_unringed." + in_ext
                )

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(MRDeGibbs, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.out_file = self.out_file
        self.process.axes = self.axes
        self.process.maxW = self.maxW
        self.process.minW = self.minW
        self.process.nshifts = self.nshifts

        return self.process.run(configuration_dict={})


class ResponseSDDhollander(ProcessMIA):
    """
    *Estimate response function(s) for spherical deconvolution using
    the Dhollander algorithm (dwi2response command)*

    Please, see the complete documentation for the `ResponseSDDhollander brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/ResponseSDDhollander.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ResponseSDDhollander, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "Input DWI image (a pathlike object"
            "string representing an existing file)"
        )
        # Optionnal inputs description
        erode_desc = (
            "Number of erosion passes to apply to initial (whole brain) mask. "
            "Set to 0 to not erode the brain mask. (an intefer, default is 3)"
        )
        fa_desc = (
            "FA threshold for crude WM versus GM-CSF separation"
            "(a float, default is 0.2)"
        )
        get_final_voxels_desc = (
            "Get an image showing the final voxel selection(s) (a boolean)"
        )
        in_mask_desc = (
            "Provide initial mask image.(a pathlike object or "
            "string representing an existing file)"
        )
        max_sh_desc = (
            "Maximum harmonic degree of response function (a single "
            "value for single-shell response and a list for multi-shell "
            "response. (a list of items which are an integer) "
        )
        wm_algo_desc = (
            "Use external dwi2response algorithm for WM single-fibre voxel "
            "selection (fa, tax, tournier)"
        )

        # Outputs description
        csf_file_desc = (
            "Output CSF response text file (a pathlike object or "
            "string representing a file) "
        )
        gm_file_desc = (
            "Output GM response text file (a pathlike object or "
            "string representing a file) "
        )

        wm_file_desc = (
            "Output WM response text file (a pathlike object or "
            "string representing a file) "
        )
        voxels_image_desc = (
            "Image showing the final voxel selection (a pathlike object or "
            "string representing a file) "
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "erode",
            Int(
                3,
                output=False,
                optional=True,
                desc=erode_desc,
            ),
        )

        self.add_trait(
            "fa",
            Float(
                0.2,
                output=False,
                optional=True,
                desc=fa_desc,
            ),
        )

        self.add_trait(
            "get_final_voxels",
            Bool(
                True,
                output=False,
                optional=True,
                desc=get_final_voxels_desc,
            ),
        )

        self.add_trait(
            "in_mask",
            Either(
                Undefined,
                File(),
                default=Undefined,
                output=False,
                optional=True,
                desc=in_mask_desc,
            ),
        )

        self.add_trait(
            "max_sh",
            List(
                Int(),
                output=False,
                optional=True,
                desc=max_sh_desc,
            ),
        )

        self.add_trait(
            "wm_algo",
            Either(
                Undefined,
                Enum("fa", "tax", "tournier"),
                default=Undefined,
                output=False,
                optional=True,
                desc=wm_algo_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "csf_file", File(output=True, optional=False, desc=csf_file_desc)
        )
        self.add_trait(
            "gm_file", File(output=True, optional=False, desc=gm_file_desc)
        )
        self.add_trait(
            "wm_file", File(output=True, optional=False, desc=wm_file_desc)
        )
        self.add_trait(
            "voxels_image",
            File(output=True, optional=False, desc=voxels_image_desc),
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.DWIDenoise")

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
        super(ResponseSDDhollander, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()
                return

            if self.output_directory:
                self.outputs["csf_file"] = os.path.join(
                    self.output_directory, fileName + "_out_csf.txt"
                )
                self.outputs["gm_file"] = os.path.join(
                    self.output_directory, fileName + "_out_gm.txt"
                )
                self.outputs["wm_file"] = os.path.join(
                    self.output_directory, fileName + "_out_wm.txt"
                )
                if self.get_final_voxels:
                    self.outputs["voxels_image"] = os.path.join(
                        self.output_directory,
                        fileName + "response_voxels." + in_ext,
                    )

        if self.outputs:
            self.inheritance_dict[self.outputs["csf_file"]] = self.in_file
            self.inheritance_dict[self.outputs["gm_file"]] = self.in_file
            self.inheritance_dict[self.outputs["wm_file"]] = self.in_file
            if self.get_final_voxels:
                self.inheritance_dict[
                    self.outputs["voxels_image"]
                ] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ResponseSDDhollander, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.algorithm = "dhollander"
        self.process.csf_file = self.csf_file
        self.process.gm_file = self.gm_file
        self.process.wm_file = self.wm_file
        if self.in_mask:
            self.process.in_mask = self.in_mask
        if self.max_sh:
            self.process.max_sh = self.max_sh

        args = ""
        if self.erode != 3:
            args += "-erode " + self.erode + " "
        if self.fa != 0.2:
            args += "-fa " + self.fa + " "
        if self.wm_algo:
            args += "-wm_algo " + self.wm_algo + " "
        if self.get_final_voxels:
            args += "-voxels " + self.voxels_image + " "
        if args:
            self.process.args = args

        return self.process.run(configuration_dict={})
