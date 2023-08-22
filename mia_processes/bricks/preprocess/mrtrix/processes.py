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
        - DWIExtract
        - DWIPreproc
        - FitTensor
        - Generate5tt
        - MRCat
        - MRConvert
        - MRDeGibbs
        - MRMath
        - MRTransform
        - MTnormalise
        - ResponseSDDhollander
        - TensorMetrics
        - Tractography - to do

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################
import os

from capsul.in_context import mrtrix
from nipype.interfaces.base import File, InputMultiPath
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

EXT = {"NIFTI_GZ": "nii.gz", "NIFTI": "nii", "MIF": "mif"}


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


class DWIExtract(ProcessMIA):
    """
    *Extract diffusion-weighted volumes, b=0 volumes,
    or certain shells from a DWI dataset*

    Please, see the complete documentation for the `DWIExtract brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/mrtrix/DWIExtract.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(DWIExtract, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "Input image (a pathlike object"
            "string representing an existing file)"
        )
        # Optional inputs description
        bzero_desc = "Extract b=0 volume (a boolean)"
        nobzero_desc = "Extract non b=0 volume (a boolean)"
        shell_desc = (
            "Specify one or more gradient shells "
            "(a list of items which are a float)"
        )
        singleshell_desc = (
            "Extract volumes with a specific shell" "(a boolean)"
        )

        # Output description
        out_file_desc = (
            "Output image (a pathlike object or string representing "
            "an existing file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs traits
        self.add_trait(
            "bzero",
            Bool(
                True,
                output=False,
                optional=True,
                desc=bzero_desc,
            ),
        )

        self.add_trait(
            "nobzero",
            Bool(
                False,
                output=False,
                optional=True,
                desc=nobzero_desc,
            ),
        )

        self.add_trait(
            "shell",
            Either(
                Undefined,
                List(Float()),
                default=Undefined,
                output=False,
                optional=True,
                desc=shell_desc,
            ),
        )

        self.add_trait(
            "singleshell",
            Bool(
                False,
                output=False,
                optional=True,
                desc=singleshell_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.DWIExtract")

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
        super(DWIExtract, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                if self.nobzero:
                    fileName += "_nobzero"
                if self.bzero:
                    fileName += "_bzero"

                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "." + in_ext
                )

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(DWIExtract, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.out_file = self.out_file
        if self.bzero:
            self.process.bzero = self.bzero
        if self.nobzero:
            self.process.nobzero = self.nobzero
        if self.shell:
            self.process.shell = self.shell
        if self.singleshell:
            self.process.singleshell = self.singleshell

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


class MRCat(ProcessMIA):
    """
    *Concatenate several images into one*

    Please, see the complete documentation for the `MRCat brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/mrtrix/MRCat.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(MRCat, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_files_desc = (
            "Input images to concatenate (a list of items which are "
            "a pathlike object or a string representing an "
            "existing file)"
        )
        # Optional inputs description
        axis_desc = (
            "Specify axis along which concatenation should be performed "
            "(an integer, default is 3)"
        )

        out_file_name_desc = "Output file name (a string)"

        # Output description
        out_file_desc = (
            "Output image (a pathlike object or string representing "
            "an existing file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_files",
            InputMultiPath(
                Either(File(), List(File())),
                output=False,
                desc=in_files_desc,
            ),
        )

        # Optional inputs traits
        self.add_trait(
            "axis",
            Int(
                3,
                output=False,
                optional=True,
                desc=axis_desc,
            ),
        )

        self.add_trait(
            "out_file_name",
            String(
                "concatenated",
                output=False,
                optional=True,
                desc=out_file_name_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.MRCat")

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
        super(MRCat, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_files:
            in_ext = ""
            for in_file in self.in_files:
                valid_ext, in_ext_file, fileName = checkFileExt(in_file, EXT)
                if not in_ext:
                    in_ext = in_ext_file
                else:
                    if in_ext != in_ext_file:
                        print(
                            "\nAll input files should have the same "
                            "extension...!"
                        )
                        return self.make_initResult()

                if not valid_ext:
                    print("\nThe input image format is not recognized...!")
                    return self.make_initResult()

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, self.out_file_name + "." + in_ext
                )

        if self.outputs:
            # FIXME : out_file inherits only from first file
            self.inheritance_dict[self.outputs["out_file"]] = self.in_files[0]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(MRCat, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.out_file = self.out_file
        if self.bzero:
            self.process.bzero = self.bzero
        if self.nobzero:
            self.process.nobzero = self.nobzero
        if self.shell:
            self.process.shell = self.shell
        if self.singleshell:
            self.process.singleshell = self.singleshell

        return self.process.run(configuration_dict={})


class MRConvert(ProcessMIA):
    """
    *Perform conversion between different file types and optionally
    extract a subset of the input image. (mrconvert command)*

    Please, see the complete documentation for the `MRConvert brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/mrtrix/MRConvert.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(MRConvert, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "Input image (a pathlike object"
            "string representing an existing file)"
        )
        # Optional inputs description
        axes_desc = (
            "Specify the axes from the input image that will be used "
            "to form the output image (a list of items which are an integer)"
        )
        coord_desc = (
            "Extract data at the specified coordinates "
            "(a list of items which are an integer)"
        )
        scaling_desc = (
            "Specify the data scaling parameters used to rescale the "
            "intensity values(a list of items which are a float)"
        )
        vox_desc = (
            "Change the voxel dimensions reported in the output image header"
            "(a list of items which are a float)"
        )
        out_file_format_desc = (
            "Format of the output image (NIFTI, NIFTI_GZ or MIF )"
        )
        suffix_desc = "Output file suffix (a string, not mandatory)"
        # FIXME: json import / export --> in mia it is not the same json as
        #        used in BIDS, see if it is usefull/ possible to add
        #        thoses option

        # Optional base inputs description
        bval_scale_desc = (
            "Specifies whether the b - values should be scaled by the square "
            "of the corresponding DW gradient norm, as often required for"
            "multishell or DSI DW acquisition schemes. "
            "(yes or no, default is yes)"
        )
        grad_file_desc = (
            "Provide the diffusion-weighted gradient scheme used in the "
            "acquisition in a text file (MRTrix format) (a pathlike object "
            "or string representing an existing file) "
        )
        # grad_fsl_desc = (
        #     "Provide the diffusion-weighted gradient scheme used in the "
        #     "acquisition in FSL bvecs/bvals format files. (a tuple of "
        #     "the form: (a pathlike object or string representing an "
        #     "existing file, a pathlike object or string representing "
        #     "an existing file), it should be (bvecs, bvals)) "
        # )
        in_bvec_desc = (
            "Bvecs file in FSL format (a pathlike object or string "
            "representing an existing file) "
        )
        in_bval_desc = (
            "Bvals file in FSL format (a pathlike object or string "
            "representing an existing file) "
        )
        export_bvec_bval_desc = (
            "Export bvec / bval files in FSL format (a boolean)"
        )

        # Output description
        out_file_desc = (
            "Output image (a pathlike object or string representing "
            "an existing file)"
        )
        out_bvec_desc = (
            "bvec file in FSL format (a pathlike object or "
            "string representing a file)"
        )
        out_bval_desc = (
            "bval file in FSL format (a pathlike object or "
            "string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs traits
        self.add_trait(
            "axes",
            Either(
                Undefined,
                List(Int()),
                default=Undefined,
                output=False,
                optional=True,
                desc=axes_desc,
            ),
        )

        self.add_trait(
            "coord",
            Either(
                Undefined,
                List(Int()),
                default=Undefined,
                output=False,
                optional=True,
                desc=coord_desc,
            ),
        )

        self.add_trait(
            "scaling",
            Either(
                Undefined,
                List(Float()),
                default=Undefined,
                output=False,
                optional=True,
                desc=scaling_desc,
            ),
        )

        self.add_trait(
            "vox",
            Either(
                Undefined,
                List(Float()),
                default=Undefined,
                output=False,
                optional=True,
                desc=vox_desc,
            ),
        )

        self.add_trait(
            "out_file_format",
            Enum(
                "MIF",
                "NIFTI",
                "NIFTI_GZ",
                output=False,
                optional=True,
                desc=out_file_format_desc,
            ),
        )

        self.add_trait(
            "export_bvec_bval",
            Bool(
                True, output=False, optional=True, desc=export_bvec_bval_desc
            ),
        )
        self.add_trait(
            "suffix",
            Either(
                Undefined,
                String(),
                default=Undefined,
                output=False,
                optional=True,
                desc=suffix_desc,
            ),
        )

        # Optional base inputs traits
        self.add_trait(
            "bval_scale",
            Enum(
                "no", "yes", output=False, optional=True, desc=bval_scale_desc
            ),
        )

        self.add_trait(
            "grad_file", File(output=False, optional=True, desc=grad_file_desc)
        )

        # FIXME: in_bvec and in_bval already in grad-fsl ? (command -fslgrad)
        # self.add_trait(
        #     "grad_fsl",
        #     Tuple(
        #         File(), File(), output=False, optional=True,
        #         desc=grad_fsl_desc
        #     ),
        # )

        self.add_trait(
            "in_bvec", File(output=False, optional=True, desc=in_bvec_desc)
        )

        self.add_trait(
            "in_bval", File(output=False, optional=True, desc=in_bval_desc)
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )
        self.add_trait(
            "out_bvec", File(output=True, optional=True, desc=out_bvec_desc)
        )

        self.add_trait(
            "out_bval", File(output=True, optional=True, desc=out_bval_desc)
        )

        self.init_default_traits()

        # Nipype command not used because only working
        # for diffusion images (with bvec / bval)

        # self.init_process("nipype.interfaces.mrtrix3.MRConvert")

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
        super(MRConvert, self).list_outputs()

        if (self.in_bval and not self.in_bvec) or (
            not self.in_bval and self.in_bvec
        ):
            print("\nIf grad_file used, do not provied bvec or bval")
            return self.make_initResult()

        if self.grad_file and self.in_bvec:
            print("\nIf grad_file used, do not provied bvec or bval")
            return self.make_initResult()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                if self.suffix:
                    fileName += "_" + self.sufix
                else:
                    if self.scaling:
                        fileName += "_scaled"
                    if self.vox:
                        fileName += "_vox"
                    if self.coord:
                        fileName += "_coord"
                    if self.axes:
                        fileName += "_axes"

                self.outputs["out_file"] = os.path.join(
                    self.output_directory,
                    fileName + "." + EXT[self.out_file_format],
                )

                if self.export_bvec_bval:
                    self.outputs["out_bvec"] = os.path.join(
                        self.output_directory, fileName + ".bvec"
                    )
                    self.outputs["out_bval"] = os.path.join(
                        self.output_directory, fileName + ".bval"
                    )

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file
            if self.export_bvec_bval:
                self.inheritance_dict[self.outputs["out_bvec"]] = self.in_file
                self.inheritance_dict[self.outputs["out_bval"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(MRConvert, self).run_process_mia()
        # Nipype command not used because only working
        # for diffusion images (with bvec / bval)

        # self.process.in_file = self.in_file
        # self.process.out_file = self.out_file
        # if self.axes:
        #     self.process.axes = self.axes
        # if self.coord:
        #     self.process.coord = self.coord
        # if self.scaling:
        #     self.process.scaling = self.scaling
        # if self.vox:
        #     self.process.vox = self.vox
        # if self.grad_file:
        #     self.process.grad_file = self.grad_file
        # if self.in_bvec:
        #     self.process.grad_fsl = (self.in_bvec, self.in_bval)
        # self.process.bval_scale = self.bval_scale
        # if self.export_bvec_bval:
        #     self.process.out_bvec = self.out_bvec
        #     self.process.out_bval = self.out_bval

        # return self.process.run(configuration_dict={})

        cmd = ["mrconvert", self.in_file]

        if self.axes:
            cmd += ["-axes", self.axes]
        if self.coord:
            cmd += ["-coord", self.coord]
        if self.scaling:
            cmd += ["-scaling", self.scaling]
        if self.vox:
            cmd += ["-vox", self.vox]
        if self.grad_file:
            cmd += ["-grad", self.axes]
        if self.in_bvec:
            cmd += ["-fslgrad", self.in_bvec, self.in_bval]
        if self.bval_scale == "yes":
            cmd += ["-bvalue_scaling", self.bval_scale]
        if self.export_bvec_bval:
            cmd += ["-export_grad_fsl", self.out_bvec, self.out_bval]

        cmd += [self.out_file]

        return mrtrix.mrtrix_call(cmd)


class MRDeGibbs(ProcessMIA):
    """
    *Remove Gibbs ringing artifacts. (mrdegibbs command)*

    Please, see the complete documentation for the `MRDeGibbs brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/MRDeGibbs.html>`_
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

        self.init_process("nipype.interfaces.mrtrix3.MRDeGibbs")

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


class MRMath(ProcessMIA):
    """
    *Compute summary statistic on image intensities along a
    specified axis of a single image*

    Please, see the complete documentation for the `MRMath brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/mrtrix/MRMath.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(MRMath, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "Input image (a pathlike object or a string representing an "
            "existing file)"
        )
        operation_desc = (
            "Operation to computer along a specified axis "
            "(mean or median or sum or product or rms or norm or var or std "
            "or min or max or absmax or magmax)"
        )
        # Optional inputs description
        axis_desc = (
            "Specify axis along which concatenation should be performed "
            "(an integer, default is 3)"
        )
        out_file_name_desc = "Output file name (a string)"
        # Output description
        out_file_desc = (
            "Output image (a pathlike object or string representing "
            "an existing file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "operation",
            Enum(
                "mean",
                "median",
                "sum",
                "product",
                "rms",
                "norm",
                "var",
                "sdt",
                "min",
                "max",
                "absmax",
                "absmin",
                "magmax",
                output=False,
                optional=False,
                desc=operation_desc,
            ),
        )

        # Optional inputs traits
        self.add_trait(
            "axis",
            Int(
                3,
                output=False,
                optional=True,
                desc=axis_desc,
            ),
        )

        self.add_trait(
            "out_file_name",
            Either(
                Undefined,
                String(),
                default=Undefined,
                output=False,
                optional=True,
                desc=out_file_name_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.MRMath")

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
        super(MRMath, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                if self.out_file_name:
                    name = self.out_file_name
                else:
                    name = fileName + "_" + self.operation

                self.outputs["out_file"] = os.path.join(
                    self.output_directory, name + "." + in_ext
                )

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(MRMath, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.operation = self.operation
        self.process.out_file = self.out_file
        self.process.axis = self.axis

        return self.process.run(configuration_dict={})


class MRTransform(ProcessMIA):
    """
    *Apply spatial transformations or reslice images (mrtransform command)*

    Please, see the complete documentation for the `MRTransform brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/MRTransform.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(MRTransform, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Nipype MRTransform process is not up to date for some option
        # so the "args" input is used for several options

        # Mandatory inputs description
        in_file_desc = (
            "Input images to be transformed (a "
            "pathlike object or string representing an existing file)"
        )
        # Optionnal inputs description
        linear_transform_desc = (
            "Specify a linear transform to apply "
            "(a pathlike object or string representing an existing file) "
        )
        inverse_desc = (
            "Invert the specified transform before using it " "(a boolean)"
        )
        flip_axes_desc = (
            "flip the specified axes " "(a list of int with 0:x, 1:y and 2:z)"
        )
        half_desc = (
            "Apply the matrix square root of the transformation" "(a boolean)"
        )
        replace_file_desc = (
            "Replace the linear transform of the original image by that "
            "specified, rather than applying it to the original image "
            "(a pathlike object or string representing an existing file)"
        )
        identity_desc = (
            "Set the header transform of the image to the identity matrix"
            "(a boolean)"
        )
        template_image_desc = (
            "Reslice the input image to match the specified template image. "
            "(a pathlike object or string representing an existing file)"
        )
        midway_space_desc = (
            "Reslice the input image to the midway space. Requires either the "
            "-template or -warp option. (a boolean)"
        )
        interpolation_desc = (
            "Set the interpolation method to use when reslicing "
            "(cubic, nearest, linear, sinc)"
        )
        oversample_factor_desc = (
            "Set the amount of over-sampling (in the target space) to perform "
            "when regridding (an integer or a list of 3 integers)"
        )
        warp_image_desc = (
            "Apply a non-linear 4D deformation field to warp the input image"
            "(a pathlike object or string representing an existing file)"
        )
        warp_full_image_desc = (
            "Warp the input image using a 5D warp file output from mrregister"
            "(a pathlike object or string representing an existing file)"
        )
        fod_modulate_desc = "Intensity modulation method for fod (fod or jac)"
        fod_directions_file_desc = (
            "Directions defining the number and orientation of the apodised"
            "point spread functions used in FOD reorientation"
            "(a pathlike object or string representing an existing file)"
        )
        fod_reorient_desc = (
            "Specify whether to perform FOD reorientation" "(a boolean)"
        )

        # Outputs description
        out_file_desc = (
            " The output image of the transformation."
            "(a pathlike object or string representing an existing file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs traits
        self.add_trait(
            "linear_transform",
            File(output=False, optional=True, desc=linear_transform_desc),
        )

        self.add_trait(
            "inverse",
            Bool(
                False,
                output=False,
                optional=True,
                desc=inverse_desc,
            ),
        )

        self.add_trait(
            "half",
            Bool(
                False,
                output=False,
                optional=True,
                desc=half_desc,
            ),
        )

        self.add_trait(
            "identity",
            Bool(
                False,
                output=False,
                optional=True,
                desc=identity_desc,
            ),
        )

        self.add_trait(
            "flip_axes",
            Either(
                Undefined,
                List(),
                default=Undefined,
                output=False,
                optional=True,
                desc=flip_axes_desc,
            ),
        )

        self.add_trait(
            "replace_file",
            File(output=False, optional=True, desc=replace_file_desc),
        )

        self.add_trait(
            "template_image",
            File(output=False, optional=True, desc=template_image_desc),
        )

        self.add_trait(
            "midway_space",
            Bool(
                False,
                output=False,
                optional=True,
                desc=midway_space_desc,
            ),
        )

        self.add_trait(
            "interpolation",
            Enum(
                "cubic",
                "nearest",
                "linear",
                "sinc",
                False,
                output=False,
                optional=True,
                desc=interpolation_desc,
            ),
        )

        self.add_trait(
            "oversample_factor",
            Either(
                Undefined,
                List(),
                List(Int()),
                default=Undefined,
                output=False,
                optional=True,
                desc=oversample_factor_desc,
            ),
        )

        self.add_trait(
            "warp_image",
            File(output=False, optional=True, desc=warp_image_desc),
        )

        self.add_trait(
            "warp_full_image",
            File(output=False, optional=True, desc=warp_full_image_desc),
        )

        self.add_trait(
            "fod_modulate",
            Either(
                Undefined,
                Enum("fod", "jac"),
                default=Undefined,
                output=False,
                optional=True,
                desc=fod_modulate_desc,
            ),
        )

        self.add_trait(
            "fod_directions_file",
            File(output=False, optional=True, desc=fod_directions_file_desc),
        )

        self.add_trait(
            "fod_reorient",
            Bool(
                False,
                output=False,
                optional=True,
                desc=fod_reorient_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=True, desc=out_file_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.MRTransform")

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
        super(MRTransform, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()
                return

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_transformed." + in_ext
                )

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(MRTransform, self).run_process_mia()
        self.process.in_files = self.in_file
        self.process.out_file = self.out_file
        if self.linear_transform:
            self.process.linear_transform = self.linear_transform
        if self.inverse:
            self.process.invert = self.inverse
        if self.template_image:
            self.process.template_image = self.template_image
        args = ""
        if self.flip_axes:
            args += "-flip " + self.flip_axes
        if self.half:
            args += "-half "
        if self.replace_file:
            args += "-replace " + self.replace_file + " "
        if self.identify:
            args += "-identity "
        if self.midway_space:
            args += "-midway_space "
        if self.interpolation:
            args += "-interp " + self.interpolation + " "
        if self.oversample_factor:
            args += "-oversample " + self.oversample_factor + " "
        if self.warp_image:
            args += "-warp " + self.warp_image + " "
        if self.warp_full_image:
            args += "-warp_full " + self.warp_full_image + " "
        if self.fod_modulate:
            args += "-modulate " + self.fod_modulate + " "
        if self.fod_directions_file:
            args += "-directions " + self.fod_directions_file + " "
        if self.fod_reoriente:
            args += "-reorient_fod "
        if args:
            self.process.args = args

        return self.process.run(configuration_dict={})


class MTnormalise(ProcessMIA):
    """
    *Multi-tissue informed log-domain intensity normalisation.
    (mtnormalise command)*

    Please, see the complete documentation for the `MTnormalise brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/MTnormalise.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(MTnormalise, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["mrtrix"]
        # process nit in nipype

        # Mandatory inputs description
        in_files_desc = (
            "Input tissue component (a list of items which are a pathlike "
            "object or a string representing an existing file)"
        )
        mask_desc = (
            "The mask defines the data used to compute the intensity "
            "normalisation (a pathlike object or "
            "string representing a file)"
        )
        # Optionnal inputs description
        order_number_desc = (
            "The maximum order of the polynomial basis used to fit the "
            "normalisation field in the log-domain. "
            "(an integer, default is 3)"
        )
        niter_number_desc = (
            "Number of iteration (an integer or a list of integer, "
            "default is [15, 7])"
        )
        reference_number_desc = (
            "The (positive) reference value to which the summed tissue "
            "compartments will be normalised (a float, default is 0.282095) "
        )
        balanced_desc = (
            "Incorporate the per-tissue balancing factors into scaling of "
            "the output images (a bollean)"
        )

        # Outputs description
        out_files_desc = (
            "Nomalised outputs images(a pathlike object or "
            "string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_files",
            InputMultiPath(
                Either(File(), List(File())),
                output=False,
                desc=in_files_desc,
            ),
        )

        self.add_trait("mask", File(output=False, desc=mask_desc))

        # Optional inputs traits
        self.add_trait(
            "order_number",
            Int(
                3,
                output=False,
                optional=True,
                desc=order_number_desc,
            ),
        )

        self.add_trait(
            "niter_number",
            Either(
                Int(),
                List(Int()),
                default=[15, 7],
                output=False,
                optional=True,
                desc=niter_number_desc,
            ),
        )

        self.add_trait(
            "reference_number",
            Float(
                0.282095,
                output=False,
                optional=True,
                desc=reference_number_desc,
            ),
        )

        self.add_trait(
            "balanced",
            Bool(
                False,
                output=False,
                optional=True,
                desc=balanced_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_files",
            InputMultiPath(
                Either(File(), List(File())),
                output=False,
                desc=out_files_desc,
            ),
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
        super(MTnormalise, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        out_files = []
        if self.in_files:
            for f in self.in_files:
                valid_ext, in_ext, fileName = checkFileExt(f, EXT)

                if not valid_ext:
                    print("\nThe input image format is not recognized...!")
                    return self.make_initResult()

                if self.output_directory:
                    out_file = os.path.join(
                        self.output_directory, fileName + "_norm." + in_ext
                    )
                    out_files.append(out_file)

        if out_files:
            self.outputs["out_files"] = out_files

        if self.outputs:
            # FIXME: out_files inherits only from the first file
            self.inheritance_dict[self.outputs["out_files"]] = self.in_files[0]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(MTnormalise, self).run_process_mia()

        cmd = ["mtnormalise"]
        i = 0
        for in_file in self.in_files:
            cmd += [in_file, self.out_files[i]]
            i += 1
        if self.order_number:
            cmd += ["-order", self.order_number]
        if self.niter_number:
            niter = ""
            for n in self.niter_number:
                niter += str(n) + ","
            cmd += ["-niter", niter[:-1]]
        if self.reference_number:
            cmd += ["-reference", self.reference_number]
        if self.balanced:
            cmd += ["-balanced"]

        return mrtrix.mrtrix_call(cmd)


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

        self.init_process("nipype.interfaces.mrtrix3.ResponseSD")

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


class TensorMetrics(ProcessMIA):
    """
    *Compute metrics from tensors (tensor2metric command)*

    Please, see the complete documentation for the `TensorMetrics brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/TensorMetrics.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(TensorMetrics, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_dti_desc = (
            "Input DTI image (a pathlike object"
            "string representing an existing file)"
        )
        # Optionnal inputs description
        component_desc = (
            "Specify the desired eigenvalue/eigenvector(s) "
            "(a list of items which are any value)"
        )

        in_mask_desc = (
            "Only perform computation within the specified binary brain mask "
            "image.(a pathlike object or "
            "string representing an existing file)"
        )
        modulate_desc = (
            "Specify how to modulate the magnitude of the eigenvectors. "
            "(FA, none or eigval)"
        )
        get_ad_desc = "Get AD file (a boolean, default is False)"
        get_adc_desc = "Get ADC file (a boolean, default is True)"
        get_cl_desc = "Get CL file (a boolean, default is False)"
        get_cp_desc = "Get CP file (a boolean, default is False)"
        get_cs_desc = "Get CS file (a boolean, default is False)"
        get_value_desc = (
            "Get selected eigenvalue(s) (a boolean, default is False)"
        )
        get_vector_desc = (
            "Get selected eigenvector(s) (a boolean, default is True)"
        )
        get_fa_desc = "Get FA file (a boolean, default is True)"
        get_rd_desc = "Get RD file (a boolean, default is False)"

        # Outputs description
        adc_file_desc = (
            "Output ADC file (a pathlike object or "
            "string representing a file) "
        )
        fa_file_desc = (
            "Output FA file (a pathlike object or "
            "string representing a file) "
        )
        ad_file_desc = (
            "Output AD file (a pathlike object or "
            "string representing a file) "
        )
        rd_file_desc = (
            "Output RD file (a pathlike object or "
            "string representing a file) "
        )
        cl_file_desc = (
            "Output CL file (a pathlike object or "
            "string representing a file) "
        )
        cp_file_desc = (
            "Output CP file (a pathlike object or "
            "string representing a file) "
        )
        cs_file_desc = (
            "Output CS file (a pathlike object or "
            "string representing a file) "
        )
        value_file_desc = (
            "Output selected eigenvalue(s) file (a pathlike object or "
            "string representing a file) "
        )
        vector_file_desc = (
            "Output selected eigenvector(s) file (a pathlike object or "
            "string representing a file) "
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_dti", File(output=False, optional=False, desc=in_dti_desc)
        )

        # Optional inputs traits
        self.add_trait(
            "in_mask", File(output=False, optional=True, desc=in_mask_desc)
        )
        self.add_trait(
            "component",
            List(
                [1],
                output=False,
                optional=True,
                desc=component_desc,
            ),
        )

        self.add_trait(
            "modulate",
            Enum(
                "FA",
                "none",
                "eigval",
                output=False,
                optional=True,
                desc=modulate_desc,
            ),
        )

        self.add_trait(
            "get_adc",
            Bool(
                True,
                output=False,
                optional=True,
                desc=get_adc_desc,
            ),
        )

        self.add_trait(
            "get_fa",
            Bool(
                True,
                output=False,
                optional=True,
                desc=get_fa_desc,
            ),
        )

        self.add_trait(
            "get_ad",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_ad_desc,
            ),
        )

        self.add_trait(
            "get_rd",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_rd_desc,
            ),
        )

        self.add_trait(
            "get_cl",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_cl_desc,
            ),
        )

        self.add_trait(
            "get_cp",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_cp_desc,
            ),
        )

        self.add_trait(
            "get_cl",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_cl_desc,
            ),
        )

        self.add_trait(
            "get_cs",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_cs_desc,
            ),
        )

        self.add_trait(
            "get_value",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_value_desc,
            ),
        )

        self.add_trait(
            "get_vevtor",
            Bool(
                True,
                output=False,
                optional=True,
                desc=get_vector_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "adc_file", File(output=True, optional=True, desc=adc_file_desc)
        )
        self.add_trait(
            "fa_file", File(output=True, optional=True, desc=fa_file_desc)
        )
        self.add_trait(
            "ad_file", File(output=True, optional=True, desc=ad_file_desc)
        )
        self.add_trait(
            "rd_file", File(output=True, optional=True, desc=rd_file_desc)
        )
        self.add_trait(
            "cl_file", File(output=True, optional=True, desc=cl_file_desc)
        )
        self.add_trait(
            "cp_file", File(output=True, optional=True, desc=cp_file_desc)
        )
        self.add_trait(
            "cs_file", File(output=True, optional=True, desc=cs_file_desc)
        )
        self.add_trait(
            "value_file",
            File(output=True, optional=True, desc=value_file_desc),
        )
        self.add_trait(
            "vector_file",
            File(output=True, optional=True, desc=vector_file_desc),
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.TensorMetrics")

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
        super(TensorMetrics, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_dti:
            valid_ext, in_ext, fileName = checkFileExt(self.in_dti, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()
                return

            if self.output_directory:
                if self.get_adc:
                    self.outputs["adc_file"] = self.in_dti.replace(
                        "_dti", "_adc"
                    )
                if self.get_fa:
                    self.outputs["fa_file"] = self.in_dti.replace(
                        "_dti", "_fa"
                    )
                if self.get_ad:
                    self.outputs["ad_file"] = self.in_dti.replace(
                        "_dti", "_ad"
                    )
                if self.get_rd:
                    self.outputs["rd_file"] = self.in_dti.replace(
                        "_dti", "_rd"
                    )
                if self.get_cl:
                    self.outputs["cl_file"] = self.in_dti.replace(
                        "_dti", "_cl"
                    )
                if self.get_cp:
                    self.outputs["cp_file"] = self.in_dti.replace(
                        "_dti", "_cp"
                    )
                if self.get_cs:
                    self.outputs["cs_file"] = self.in_dti.replace(
                        "_dti", "_cs"
                    )
                if self.get_cp:
                    self.outputs["cp_file"] = self.in_dti.replace(
                        "_dti", "_cp"
                    )
                if self.get_value:
                    self.outputs["value_file"] = self.in_dti.replace(
                        "_dti", "_value"
                    )
                if self.get_vector:
                    self.outputs["vector_file"] = self.in_dti.replace(
                        "_dti", "_vector"
                    )

        if self.outputs:
            if self.get_adc:
                self.inheritance_dict[self.outputs["adc_file"]] = self.in_dti
            if self.get_fa:
                self.inheritance_dict[self.outputs["fa_file"]] = self.in_dti
            if self.get_ad:
                self.inheritance_dict[self.outputs["ad_file"]] = self.in_dti
            if self.get_rd:
                self.inheritance_dict[self.outputs["rd_file"]] = self.in_dti
            if self.get_cl:
                self.inheritance_dict[self.outputs["cl_file"]] = self.in_dti
            if self.get_cp:
                self.inheritance_dict[self.outputs["cp_file"]] = self.in_dti
            if self.get_cs:
                self.inheritance_dict[self.outputs["cs_file"]] = self.in_dti
            if self.get_value:
                self.inheritance_dict[self.outputs["value_file"]] = self.in_dti
            if self.get_vector:
                self.inheritance_dict[
                    self.outputs["vector_file"]
                ] = self.in_dti

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(TensorMetrics, self).run_process_mia()
        self.process.in_file = self.in_dti
        self.process.component = self.component
        self.process.modulate = self.modulate
        if self.in_mask:
            self.process.in_mask = self.in_mask
        if self.get_adc:
            self.process.out_adc = self.adc_file
        if self.get_fa:
            self.process.out_fa = self.fa_file
        if self.get_ad:
            self.process.out_ad = self.ad_file
        if self.get_rd:
            self.process.out_rd = self.rd_file
        if self.get_cl:
            self.process.out_cl = self.cl_file
        if self.get_cp:
            self.process.out_cp = self.cp_file
        if self.get_cs:
            self.process.out_cs = self.cs_file
        if self.get_value:
            self.process.out_eval = self.value_file
        if self.get_vector:
            self.process.out_evec = self.vector_file

        return self.process.run(configuration_dict={})


class Tractography(ProcessMIA):
    """
    *Performs streamlines tractography after selecting the appropriate
    algorithm. (tckgen command)*

    Please, see the complete documentation for the `Tractography brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/Tractography.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Tractography, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "Input file to be processed (a pathlike object"
            "string representing an existing file)"
        )
        # Optionnal inputs description
        algorithm_desc = (
            "Tractography algorithm to be used (iFOD2, FACT, iFOD1, "
            "Nulldist, SD_Stream, Tensor_Det or Tensor_Prob)"
        )
        angle_desc = "Set the maximum angle between successive steps (a float)"
        cutoff_desc = (
            "Set the FA or FOD amplitude cutoff for terminating tracks "
        )
        downsample_desc = (
            "Downsample the generated streamlines to reduce output "
            "file size (a float)"
        )
        max_length_desc = " Set the max length of any track in mm (a float)"
        min_length_desc = (
            " Set the minimum length of any track in mm (a float)"
        )
        noprecompt_desc = (
            "Do NOT pre-compute legendre polynomial values (a boolean)"
        )
        select_desc = "Set the desired number of tracks. (an integer)"
        step_size_desc = "Set the step size of the algorithm in mm (a float)"
        trials_desc = (
            "Set the maximum number of sampling trials at each point "
            "(only used for probabilistic tracking) (an integer)"
        )
        stop_desc = (
            "Stop propagating a streamline once it has traversed all include "
            "regions (a boolean)"
        )
        use_rk4_desc = "Use 4th-order Runge-Kutta integration (a boolean)"
        act_file_desc = (
            " Use the Anatomically-Constrained Tractography "
            "framework during tracking (provided image must be "
            "in the 5TT ie five tissue type format(a pathlike object"
            "string representing an existing file)"
        )

        # Outputs description
        out_file_desc = (
            "Output file containing tracks (a pathlike object or "
            "string representing a file) "
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs traits
        self.add_trait(
            "algorithm",
            Enum(
                "iFOD2",
                "FACT",
                "iFOD1",
                "Nulldist",
                "SD_Stream",
                "Tensor_Det",
                "Tensor_Prob",
                output=False,
                optional=True,
                desc=algorithm_desc,
            ),
        )

        self.add_trait(
            "angle",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=angle_desc,
            ),
        )

        self.add_trait(
            "cutoff",
            Float(
                0.1,
                output=False,
                optional=True,
                desc=cutoff_desc,
            ),
        )

        self.add_trait(
            "downsample",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=downsample_desc,
            ),
        )

        self.add_trait(
            "max_length",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=max_length_desc,
            ),
        )

        self.add_trait(
            "min_length",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=min_length_desc,
            ),
        )

        self.add_trait(
            "noprecompt",
            Bool(
                False,
                output=False,
                optional=True,
                desc=noprecompt_desc,
            ),
        )

        self.add_trait(
            "select",
            Int(
                5000,
                output=False,
                optional=True,
                desc=select_desc,
            ),
        )

        self.add_trait(
            "step_size",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=step_size_desc,
            ),
        )

        self.add_trait(
            "stop",
            Bool(
                False,
                output=False,
                optional=True,
                desc=stop_desc,
            ),
        )

        self.add_trait(
            "trials",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=trials_desc,
            ),
        )

        self.add_trait(
            "use_rk4",
            Bool(
                False,
                output=False,
                optional=True,
                desc=use_rk4_desc,
            ),
        )

        self.add_trait(
            "act_file", File(output=False, optional=True, desc=act_file_desc)
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=True, desc=out_file_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.Tractograohy")

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
        super(Tractography, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()
                return

            if self.output_directory:
                self.outputs["out_file"] = ""

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Tractography, self).run_process_mia()

        return self.process.run(configuration_dict={})
