# -*- coding: utf-8 -*-

"""The mrtrix (mrtrix3) preprocess library of the mia_processes package.

The purpose of this module is to customise the main mrtrix preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - ConstrainedSphericalDeconvolution
        - DWIBiasCorrect
        - DWIBrainMask
        - DWIDenoise
        - DWIExtract
        - DWIPreproc
        - EditingTrack
        - FitTensor
        - FilteringTrack
        - Generate5ttfsl
        - Generate5tt2gmwmi
        - MRCat
        - MRConvert
        - MRDeGibbs
        - MRMath
        - MRTransform
        - MTnormalise
        - ResponseSDDhollander
        - SphericalHarmonicExtraction
        - TensorMetrics
        - Tractography
        - TransformFSLConvert
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

# for diffusion, in order to simplify the bricks,
# we only accept mif data (with bvec/bval inside the .mif)
EXT_DWI = {"MIF": "mif"}
EXT = {"NIFTI_GZ": "nii.gz", "NIFTI": "nii", "MIF": "mif"}
EXT_TCK = {"TCK": "tck"}


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
        wm_txt_desc = (
            "WM response text file."
            "(a pathlike object or string representing an existing file)"
        )
        csf_txt_desc = (
            "CSF response text file."
            "(a pathlike object or string representing an existing file)"
        )
        gm_txt_desc = (
            "GM response text file."
            "(a pathlike object or string representing an existing file)"
        )
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
            "Maximum harmonic degree of response function (one by shell) "
            "(a list of items which are an integer) "
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
        gm_odf_desc = (
            "Output GM ODF."
            "(a pathlike object or string representing a file)"
        )
        wm_odf_desc = (
            "Output WM ODF."
            "(a pathlike object or string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "wm_txt", File(output=False, optional=False, desc=wm_txt_desc)
        )

        # Optionnal inouts traits
        self.add_trait(
            "gm_txt", File(output=False, optional=True, desc=gm_txt_desc)
        )

        self.add_trait(
            "csf_txt", File(output=False, optional=True, desc=csf_txt_desc)
        )

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
            "gm_odf", File(output=True, optional=True, desc=gm_odf_desc)
        )

        self.add_trait(
            "wm_odf", File(output=True, optional=False, desc=wm_odf_desc)
        )

        self.add_trait(
            "predicted_signal_file",
            File(output=True, optional=True, desc=predicted_signal_file_desc),
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
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT_DWI)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["csf_odf"] = os.path.join(
                    self.output_directory, fileName + "_csf_odf." + in_ext
                )
                self.outputs["gm_odf"] = os.path.join(
                    self.output_directory, fileName + "_gm_odf." + in_ext
                )
                self.outputs["wm_odf"] = os.path.join(
                    self.output_directory, fileName + "_wm_odf." + in_ext
                )
                if self.get_predicted_signal:
                    self.outputs["predicted_signal_file"] = os.path.join(
                        self.output_directory,
                        fileName + "_odf_predicted_signal." + in_ext,
                    )

        if self.outputs:
            for k in ["csf_odf", "gm_odf", "wm_odf"]:
                self.tags_inheritance(
                    in_file=self.in_file, out_file=self.outputs[k]
                )

            if self.get_predicted_signal:
                self.tags_inheritance(
                    in_file=self.in_file,
                    out_file=self.outputs["predicted_signal_file"],
                )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ConstrainedSphericalDeconvolution, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.algorithm = self.algorithm
        self.process.wm_odf = self.wm_odf
        self.process.wm_txt = self.wm_txt
        if self.gm_txt:
            self.process.gm_odf = self.gm_odf
            self.process.gm_txt = self.gm_txt
        if self.csf_txt:
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
    (dwibiascorrect command)*

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
            "The output bias corrected DWI image (a pathlike object or "
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
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT_DWI)

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
            self.tags_inheritance(
                in_file=self.in_file, out_file=self.outputs["out_file"]
            )

            if self.bias:
                self.tags_inheritance(
                    in_file=self.in_file,
                    out_file=self.outputs["bias_field_map"],
                )

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


class DWIBrainMask(ProcessMIA):
    """
    *Generates a whole brain mask from a DWI image.
    (dwi2mask command)*

    Please, see the complete documentation for the `DWIBrainMask brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/DWIBrainMask.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(DWIBrainMask, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "Input DWI image (a pathlike object"
            "or a string representing an existing file)"
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
        super(DWIBrainMask, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT_DWI)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_brainmask." + in_ext
                )

        if self.outputs:
            self.tags_inheritance(
                in_file=self.in_file, out_file=self.outputs["out_file"]
            )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(DWIBrainMask, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.out_file = self.out_file

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
            Either(
                Undefined,
                Tuple(Int, Int, Int),
                default=Undefined,
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
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT_DWI)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

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
            self.tags_inheritance(
                in_file=self.in_file, out_file=self.outputs["out_file"]
            )

            if self.noise:
                self.tags_inheritance(
                    in_file=self.in_file, out_file=self.outputs["noise_map"]
                )

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
    or certain shells from a DWI dataset (dwiextract command)*

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
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT_DWI)

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
            self.tags_inheritance(
                in_file=self.in_file, out_file=self.outputs["out_file"]
            )

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
    tool if possible (dwifslpreproc command)*

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
            "Keep ALL outputs generated by eddy (including images), and the "
            "output of eddy_qc (if installed)"
            "(a boolean)"
        )
        eddyqc_text_desc = (
            "Keep the various text-based statistical outputs generated "
            "by eddy, and the output of eddy_qc (if installed), "
            "(a boolean)"
        )
        se_epi_corr_desc = (
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
                output=False,
                optional=True,
                desc=rpe_otions_desc,
            ),
        )

        # Optionnal inputs traits
        self.add_trait(
            "se_epi_corr",
            File(output=False, optional=True, desc=se_epi_corr_desc),
        )

        self.add_trait(
            "pe_dir",
            Enum(
                "ap",
                "pa",
                "lr",
                "rl",
                output=False,
                optional=True,
                desc=pe_dir_desc,
            ),
        )

        self.add_trait(
            "ro_time",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=ro_time_desc,
            ),
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
                desc=eddyqc_text_desc,
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
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT_DWI)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory,
                    fileName + "_dwifslpreproc." + in_ext,
                )

        if self.outputs:
            self.tags_inheritance(
                in_file=self.in_file, out_file=self.outputs["out_file"]
            )

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


class EditingTrack(ProcessMIA):
    """
    *Perform various editing operations on track files.
    (tckedit command)*

    Please, see the complete documentation for the `EditingTrack brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/EditingTrack.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(EditingTrack, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["mrtrix"]

        # Mandatory inputs description
        in_tracks_desc = (
            "Inputs track file(s) (a list of items which are a pathlike "
            "object or a string representing an existing file)"
        )
        suffix_desc = "Output file suffix (a string)"
        roi_excl_desc = (
            "Specify an exclusion region of interest, streamlines that enter "
            "ANY exclude region will be discarded. (a pathlike "
            "object or string representing an existing file or a tuple of the "
            "form: (a float, a float, a float, a float))"
        )
        roi_incl_desc = (
            "Specify an inclusion region of interest, streamlines must "
            "traverse ALL inclusion regions to be accepted. (a pathlike "
            "object or string representing an existing file or a tuple of the "
            "form: (a float, a float, a float, a float))"
        )
        roi_incl_ordered_desc = (
            "Specify an inclusion region of interest, streamlines must "
            "traverse ALL inclusion_ordered regions in the order "
            "they are specified in order to be accepted. (a pathlike "
            "object or string representing an existing file or a tuple of the "
            "form: (a float, a float, a float, a float))"
        )
        roi_mask_desc = (
            "Specify a masking region of interest. If defined,streamlines "
            "exiting the mask will be truncated. (a pathlike "
            "object or string representing an existing file or a tuple of the "
            "form: (a float, a float, a float, a float))"
        )
        maxlength_desc = "The maximum length of any streamline in mm (a float)"
        minlength_desc = "The minimum length of any streamline in mm (a float)"
        number_desc = (
            "The desired number of selected streamlines to be propagated to "
            "the output file (an integer)"
        )
        skip_desc = (
            "Omit this number of selected streamlines before commencing "
            "writing to the output file (an integer)"
        )
        maxweight_desc = "The maximum weight of any streamline (an integer)"
        minweight_desc = "The minimum weight of any streamline (an integer)"
        inverse_desc = (
            "Output the inverse selection of streamlines based on the "
            "criteria provided; i.e. only those streamlines that fail at "
            "least one selection criterion, and/or vertices that are outside "
            "masks if provided, will be written to file (a boolean)"
        )
        ends_only_desc = (
            "Only test the ends of each streamline against the provided "
            "include/exclude ROIs (a boolean)"
        )
        tck_weights_in_desc = (
            "Specify a text scalar file containing the streamline weights"
            "(a pathlike object or string representing an existing file)"
        )
        get_tck_weights_out_desc = (
            "Get  an output text scalar file containing streamline weights "
            "(a boolean)"
        )
        # Output descriptions
        tracks_out_desc = (
            "The output track file (a pathlike object or "
            "string representing a file)"
        )
        tck_weights_out_desc = (
            "Output text scalar file containing streamline weights"
            "(a pathlike object or string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_tracks",
            InputMultiPath(
                Either(File(), List(File())),
                output=False,
                desc=in_tracks_desc,
            ),
        )

        self.add_trait(
            "suffix",
            String(
                "edited",
                output=False,
                optional=True,
                desc=suffix_desc,
            ),
        )

        self.add_trait(
            "roi_excl",
            Either(
                File(),
                Tuple(Float, Float, Float, Float),
                default=Undefined,
                output=False,
                optional=True,
                desc=roi_excl_desc,
            ),
        )

        self.add_trait(
            "roi_incl",
            Either(
                File(),
                Tuple(Float, Float, Float, Float),
                output=False,
                optional=True,
                desc=roi_incl_desc,
            ),
        )

        self.add_trait(
            "roi_incl_ordered",
            Either(
                File(),
                Tuple(Float, Float, Float, Float),
                output=False,
                optional=True,
                desc=roi_incl_ordered_desc,
            ),
        )

        self.add_trait(
            "roi_mask",
            Either(
                File(),
                Tuple(Float, Float, Float, Float),
                output=False,
                optional=True,
                desc=roi_mask_desc,
            ),
        )

        self.add_trait(
            "maxlenght",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=maxlength_desc,
            ),
        )

        self.add_trait(
            "minlenght",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=minlength_desc,
            ),
        )

        self.add_trait(
            "number",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=number_desc,
            ),
        )

        self.add_trait(
            "skip",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=skip_desc,
            ),
        )

        self.add_trait(
            "maxweight",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=maxweight_desc,
            ),
        )

        self.add_trait(
            "minweight",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=minweight_desc,
            ),
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
            "ends_only",
            Bool(
                False,
                output=False,
                optional=True,
                desc=ends_only_desc,
            ),
        )

        self.add_trait(
            "tck_weights_in",
            File(output=False, optional=True, desc=tck_weights_in_desc),
        )

        self.add_trait(
            "get_tck_weights_out",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_tck_weights_out_desc,
            ),
        )

        # Outputs
        self.add_trait(
            "tracks_out",
            File(output=True, optional=False, desc=tracks_out_desc),
        )

        self.add_trait(
            "tck_weights_out",
            File(output=True, optional=True, desc=tck_weights_out_desc),
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
        super(EditingTrack, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_tracks:
            for in_file in self.in_tracks:
                valid_ext, in_ext, fileName = checkFileExt(in_file, EXT_TCK)
                if not valid_ext:
                    print("\nThe input image format is not recognized...!")
                    return self.make_initResult()

            if self.output_directory:
                # Name of the first track used for output name
                valid_ext, in_ext, fileName = checkFileExt(
                    self.in_tracks[0], EXT
                )
                self.outputs["tracks_out"] = os.path.join(
                    self.output_directory,
                    fileName + "_" + self.suffix + "." + in_ext,
                )

                self.outputs["tck_weights_out"] = os.path.join(
                    self.output_directory, fileName + "_weight.txt"
                )

        if self.outputs:
            # First track file used for inheritance
            self.tags_inheritance(
                in_file=self.in_tracks[0], out_file=self.outputs["tracks_out"]
            )
            self.tags_inheritance(
                in_file=self.in_tracks[0],
                out_file=self.outputs["tck_weights_out"],
            )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(EditingTrack, self).run_process_mia()
        cmd = ["tckedit"]

        for in_file in self.in_tracks:
            cmd += [in_file]

        cmd += [self.tracks_out]

        if self.roi_excl:
            cmd += ["-exclude", self.roi_excl]
        if self.roi_incl:
            cmd += ["-include", self.roi_incl]
        if self.roi_incl_ordered:
            cmd += ["-include_ordered", self.roi_incl_ordered]
        if self.roi_mask:
            cmd += ["-mask", self.roi_mask]
        if self.maxlenght:
            cmd += ["-maxlenght", self.maxlenght]
        if self.minlength:
            cmd += ["-minlength", self.minlength]
        if self.number:
            cmd += ["-number", self.number]
        if self.skip:
            cmd += ["-skip", self.skip]
        if self.maxweight:
            cmd += ["-maxweight", self.maxweight]
        if self.minweight:
            cmd += ["-minweight", self.minweight]
        if self.inverse:
            cmd += ["-inverse"]
        if self.ends_only:
            cmd += ["-minweight"]
        if self.tck_weights_in:
            cmd += ["-tck_weights_in", self.tck_weights_in]
        if self.get_tck_weights_out:
            cmd += ["-tck_weights_out", self.tck_weights_out]

        return mrtrix.mrtrix_call(cmd)


class FilteringTrack(ProcessMIA):
    """
    *Filter a whole-brain fibre-tracking data set such that the streamline
    densities match the FOD lobe integrals.
    (tcksift command)*

    Please, see the complete documentation for the `FilteringTrack brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/FilteringTrack.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(FilteringTrack, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["mrtrix"]

        # Mandatory inputs description
        in_tracks_desc = (
            "Input track file (a pathlike "
            "object or a string representing an existing file)"
        )
        suffix_desc = "Output file suffix (a string)"
        in_fod_desc = (
            "Input image containing the spherical harmonics of the fibre "
            "orientation distributions (a pathlike "
            "object or a string representing an existing file)"
        )
        proc_mask_desc = (
            "Provide an image containing the processing mask weights for "
            "the model (a pathlike object or a string representing an "
            "existing file)"
        )
        act_image_desc = (
            "Provide an ACT five-tissue-type segmented anatomical image "
            "to derive the processing mask (a pathlike object or a string "
            "representing an existing file)"
        )
        fd_scale_gm_desc = (
            "Provide this option (in conjunction with -act) to heuristically "
            "downsize the fibre density estimates based on the presence of GM "
            "in the vox (a boolean)"
        )
        no_dilate_lut_desc = (
            "Do NOT dilate FOD lobe lookup tables; only map streamlines to "
            "FOD lobes if the precise tangent lies within the angular spread "
            "of that lobe (a boolean)"
        )
        make_null_lobes_desc = (
            "Add an additional FOD lobe to each voxel, with zero integral, "
            "that covers all directions with zero / negative FOD amplitudes "
            "(a boolean)"
        )
        remove_untracked_desc = (
            "Remove FOD lobes that do not have any streamline density "
            "attributed to them (a boolean)"
        )
        fd_thresh_value_desc = "Fibre density threshold (an integer)"
        term_number_value_desc = (
            "Number of streamlines - continue filtering until this number "
            "of streamlines remain (an integer)"
        )
        term_ratio_value_desc = (
            "Termination ratio - defined as the ratio between reduction in "
            "cost function, and reduction in density of streamlines. "
            "(an integer)"
        )
        term_mu_value_desc = (
            "Terminate filtering once the SIFT proportionality coefficient "
            "reaches a given value (an integer)"
        )
        get_csv_file_desc = (
            "Output statistics of execution per iteration to a .csv file "
            "(a boolean)"
        )
        get_mu_file_desc = (
            "Output the final value of SIFT proportionality coefficient mu "
            "to a text file (a boolean)"
        )
        get_out_selection_file_desc = (
            "Output a text file containing the binary selection of streamlines"
            "(a boolean)"
        )
        # Output descriptions
        tracks_out_desc = (
            "The output filtered tracks file (a pathlike object or "
            "string representing a file)"
        )
        csv_file_out_desc = (
            "A csv file with the output statistics of execution per iteration "
            "(a pathlike object or string representing a file)"
        )
        mu_file_out_desc = (
            "The final value of SIFT proportionality coefficient mu in a text "
            "file (a pathlike object or string representing a file)"
        )
        selection_file_out_desc = (
            "A text file containing the binary selection of streamlines "
            "(a pathlike object or string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_tracks",
            InputMultiPath(
                Either(File(), List(File())),
                output=False,
                desc=in_tracks_desc,
            ),
        )

        self.add_trait(
            "in_fod",
            File(
                output=False,
                desc=in_fod_desc,
            ),
        )

        # Optional inputs traits
        self.add_trait(
            "suffix",
            String(
                "sift",
                output=False,
                optional=True,
                desc=suffix_desc,
            ),
        )

        self.add_trait(
            "proc_mask",
            File(
                output=False,
                optional=True,
                desc=proc_mask_desc,
            ),
        )

        self.add_trait(
            "act_image",
            File(
                output=False,
                optional=True,
                desc=act_image_desc,
            ),
        )

        self.add_trait(
            "fd_scale_gm",
            Bool(
                False,
                output=False,
                optional=True,
                desc=fd_scale_gm_desc,
            ),
        )

        self.add_trait(
            "no_dilate_lut",
            Bool(
                False,
                output=False,
                optional=True,
                desc=no_dilate_lut_desc,
            ),
        )

        self.add_trait(
            "make_null_lobes",
            Bool(
                False,
                output=False,
                optional=True,
                desc=make_null_lobes_desc,
            ),
        )

        self.add_trait(
            "remove_untracked",
            Bool(
                False,
                output=False,
                optional=True,
                desc=remove_untracked_desc,
            ),
        )

        self.add_trait(
            "fd_thresh_value",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=fd_thresh_value_desc,
            ),
        )

        self.add_trait(
            "term_number_value",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=term_number_value_desc,
            ),
        )

        self.add_trait(
            "term_ratio_value",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=term_ratio_value_desc,
            ),
        )

        self.add_trait(
            "term_mu_value",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=term_mu_value_desc,
            ),
        )

        self.add_trait(
            "get_csv_file",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_csv_file_desc,
            ),
        )

        self.add_trait(
            "get_mu_file",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_mu_file_desc,
            ),
        )

        self.add_trait(
            "get_out_selection_file",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_out_selection_file_desc,
            ),
        )

        # Outputs
        self.add_trait(
            "tracks_out",
            File(output=True, optional=False, desc=tracks_out_desc),
        )

        self.add_trait(
            "csv_file_out",
            File(output=True, optional=True, desc=csv_file_out_desc),
        )

        self.add_trait(
            "mu_file_out",
            File(output=True, optional=True, desc=mu_file_out_desc),
        )

        self.add_trait(
            "selection_file_out",
            File(output=True, optional=True, desc=selection_file_out_desc),
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
        super(FilteringTrack, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_tracks:
            valid_ext, in_ext, fileName = checkFileExt(self.in_tracks, EXT_TCK)
            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["tracks_out"] = os.path.join(
                    self.output_directory,
                    fileName + "_" + self.suffix + "." + in_ext,
                )
                if self.get_csv_file:
                    self.outputs["csv_file_out"] = os.path.join(
                        self.output_directory, fileName + "_tcksift_stats.csv"
                    )
                if self.get_mu_file:
                    self.outputs["mu_file_out"] = os.path.join(
                        self.output_directory, fileName + "_tcksift_mu.txt"
                    )
                if self.get_out_selection_file:
                    self.outputs["selection_file_out"] = os.path.join(
                        self.output_directory,
                        fileName + "_tcksift_selection.txt",
                    )

        if self.outputs:
            self.tags_inheritance(
                in_file=self.in_tracks, out_file=self.outputs["tracks_out"]
            )
            if self.get_csv_file:
                self.tags_inheritance(
                    in_file=self.in_tracks,
                    out_file=self.outputs["csv_file_out"],
                )
            if self.get_mu_file:
                self.tags_inheritance(
                    in_file=self.in_tracks,
                    out_file=self.outputs["mu_file_out"],
                )
            if self.get_out_selection_file:
                self.tags_inheritance(
                    in_file=self.in_tracks,
                    out_file=self.outputs["selection_file_out"],
                )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(FilteringTrack, self).run_process_mia()
        cmd = ["tcksift"]

        if self.proc_mask:
            cmd += ["-proc_mask", self.proc_mask]
        if self.act_image:
            cmd += ["-act", self.act_image]
        if self.fd_scale_gm:
            cmd += ["-fd_scale_gm"]
        if self.no_dilate_lut:
            cmd += ["-no_dilate_lut"]
        if self.make_null_lobes:
            cmd += ["-make_null_lobes"]
        if self.remove_untracked:
            cmd += ["-remove_untracked"]
        if self.fd_thresh_value:
            cmd += ["-fd_thresh", self.fd_thresh_value]
        if self.term_number_value:
            cmd += ["-term_number", self.term_number_value]
        if self.term_ratio_value:
            cmd += ["-term_ratio", self.term_ratio_value]
        if self.term_mu_value:
            cmd += ["-term_mu", self.term_mu_value]
        if self.get_csv_file:
            cmd += ["-csv", self.csv_file_out]
        if self.get_mu_file:
            cmd += ["-out_mu", self.mu_file_out]
        if self.get_out_selection_file:
            cmd += ["-out_selection", self.selection_file_out]

        cmd += [self.in_tracks, self.in_fod, self.tracks_out]

        return mrtrix.mrtrix_call(cmd)


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

        # Nipype FitTransform process is not up to date for some option
        # so the "args" input is used for several options

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
        estimate_dkt_desc = "Estimate diffusion kurtosis (a boolean)"
        get_predicted_signal_desc = (
            "Get a file with the predicted signal from the tensor fits "
            "(a boolean) "
        )
        get_output_b0_desc = "Get the output b0 (a boolean)"
        ols_option_desc = (
            "Perform initial fit using an ordinary least-squares (OLS) fit"
            "(a boolean)"
        )
        number_of_iter_desc = (
            "Number of iterative reweightings for IWLS algorithm"
            "(an integer, default value is 2)"
        )
        # Output descriptions
        out_file_desc = (
            "Output tensor image (a pathlike object or "
            "string representing a file)"
        )
        out_dkt_desc = (
            "Out diffusion kurtosis (a pathlike object or "
            "string representing a file)"
        )
        out_b0_desc = (
            "Out b0 image (a pathlike object or " "string representing a file)"
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
            "estimate_dkt",
            Bool(
                False,
                output=False,
                optional=True,
                desc=estimate_dkt_desc,
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

        self.add_trait(
            "get_output_b0",
            Bool(
                False,
                output=False,
                optional=True,
                desc=get_output_b0_desc,
            ),
        )

        self.add_trait(
            "ols_option",
            Bool(
                False,
                output=False,
                optional=True,
                desc=ols_option_desc,
            ),
        )

        self.add_trait(
            "number_of_iter",
            Int(
                2,
                output=False,
                optional=True,
                desc=number_of_iter_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )

        self.add_trait(
            "out_dkt",
            File(output=True, optional=True, desc=out_dkt_desc),
        )

        self.add_trait(
            "out_b0",
            File(output=True, optional=True, desc=out_b0_desc),
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
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT_DWI)

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
                if self.estimate_dkt:
                    self.outputs["out_dkt"] = os.path.join(
                        self.output_directory,
                        fileName + "_dki." + in_ext,
                    )
                if self.get_output_b0:
                    self.outputs["out_b0"] = os.path.join(
                        self.output_directory,
                        fileName + "_b0." + in_ext,
                    )

        if self.outputs:
            self.tags_inheritance(
                in_file=self.in_file, out_file=self.outputs["out_file"]
            )
            if self.get_predicted_signal:
                self.tags_inheritance(
                    in_file=self.in_file,
                    out_file=self.outputs["predicted_signal_file"],
                )
            if self.estimate_dkt:
                self.tags_inheritance(
                    in_file=self.in_file,
                    out_file=self.outputs["out_dkt"],
                )
            if self.get_output_b0:
                self.tags_inheritance(
                    in_file=self.in_file,
                    out_file=self.outputs["out_b0"],
                )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(FitTensor, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.out_file = self.out_file
        if self.in_mask:
            self.process.in_mask = self.in_mask
        if self.get_predicted_signal:
            self.process.predicted_signal = self.predicted_signal_file
        args = ""
        if self.ols_option:
            args += "-ols "
        if self.number_of_iter:
            args += "-iter " + self.number_of_iter + " "
        if self.estimate_dkt:
            args += "-dkt " + self.out_dkt + " "
        if self.get_output_b0:
            args += "-b0 " + self.out_b0 + " "
        if args:
            self.process.args = args

        return self.process.run(configuration_dict={})


class Generate5ttfsl(ProcessMIA):
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
        super(Generate5ttfsl, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix", "fsl"]

        # Mandatory inputs description
        in_file_desc = (
            "Input T1 image (a pathlike object"
            "string representing an existing file)"
        )
        # Optionnal inputs description
        in_mask_desc = (
            "Manually provide a brain mask, rather than deriving one in the "
            "script (a pathlike object string representing an existing file)"
        )
        t2_image_desc = (
            "Provide a T2-weighted image in addition to T1w image"
            "(a pathlike object string representing an existing file)"
        )
        no_crop_desc = (
            "Do NOT crop the resulting 5TT image to reduce its size "
            "(keep the same dimensions as the input image)"
            "(a boolean)"
        )
        sgm_amyg_hipp_desc = (
            "Represent the amygdalae and hippocampi as sub-cortical grey "
            "matter in the 5TT image (a boolean)"
        )
        premasked_desc = (
            "Indicate that brain masking has already been applied to the "
            "input image(a boolean)"
        )
        # Output descriptions
        out_file_desc = (
            "Output image (a pathlike object or string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs
        self.add_trait(
            "in_mask", File(output=False, optional=True, desc=in_mask_desc)
        )

        self.add_trait(
            "t2_image", File(output=False, optional=True, desc=t2_image_desc)
        )

        self.add_trait(
            "no_crop",
            Bool(
                False,
                output=False,
                optional=True,
                desc=no_crop_desc,
            ),
        )

        self.add_trait(
            "sgm_amyg_hipp",
            Bool(
                False,
                output=False,
                optional=True,
                desc=sgm_amyg_hipp_desc,
            ),
        )

        self.add_trait(
            "premasked",
            Bool(
                False,
                output=False,
                optional=True,
                desc=premasked_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.Generate5tt")

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
        super(Generate5ttfsl, self).list_outputs()

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
            self.tags_inheritance(
                in_file=self.in_file, out_file=self.outputs["out_file"]
            )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Generate5ttfsl, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.algorithm = "fsl"
        self.process.out_file = self.out_file
        args = ""
        if self.in_mask:
            args += "-mask " + self.in_mask + " "
        if self.t2_image:
            args += "-t2 " + self.t2_image + " "
        if self.premasked:
            args += "-premasked "
        if self.no_crop:
            args += "-nocrop "
        if self.sgm_amyg_hipp:
            args += "-sgm_amyg_hipp "
        if args:
            self.process.args = args

        return self.process.run(configuration_dict={})


class Generate5tt2gmwmi(ProcessMIA):
    """
    *Generate a mask image appropriate for seeding streamlines
    on the grey matter-white matter interface (5tt2gmwmi command)*

    Please, see the complete documentation for the `Generate5tt2gmwmi brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/Generate5tt2gmwmi.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Generate5tt2gmwmi, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix", "fsl", "freesurfer"]

        # Mandatory inputs description
        in_file_desc = (
            "The input 5TT segmented anatomical image (a pathlike object"
            "string representing an existing file)"
        )
        # Optionnal inputs description
        in_mask_desc = (
            "Filter an input mask image according to those voxels that lie "
            "upon the grey matter - white matter boundary (a pathlike object"
            "string representing an existing file)"
        )
        # Outputs description
        out_file_desc = (
            "The output mask image (a pathlike object or "
            "string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs
        self.add_trait(
            "in_mask", File(output=False, optional=True, desc=in_mask_desc)
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )

        # No nipype command for this process
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
        super(Generate5tt2gmwmi, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_gmwmSeed." + in_ext
                )

        if self.outputs:
            self.tags_inheritance(
                in_file=self.in_file, out_file=self.outputs["out_file"]
            )
        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Generate5tt2gmwmi, self).run_process_mia()
        cmd = ["5tt2gmwmi"]
        if self.in_mask:
            cmd += ["-mask_in", self.in_mask]
        cmd += [self.in_file, self.out_file]

        return mrtrix.mrtrix_call(cmd)


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
            self.tags_inheritance(
                in_file=self.in_files[0], out_file=self.outputs["out_file"]
            )

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
            "(a list of items which are a tuple of the form (axis, selection))"
        )
        scaling_desc = (
            "Specify the data scaling parameters used to rescale the "
            "intensity values(a list of items which are a float)"
        )
        vox_desc = (
            "Change the voxel dimensions reported in the output image header"
            "(a list of items which are a float or an integer)"
        )
        out_file_format_desc = (
            "Format of the output image (NIFTI, NIFTI_GZ or MIF)"
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
            "(yes or no)"
        )
        grad_file_desc = (
            "Provide the diffusion-weighted gradient scheme used in the "
            "acquisition in a text file (MRTrix format) (a pathlike object "
            "or string representing an existing file) "
        )
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
                Tuple(Int(), String()),
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
                List(),
                default=Undefined,
                output=False,
                optional=True,
                desc=vox_desc,
            ),
        )

        self.add_trait(
            "export_bvec_bval",
            Bool(
                False, output=False, optional=True, desc=export_bvec_bval_desc
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
            Either(
                Undefined,
                Enum("no", "yes"),
                output=False,
                optional=True,
                desc=bval_scale_desc,
            ),
        )

        self.add_trait(
            "grad_file", File(output=False, optional=True, desc=grad_file_desc)
        )

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
            self.tags_inheritance(
                in_file=self.in_file, out_file=self.outputs["out_file"]
            )
            if self.export_bvec_bval:
                self.tags_inheritance(
                    in_file=self.in_file, out_file=self.outputs["out_bvec"]
                )
                self.tags_inheritance(
                    in_file=self.in_file, out_file=self.outputs["out_bval"]
                )

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
            axes = ""
            for i in self.axes:
                axes += str(i) + ","
            cmd += ["-axes", axes[:-1]]
        if self.coord:
            axis, selection = self.coord
            coord = str(axis) + " " + str(selection)
            cmd += ["-coord", coord]
        if self.scaling:
            scaling = ""
            for i in self.scaling:
                scaling += str(i) + ","
            cmd += ["-scaling", scaling[:-1]]
        if self.vox:
            vox = ""
            for i in self.vox:
                vox += str(i) + ","
            cmd += ["-vox", vox[:-1]]
        if self.grad_file:
            cmd += ["-grad", self.grad_file]
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
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT_DWI)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_unringed." + in_ext
                )

        if self.outputs:
            self.tags_inheritance(
                in_file=self.in_file, out_file=self.outputs["out_file"]
            )

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
    specified axis of a single image (mrmath command)*

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
                optional=True,
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
            self.tags_inheritance(
                in_file=self.in_file, out_file=self.outputs["out_file"]
            )

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
            "Invert the specified transform before using it (a boolean)"
        )
        flip_axes_desc = (
            "Flip the specified axes (a list of int with 0:x, 1:y and 2:z)"
        )
        half_desc = (
            "Apply the matrix square root of the transformation (a boolean)"
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
            "Specify whether to perform FOD reorientation (a boolean)"
        )

        # Outputs description
        out_file_desc = (
            "The output image of the transformation."
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

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_transformed." + in_ext
                )

        if self.outputs:
            self.tags_inheritance(
                in_file=self.in_file, out_file=self.outputs["out_file"]
            )

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
            flip = ""
            for i in self.flip_axes:
                flip += str(i) + ","
            args += "-flip " + flip[:-1]
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
            "the output images (a boolean)"
        )

        # Outputs description
        out_files_desc = (
            "Nomalised outputs images (a list of items which are a pathlike "
            "object or string representing a file)"
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
                output=True,
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
                valid_ext, in_ext, fileName = checkFileExt(f, EXT_DWI)

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
            self.tags_inheritance(
                in_file=self.in_files[0], out_file=self.outputs["out_files"]
            )

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
        fa_thresh_desc = (
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
            "fa_thresh",
            Float(
                0.2,
                output=False,
                optional=True,
                desc=fa_thresh_desc,
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
            "csf_file", File(output=True, optional=True, desc=csf_file_desc)
        )
        self.add_trait(
            "gm_file", File(output=True, optional=True, desc=gm_file_desc)
        )
        self.add_trait(
            "wm_file", File(output=True, optional=True, desc=wm_file_desc)
        )
        self.add_trait(
            "voxels_image",
            File(output=True, optional=True, desc=voxels_image_desc),
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
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT_DWI)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["csf_file"] = os.path.join(
                    self.output_directory, fileName + "_response_csf.txt"
                )
                self.outputs["gm_file"] = os.path.join(
                    self.output_directory, fileName + "_response_gm.txt"
                )
                self.outputs["wm_file"] = os.path.join(
                    self.output_directory, fileName + "_response_wm.txt"
                )
                if self.get_final_voxels:
                    self.outputs["voxels_image"] = os.path.join(
                        self.output_directory,
                        fileName + "_response_voxels." + in_ext,
                    )

        if self.outputs:
            for k in ["csf_file", "gm_file", "wm_file"]:
                self.tags_inheritance(
                    in_file=self.in_file, out_file=self.outputs[k]
                )
            if self.get_final_voxels:
                self.tags_inheritance(
                    in_file=self.in_file, out_file=self.outputs["voxels_image"]
                )

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
        if self.fa_thresh != 0.2:
            args += "-fa " + self.fa_thresh + " "
        if self.wm_algo:
            args += "-wm_algo " + self.wm_algo + " "
        if self.get_final_voxels:
            args += "-voxels " + self.voxels_image + " "
        if args:
            self.process.args = args

        return self.process.run(configuration_dict={})


class SphericalHarmonicExtraction(ProcessMIA):
    """
    *Extract the peaks of a spherical harmonic function in each voxel.
    (sh2peaks command)*

    Please, see the complete documentation for the
    `SphericalHarmonicExtraction brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/SphericalHarmonicExtraction.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(SphericalHarmonicExtraction, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["mrtrix"]

        # Mandatory inputs description
        in_SH_coeff_desc = (
            "The input image of SH coefficients (a pathlike "
            "object or a string representing an existing file)"
        )
        suffix_desc = "Output file suffix (a string)"
        num_desc = "The number of peaks to extract (an integer)"
        direction_desc = (
            "The direction of a peak to estimate (phi, theta) "
            "(a tuple of the form (a Float, a Float)) "
        )
        peaks_image_desc = (
            "The program will try to find the peaks that most closely match "
            "those in the image provided.  (a pathlike "
            "object or a string representing an existing file)"
        )
        thresh_value_desc = (
            "Only peak amplitudes greater than the threshold will be "
            "considered. (a float)"
        )
        seeds_file_desc = (
            "Specify a set of directions from which to start the multiple "
            "restarts of the optimisation (a pathlike object or a string "
            "representing an existing file)"
        )
        mask_image_desc = (
            "Only perform computation within the specified binary brain "
            "mask image. (a pathlike object or a string representing "
            "an existing file)"
        )
        fast_desc = (
            "Use lookup table to compute associated Legendre polynomials "
            "(faster, but approximate). (a boolean)"
        )
        # Output descriptions
        output_image_desc = (
            "The output image. Each volume corresponds to the x, y & z "
            "component of each peak direction vector in turn. "
            "(a pathlike object or string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_SH_coeff",
            File(
                output=False,
                desc=in_SH_coeff_desc,
            ),
        )

        # Optional inputs traits
        self.add_trait(
            "suffix",
            String(
                "peaks",
                output=False,
                optional=True,
                desc=suffix_desc,
            ),
        )

        self.add_trait(
            "num",
            Int(
                3,
                output=False,
                optional=True,
                desc=num_desc,
            ),
        )

        self.add_trait(
            "direction",
            Either(
                Undefined,
                Tuple(Float, Float),
                default=Undefined,
                output=False,
                optional=True,
                desc=direction_desc,
            ),
        )

        self.add_trait(
            "peaks_image",
            File(output=False, optional=False, desc=peaks_image_desc),
        )

        self.add_trait(
            "thresh_value",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=thresh_value_desc,
            ),
        )

        self.add_trait(
            "seeds_file",
            File(output=False, optional=False, desc=seeds_file_desc),
        )

        self.add_trait(
            "mask_image",
            File(output=False, optional=False, desc=mask_image_desc),
        )

        self.add_trait(
            "fast",
            Bool(
                False,
                output=False,
                optional=True,
                desc=fast_desc,
            ),
        )

        # Outputs
        self.add_trait(
            "output_image",
            File(output=True, optional=False, desc=output_image_desc),
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
        super(SphericalHarmonicExtraction, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_tracks:
            valid_ext, in_ext, fileName = checkFileExt(
                self.in_SH_coeff, EXT_DWI
            )
            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                self.outputs["output_image"] = os.path.join(
                    self.output_directory,
                    fileName + "_" + self.suffix + "." + in_ext,
                )

        if self.outputs:
            self.tags_inheritance(
                in_file=self.in_SH_coeff, out_file=self.outputs["output_image"]
            )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(SphericalHarmonicExtraction, self).run_process_mia()
        cmd = ["sh2peaks"]

        if self.num:
            cmd += ["-num", self.num]
        if self.direction:
            direction = str(self.direction[0]) + str(self.direction[1])
            cmd += ["-direction", direction]
        if self.peaks_image:
            cmd += ["-peaks", self.peaks_image]
        if self.thresh_value:
            cmd += ["-threshold", self.thresh_value]
        if self.seeds_file:
            cmd += ["-seeds", self.seeds_file]
        if self.mask_image:
            cmd += ["-mask", self.mask_image]
        if self.fast:
            cmd += ["-fast"]

        cmd += [self.in_SH_coeff, self.output_image]

        return mrtrix.mrtrix_call(cmd)


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
            valid_ext, in_ext, fileName = checkFileExt(self.in_dti, EXT_DWI)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

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
                self.tags_inheritance(
                    in_file=self.in_dti, out_file=self.outputs["adc_file"]
                )
            if self.get_fa:
                self.tags_inheritance(
                    in_file=self.in_dti, out_file=self.outputs["fa_file"]
                )
            if self.get_ad:
                self.tags_inheritance(
                    in_file=self.in_dti, out_file=self.outputs["ad_file"]
                )
            if self.get_rd:
                self.tags_inheritance(
                    in_file=self.in_dti, out_file=self.outputs["rd_file"]
                )
            if self.get_cl:
                self.tags_inheritance(
                    in_file=self.in_dti, out_file=self.outputs["cl_file"]
                )
            if self.get_cp:
                self.tags_inheritance(
                    in_file=self.in_dti, out_file=self.outputs["cp_file"]
                )
            if self.get_cs:
                self.tags_inheritance(
                    in_file=self.in_dti, out_file=self.outputs["cs_file"]
                )
            if self.get_value:
                self.tags_inheritance(
                    in_file=self.in_dti, out_file=self.outputs["value_file"]
                )
            if self.get_vector:
                self.tags_inheritance(
                    in_file=self.in_dti, out_file=self.outputs["vector_file"]
                )

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
        downsample_factor_desc = (
            "Downsample the generated streamlines to reduce output "
            "file size (a float)"
        )
        max_length_desc = "Set the max length of any track in mm (a float)"
        min_length_desc = "Set the minimum length of any track in mm (a float)"
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
        seed_dynamic_desc = (
            "Determine seed points dynamically using the SIFT model "
            "(must not provide any other seeding mechanism)."
            "(a pathlike object string representing an existing file)"
        )
        seed_gmwmi_desc = (
            "Seed from the grey matter - white matter interface "
            "(only valid if using ACT framework)."
            "(a pathlike object string representing an existing file)"
        )
        seed_grid_voxel_desc = (
            "Seed a fixed number of streamlines per voxel in a mask image; "
            "place seeds on a 3D mesh grid (a tuple of the form: (a pathlike "
            "object or string representing an existing file, an integer)) "
        )
        seed_image_desc = (
            "Seed streamlines entirely at random within a mask image"
            "(a pathlike object string representing an existing file)"
        )
        seed_rejection_desc = (
            "Seed from an image using rejection sampling "
            "(a pathlike object string representing an existing file)"
        )
        seed_rnd_voxel_desc = (
            "Seed a fixed number of streamlines per voxel in a mask image; "
            "random placement of seeds in each voxel."
            "(a tuple of the form: (a pathlike object or string representing "
            "an existing file, an integer)) "
        )
        seed_sphere_desc = (
            "Spherical seed. (a tuple of the form: (a float, a float, "
            "a float, a float)) "
        )
        tracto_seeds_number_desc = (
            "Set the number of seeds that tckgen will attempt to track from. "
            "(an integer)"
        )
        tracto_max_attempts_per_seed_number_desc = (
            "Set the maximum number of times that the tracking algorithm "
            "should attempt to find an appropriate tracking direction from a "
            "given seed point. (an integer, default is 1000) "
        )
        tracto_seed_cutoff_desc = (
            "Set the minimum FA or FOD amplitude for seeding tracks "
            "(an integer, default is 0.1)."
        )
        tracto_seed_unidirectional_desc = (
            "Track from the seed point in one direction only" "(a boolean)"
        )
        tracto_seed_direction_desc = (
            "Specify a seeding direction for the tracking "
            "(a tuple of the form: (a float, a float, a float))"
        )
        tracto_get_output_seeds_desc = (
            "Get the seed location of all successful streamlines into a file"
            "(a boolean)"
        )
        roi_excl_desc = (
            "Specify an exclusion region of interest, streamlines that enter "
            "ANY exclude region will be discarded. (a pathlike "
            "object or string representing an existing file or a tuple of the "
            "form: (a float, a float, a float, a float))"
        )
        roi_incl_desc = (
            "Specify an inclusion region of interest, streamlines must "
            "traverse ALL inclusion regions to be accepted. (a pathlike "
            "object or string representing an existing file or a tuple of the "
            "form: (a float, a float, a float, a float))"
        )
        roi_incl_ordered_desc = (
            "Specify an inclusion region of interest, streamlines must "
            "traverse ALL inclusion_ordered regions in the order "
            "they are specified in order to be accepted. (a pathlike "
            "object or string representing an existing file or a tuple of the "
            "form: (a float, a float, a float, a float))"
        )
        roi_mask_desc = (
            "Specify a masking region of interest. If defined,streamlines "
            "exiting the mask will be truncated. (a pathlike "
            "object or string representing an existing file or a tuple of the "
            "form: (a float, a float, a float, a float))"
        )
        act_image_desc = (
            " Use the Anatomically-Constrained Tractography "
            "framework during tracking (provided image must be "
            "in the 5TT ie five tissue type format(a pathlike object"
            "string representing an existing file)"
        )
        backtrack_desc = "Allow tracks to be truncated. (a boolean)"
        crop_at_gmwmi_desc = (
            "Crop streamline endpoints more precisely as they cross the "
            "GM-WM interface (a boolean)"
        )
        iFOD_power_desc = (
            "Raise the FOD to the power specified (default is 1/nsamples) "
            "(an integer)"
        )
        iFOD2_n_samples_desc = (
            "Set the number of FOD samples to take per"
            "step for the 2nd order (iFOD2) method. "
        )

        # Outputs description
        out_file_desc = (
            "Output file containing tracks (a pathlike object or "
            "string representing a file) "
        )
        output_seeds_desc = (
            "Out seed location of all successful streamlines ."
            "(a pathlike object or string representing a file) "
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
            "downsample_factor",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=downsample_factor_desc,
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
            Either(
                Undefined,
                Int(),
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
            "seed_dynamic",
            File(output=False, optional=True, desc=seed_dynamic_desc),
        )

        self.add_trait(
            "seed_gmwmi",
            File(output=False, optional=True, desc=seed_gmwmi_desc),
        )

        self.add_trait(
            "seed_grid_voxel",
            Tuple(
                File,
                Int,
                output=False,
                optional=True,
                desc=seed_grid_voxel_desc,
            ),
        )

        self.add_trait(
            "seed_image",
            File(output=False, optional=True, desc=seed_image_desc),
        )

        self.add_trait(
            "seed_rejection",
            File(output=False, optional=True, desc=seed_rejection_desc),
        )

        self.add_trait(
            "seed_rnd_voxel",
            Either(
                Undefined,
                Tuple(File, Int),
                default=Undefined,
                output=False,
                optional=True,
                desc=seed_rnd_voxel_desc,
            ),
        )

        self.add_trait(
            "seed_sphere",
            Either(
                Undefined,
                Tuple(Float, Float, Float, Float),
                default=Undefined,
                output=False,
                optional=True,
                desc=seed_sphere_desc,
            ),
        )

        self.add_trait(
            "tracto_seeds_number",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=tracto_seeds_number_desc,
            ),
        )

        self.add_trait(
            "tracto_max_attempts_per_seed_number",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=tracto_max_attempts_per_seed_number_desc,
            ),
        )

        self.add_trait(
            "tracto_seed_cutoff",
            Either(
                Undefined,
                Float(),
                output=False,
                optional=True,
                desc=tracto_seed_cutoff_desc,
            ),
        )

        self.add_trait(
            "tracto_seed_unidirectional",
            Bool(
                False,
                output=False,
                optional=True,
                desc=tracto_seed_unidirectional_desc,
            ),
        )

        self.add_trait(
            "tracto_seed_direction",
            Either(
                Undefined,
                Tuple(Float, Float, Float),
                default=Undefined,
                output=False,
                optional=True,
                desc=tracto_seed_direction_desc,
            ),
        )

        self.add_trait(
            "tracto_get_output_seeds",
            Bool(
                True,
                output=False,
                optional=True,
                desc=tracto_get_output_seeds_desc,
            ),
        )

        self.add_trait(
            "roi_excl",
            Either(
                File(),
                Tuple(Float, Float, Float, Float),
                default=Undefined,
                output=False,
                optional=True,
                desc=roi_excl_desc,
            ),
        )

        self.add_trait(
            "roi_incl",
            Either(
                File(),
                Tuple(Float, Float, Float, Float),
                output=False,
                optional=True,
                desc=roi_incl_desc,
            ),
        )

        self.add_trait(
            "roi_incl_ordered",
            Either(
                File(),
                Tuple(Float, Float, Float, Float),
                output=False,
                optional=True,
                desc=roi_incl_ordered_desc,
            ),
        )

        self.add_trait(
            "roi_mask",
            Either(
                File(),
                Tuple(Float, Float, Float, Float),
                output=False,
                optional=True,
                desc=roi_mask_desc,
            ),
        )

        self.add_trait(
            "act_image", File(output=False, optional=True, desc=act_image_desc)
        )

        self.add_trait(
            "backtrack",
            Bool(
                False,
                output=False,
                optional=True,
                desc=backtrack_desc,
            ),
        )

        self.add_trait(
            "crop_at_gmwmi",
            Bool(
                False,
                output=False,
                optional=True,
                desc=crop_at_gmwmi_desc,
            ),
        )

        self.add_trait(
            "iFOD_power",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=iFOD_power_desc,
            ),
        )

        self.add_trait(
            "iFOD2_n_samples",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=iFOD2_n_samples_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=False, desc=out_file_desc)
        )
        self.add_trait(
            "output_seeds",
            File(output=True, optional=True, desc=output_seeds_desc),
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.Tractography")

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
            valid_ext, in_ext, file_name = checkFileExt(self.in_file, EXT_DWI)
            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return self.make_initResult()

            if self.output_directory:
                # TODO: add number of streamline in name ?
                out_file_name = file_name + +"_tracto.tck"
                self.outputs["out_file"] = os.path.join(
                    self.output_directory, out_file_name
                )
                if self.tracto_get_output_seeds:
                    out_seed_file_name = (
                        file_name + "_tracto_out_seeds" + in_ext
                    )
                    self.outputs["output_seeds"] = os.path.join(
                        self.output_directory, out_seed_file_name
                    )

        if self.outputs:
            self.tags_inheritance(
                in_file=self.in_dti, out_file=self.outputs["out_file"]
            )
            if self.tracto_get_output_seeds:
                self.tags_inheritance(
                    in_file=self.in_dti, out_file=self.outputs["output_seeds"]
                )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Tractography, self).run_process_mia()
        # Nipype Tractoraphy process is not up to date for some option
        # so the "args" input is used for several options
        self.process.in_file = self.in_file
        self.process.out_file = self.out_file
        self.process.algorithm = self.algorithm
        self.process.cutoff = self.cutoff
        self.process.noprecompt = self.noprecompt
        self.process.select = self.select
        self.process.use_rk4 = self.use_rk4
        self.process.backtrack = self.backtrack
        self.process.crop_at_gmwmi = self.crop_at_gmwmi

        if self.angle:
            self.process.angle = self.angle
        if self.downsample_factor:
            self.process.downsample = self.downsample_factor
        if self.max_length:
            self.process.max_length = self.max_length
        if self.min_length:
            self.process.min_length = self.min_length
        if self.trials:
            self.process.trials = self.trials
        if self.step_size:
            self.process.step_size = self.step_size
        if self.seed_grid_voxel:
            self.process.seed_grid_voxel = self.seed_grid_voxel
        if self.seed_dynamic:
            self.process.seed_dynamic = self.seed_dynamic
        if self.seed_gmwmi:
            self.process.seed_gmwmi = self.seed_gmwmi
        if self.seed_image:
            self.process.seed_image = self.seed_image
        if self.seed_rejection:
            self.process.seed_rejection = self.seed_rejection
        if self.seed_rnd_voxel:
            self.process.seed_rnd_voxel = self.seed_rnd_voxel
        if self.seed_sphere:
            self.process.seed_sphere = self.seed_sphere
        args = ""
        if self.tracto_seeds_number:
            args += "-seeds" + self.tracto_seeds_number + ""
        if self.tracto_max_attempts_per_seed_number:
            args += (
                "-max_attempts_per_seed"
                + self.tracto_max_attempts_per_seed_number
                + ""
            )
        if self.tracto_seed_cutoff:
            args += "-seed_cutoff" + self.tracto_seed_cutoff + ""
        if self.tracto_seed_unidirectional:
            args += "-seed_unidirectional"
        if self.tracto_seed_direction:
            args += "-seed_direction" + self.tracto_seed_direction + ""
        if self.roi_excl:
            self.process.roi_excl = self.roi_excl
        if self.roi_incl:
            self.process.roi_incl = self.roi_incl
        if self.roi_incl_ordered:
            args += "" + self.roi_incl + " "
        if self.roi_mask:
            self.process.roi_mask = self.roi_mask
        if self.act_image:
            self.process.act_file = self.act_image
        if self.iFOD_power:
            self.process.power = self.iFOD_power
        if self.iFOD2_n_samples:
            self.process.in_samples = self.iFOD2_n_samples
        if self.tracto_get_output_seeds:
            self.process.output_seeds = self.output_seeds

        if args:
            self.process.args = args

        return self.process.run(configuration_dict={})


class TransformFSLConvert(ProcessMIA):
    """
    *Perform conversion between FSL’s transformation matrix format
    to mrtrix3’s. (transformconvert command)*

    Please, see the complete documentation for the `TransformFSLConvert brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/mrtrix/TransformFSLConvert.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(TransformFSLConvert, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "mrtrix"]

        # Mandatory inputs description
        in_file_desc = (
            "FLIRT input image (a pathlike object"
            "string representing an existing file)"
        )
        reference_desc = (
            "FLIRT reference image (a pathlike object"
            "string representing an existing file)"
        )
        in_transform_desc = (
            "FLIRT output transformation matrix (a pathlike object"
            "string representing an existing file)"
        )

        # Outputs description
        out_transform_desc = (
            "Output transformed affine in mrtrix’s format (a pathlike object "
            "or string representing a file) "
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "reference",
            File(output=False, optional=False, desc=reference_desc),
        )

        self.add_trait(
            "in_transform",
            File(output=False, optional=False, desc=in_transform_desc),
        )

        # Outputs traits
        self.add_trait(
            "out_transform",
            File(output=True, optional=False, desc=out_transform_desc),
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.mrtrix3.TransformFSLConvert")

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
        super(TransformFSLConvert, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_transform:
            valid_ext, in_ext, file_name = checkFileExt(
                self.in_file, {"MAT": "mat"}
            )
            if not valid_ext:
                print("\nThe transform matrice format is not recognized...!")
                return self.make_initResult()
            if self.output_directory:
                self.outputs["out_transform"] = os.path.join(
                    self.output_directory, file_name + "_mrtrix.txt"
                )

        if self.outputs:
            self.tags_inheritance(
                in_file=self.in_dti, out_file=self.outputs["out_transform"]
            )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(TransformFSLConvert, self).run_process_mia()
        self.process.flirt_import = True
        self.process.in_file = self.in_file
        self.process.in_transform = self.in_transform
        self.process.reference = self.reference
        self.process.out_transform = self.out_transform

        return self.process.run(configuration_dict={})
