# -*- coding: utf-8 -*-

"""
The mrtrix tools

:Contains:
    :Class:
        - DWIExtract
        - MRCat
        - MRConvert
        - MRMath

"""

import os

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################
from capsul.in_context import mrtrix
from nipype.interfaces.base import File, InputMultiPath

# populse_mia and mia_processes import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA
from traits.api import Bool, Either, Enum, Float, Int, List, String, Undefined

from mia_processes.utils import checkFileExt

EXT = {"NIFTI_GZ": "nii.gz", "NIFTI": "nii", "MIF": "mif"}


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
