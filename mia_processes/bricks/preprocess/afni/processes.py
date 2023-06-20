# -*- coding: utf-8 -*-

"""The afni preprocess library of the mia_processes package.

The purpose of this module is to customise the main afni preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - Automask
        - Calc
        - CalcDropTRs
        - Despike
        - FWHMx
        - GCOR
        - OutlierCount
        - QualityIndex
        - RefitDeoblique
        - SkullStrip
        - TShift
        - TStatMean
        - Volreg

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

import nibabel as nib

# nipype imports
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


class Automask(ProcessMIA):
    """
    *Create a brain-only mask of the image using AFNI 3dAutomask command*

    Please, see the complete documentation for the `Automask brick in
    the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/Automask.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Automask, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Mandatory inputs description
        in_file_desc = (
            "Input file (a pathlike object or string " "representing a file)."
        )
        # Optional inputs with default value description
        clfrac_desc = (
            "Sets the clip level fraction (must be 0.1-0.9)."
            "A small value will tend to make the mask larger "
            "(a float)"
        )
        out_brain_suffix_desc = "Suffix of the brain masked image (a string)"
        output_type_desc = (
            "Typecodes of the output image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the output image file "
            "(a string)."
        )
        # Optional inputs description
        dilate_desc = "Dilate the mask outwards (an integer)"
        erode_desc = "Erode the mask inwards (an integer)"

        # Outputs description
        out_file_desc = (
            "The brain mask file (a pathlike object or a "
            "string representing a file)."
        )

        brain_file_desc = (
            "The masked brain file (a pathlike object or a "
            "string representing a file)."
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits
        self.add_trait(
            "clfrac", Float(0.5, output=False, optional=True, desc=clfrac_desc)
        )

        self.add_trait(
            "out_brain_suffix",
            String(
                "_masked",
                output=False,
                optional=True,
                desc=out_brain_suffix_desc,
            ),
        )

        self.add_trait(
            "output_type",
            Enum(
                "NIFTI",
                "NIFTI_GZ",
                output=False,
                optional=True,
                desc=output_type_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String(
                "automask_", output=False, optional=True, desc=out_prefix_desc
            ),
        )

        # Optional inputs traits
        self.add_trait(
            "erode",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=erode_desc,
            ),
        )

        self.add_trait(
            "dilate",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=dilate_desc,
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.add_trait(
            "brain_file",
            File(output=True, optional=True, desc=brain_file_desc),
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.Automask")

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
        super(Automask, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.out_prefix == Undefined:
                self.out_prefix = "automask_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "automask" ...'
                )

            if self.out_brain_suffix == Undefined:
                self.out_brain_suffix = "_masked"
                print(
                    "The out_brain_suffix parameter is undefined."
                    'Automatically set to "_masked" ...'
                )

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is not recognized ...!")
                    return
                else:
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        self.out_prefix
                        + fileName
                        + "."
                        + EXT[self.output_type],
                    )

                    self.outputs["brain_file"] = os.path.join(
                        self.output_directory,
                        fileName
                        + self.out_brain_suffix
                        + "."
                        + EXT[self.output_type],
                    )

                    self.inheritance_dict[
                        self.outputs["out_file"]
                    ] = self.in_file

                    self.inheritance_dict[
                        self.outputs["brain_file"]
                    ] = self.in_file

            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Automask, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.outputtype = self.output_type
        self.process.out_file = self.out_file
        self.process.brain_file = self.brain_file
        self.process.clfrac = self.clfrac
        if self.erode:
            self.process.erode = self.erode
        if self.dilate:
            self.process.dilate = self.dilate

        return self.process.run(configuration_dict={})


class Calc(ProcessMIA):
    """
    *Voxel-by-voxel arithmetic on 3D datasets*

    Please, see the complete documentation for the `Calc brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/Calc.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Calc, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Inputs description
        in_file_a_desc = (
            "First input 3D file (a pathlike object or string "
            "representing a file)."
        )
        in_file_b_desc = (
            "Second input 3D file (a pathlike object or string "
            "representing a file)."
        )
        in_file_c_desc = (
            "Third input 3D file (a pathlike object or string "
            "representing a file)."
        )
        expr_desc = (
            "Arithmetic expression to apply between a, b and c " "(a string)."
        )
        output_type_desc = (
            "Typecodes of the output image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the output image file "
            "(a string)."
        )
        single_idx_desc = (
            "Volume index for in_file_a." "(an integer or Undefined)"
        )
        start_idx_desc = (
            "Start index for in_file_a (an integer"
            "or Undefined). Requires inputs: stop_idx"
        )
        stop_idx_desc = (
            "Stop index for in_file_a (an integer or Undefined)."
            "Requires inputs: start_idx."
        )

        # Outputs description
        out_file_desc = (
            "The calculated files (a pathlike object or a "
            "string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file_a",
            File(output=False, optional=False, desc=in_file_a_desc),
        )

        self.add_trait(
            "in_file_b",
            File(
                value=Undefined,
                output=False,
                optional=True,
                desc=in_file_b_desc,
            ),
        )

        self.add_trait(
            "in_file_c",
            File(
                value=Undefined,
                output=False,
                optional=True,
                desc=in_file_c_desc,
            ),
        )

        # Optional inputs with default value traits
        self.add_trait(
            "expr",
            String(Undefined, output=False, optional=True, desc=expr_desc),
        )
        self.expr = "a*step(b)"

        self.add_trait(
            "output_type",
            Enum(
                "NIFTI",
                "NIFTI_GZ",
                output=False,
                optional=True,
                desc=output_type_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String("c_", output=False, optional=True, desc=out_prefix_desc),
        )

        # Optional inputs traits
        self.add_trait(
            "single_idx",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=single_idx_desc,
            ),
        )

        self.add_trait(
            "start_idx",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=start_idx_desc,
            ),
        )
        self.add_trait(
            "stop_idx",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=stop_idx_desc,
            ),
        )
        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.Calc")

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
        super(Calc, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if (self.start_idx != Undefined and self.stop_idx == Undefined) or (
            self.start_idx == Undefined and self.stop_idx != Undefined
        ):
            print(
                '\nInitialisation failed. "start_idx" parameter required'
                '"stop_idx" parameters and vice versa'
            )
            return
        if self.in_file_a:
            if self.out_prefix == Undefined:
                self.out_prefix = "c_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "c" ...'
                )

            if self.output_directory:
                valid_ext, in_ext, fileName_a = checkFileExt(
                    self.in_file_a, EXT
                )

                if not valid_ext:
                    print("\nThe input image format is not recognized ...!")
                    return
                else:
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        self.out_prefix
                        + fileName_a
                        + "."
                        + EXT[self.output_type],
                    )

                    self.inheritance_dict[
                        self.outputs["out_file"]
                    ] = self.in_file_a

            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Calc, self).run_process_mia()
        self.process.in_file_a = self.in_file_a
        self.process.in_file_b = self.in_file_b
        self.process.in_file_c = self.in_file_c
        self.process.expr = self.expr
        self.process.outputtype = self.output_type
        self.process.out_file = self.out_file
        if self.single_idx:
            self.process.single_idx = self.single_idx
        if self.start_idx:
            self.process.start_idx = self.start_idx
        if self.stop_idx:
            self.process.stop_idx = self.stop_idx

        return self.process.run(configuration_dict={})


class CalcDropTRs(ProcessMIA):
    """
    *DropTRs of bold datasets (using AFNI 3dCalc command)*

    Please, see the complete documentation for the `CalcDropTRs brick in
    the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/CalcDropTRs.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(CalcDropTRs, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Inputs description
        in_file_desc = (
            "Input 3D file (a pathlike object or string "
            "representing a file)."
        )
        output_type_desc = (
            "Typecodes of the output image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the output image file "
            "(a string)."
        )
        start_idx_desc = "start index (inclusive) for in_file (an Int)."
        stop_idx_desc = "stop index (exclusive) for in_file (an Int)."
        # Outputs description
        out_file_desc = (
            "The TR cropped file (a pathlike object or a "
            "string representing a file)."
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits
        self.add_trait(
            "output_type",
            Enum(
                "NIFTI",
                "NIFTI_GZ",
                output=False,
                optional=True,
                desc=output_type_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String(
                "cropped_", output=False, optional=True, desc=out_prefix_desc
            ),
        )

        self.add_trait(
            "start_idx",
            Int(0, output=False, optional=True, desc=start_idx_desc),
        )

        self.add_trait(
            "stop_idx",
            Int(-1, output=False, optional=True, desc=stop_idx_desc),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.Calc")

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
        super(CalcDropTRs, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if not self.stop_idx or self.stop_idx == -1:
                print(
                    "\nWarning: stop_idx will be automatically set to"
                    "the length of input file"
                )
            elif self.stop_idx <= self.start_idx:
                print(
                    "\nError: stop_idx cannot be lower than or equal to"
                    "start_idx"
                )
                return

            if self.start_idx == 0 and self.stop_idx == -1:
                self.outputs["out_file"] = self.in_file
            else:
                if self.out_prefix == Undefined:
                    self.out_prefix = "cropped_"
                    print(
                        "The out_prefix parameter is undefined."
                        'Automatically set to "cropped" ...'
                    )

                if self.output_directory:
                    valid_ext, in_ext, fileName = checkFileExt(
                        self.in_file, EXT
                    )

                    if not valid_ext:
                        print(
                            "\nThe input image format is" " not recognized...!"
                        )
                        return
                    else:
                        self.outputs["out_file"] = os.path.join(
                            self.output_directory,
                            self.out_prefix
                            + fileName
                            + "."
                            + EXT[self.output_type],
                        )
                else:
                    print("No output_directory was found...!\n")
                    return

            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(CalcDropTRs, self).run_process_mia()

        if self.in_file == self.out_file:
            return

        try:
            imnii = nib.load(self.in_file)
            nb_volumes = imnii.dataobj.shape[3]

        except (
            nib.filebasedimages.ImageFileError,
            FileNotFoundError,
            TypeError,
        ) as e:
            print("\nError while opening input file" ": ", e)
            return

        self.process.in_file_a = self.in_file
        self.process.expr = "a"
        self.process.outputtype = self.output_type
        self.process.start_idx = self.start_idx
        if self.stop_idx == -1:
            self.process.stop_idx = nb_volumes - 1
        else:
            self.process.stop_idx = self.stop_idx
        self.process.out_file = self.out_file

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


class Despike(ProcessMIA):
    """
    *Removes ‘spikes’ from the 3D+time input dataset*

    Please, see the complete documentation for the `Despike brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/Despike.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Despike, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Inputs description
        in_file_desc = (
            "Input 3D file (a pathlike object or string "
            "representing a file)."
        )
        despike_desc = "Despike dataset only if true (boolean)."
        output_type_desc = (
            "Typecodes of the output image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the output image file "
            "(a string)."
        )

        # Outputs description
        out_file_desc = (
            "The despiked file (a pathlike object or a "
            "string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "despike",
            Bool(True, output=False, optional=True, desc=despike_desc),
        )

        self.add_trait(
            "output_type",
            Enum(
                "NIFTI",
                "NIFTI_GZ",
                output=False,
                optional=True,
                desc=output_type_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String("d_", output=False, optional=True, desc=out_prefix_desc),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.Despike")

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
        super(Despike, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.despike:
                if self.out_prefix == Undefined:
                    self.out_prefix = "d_"
                    print(
                        "The out_prefix parameter is undefined."
                        'Automatically set to "d" ...'
                    )

                if self.output_directory:
                    valid_ext, in_ext, fileName = checkFileExt(
                        self.in_file, EXT
                    )

                    if not valid_ext:
                        print(
                            "\nThe input image format is" " not recognized...!"
                        )
                        return
                    else:
                        self.outputs["out_file"] = os.path.join(
                            self.output_directory,
                            self.out_prefix
                            + fileName
                            + "."
                            + EXT[self.output_type],
                        )

                else:
                    print("No output_directory was found...!\n")
                    return
            else:
                self.outputs["out_file"] = self.in_file

            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Despike, self).run_process_mia()
        if not self.despike:
            return

        self.process.in_file = self.in_file
        self.process.outputtype = self.output_type
        self.process.out_file = self.out_file

        return self.process.run(configuration_dict={})


class FWHMx(ProcessMIA):
    """
    *Computes FWHMs for all sub-bricks in the input dataset, each one \
separately*

    Please, see the complete documentation for the `FWHMx brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/FWHMx.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(FWHMx, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Inputs description
        in_file_desc = (
            "Input image (a pathlike object or string representing a file)."
        )
        mask_file_desc = (
            "Mask image (a pathlike object or string representing a file)."
        )
        combine_desc = (
            "Combine the final measurements along each axis (a bool)."
        )
        detrend_desc = "Detrend to the specified order (a bool or an int)."
        out_prefix_desc = (
            "Specify the string to be prepended to the filenames of the "
            "output image file (a string)."
        )
        args_desc = ""

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
            "mask_file",
            File(Undefined, output=False, optional=True, desc=mask_file_desc),
        )

        self.add_trait(
            "combine",
            Bool(True, optional=True, output=False, desc=combine_desc),
        )

        self.add_trait(
            "detrend",
            Either(
                Bool(),
                Int(),
                default=True,
                output=False,
                optional=True,
                desc=detrend_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String("fwhm_", output=False, optional=True, desc=out_prefix_desc),
        )

        self.add_trait(
            "args", String("", output=False, optional=True, desc=args_desc)
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.FWHMx")

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
                self.out_prefix = "fwhm_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "fwhm" ...'
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
        super(FWHMx, self).run_process_mia()

        self.process.in_file = self.in_file
        if self.mask_file:
            self.process.mask = self.mask_file
        self.process.combine = self.combine
        self.process.detrend = self.detrend
        self.process.out_file = self.out_file
        if self.args:
            self.process.args = self.args

        return self.process.run(configuration_dict={})


class GCOR(ProcessMIA):
    """
    *Computes the average correlation between every voxel and every other \
voxel, over any given mask*

    Please, see the complete documentation for the `GCOR brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/GCOR.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(GCOR, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Inputs description
        in_file_desc = (
            "Input image (a pathlike object or string " "representing a file)."
        )
        mask_file_desc = (
            "Mask image (a pathlike object or string " "representing a file)."
        )
        nfirst_desc = "Specify number of initial TRs to ignore"
        no_demean_desc = "Do not (need to) demean as first step."

        # Outputs description
        out_desc = "Global correlation value (a float)."

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "mask_file",
            File(Undefined, output=False, optional=True, desc=mask_file_desc),
        )

        self.add_trait(
            "nfirst",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=nfirst_desc,
            ),
        )

        self.add_trait(
            "no_demean",
            Bool(False, output=False, optional=True, desc=no_demean_desc),
        )

        # Outputs traits
        self.add_trait(
            "out",
            Float(output=True, nipype_process_name="_out", desc=out_desc),
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.GCOR")

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
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is not recognized ...!")
                    return

            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(GCOR, self).run_process_mia()

        self.process.in_file = self.in_file
        self.process.mask = self.mask_file
        self.process.nfirst = self.nfirst
        self.process.no_demean = self.no_demean

        return self.process.run(configuration_dict={})


class OutlierCount(ProcessMIA):
    """
    *Computes outliers for all sub-bricks in the input dataset, each one \
separately*

    Please, see the complete documentation for the `OutlierCount brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/OutlierCount.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(OutlierCount, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Inputs description
        in_file_desc = (
            "Input image (a pathlike object or string " "representing a file)."
        )
        automask_desc = (
            "Clip off small voxels (a boolean). "
            "Mutually exclusive with mask"
        )
        autoclip_desc = (
            "Clip off small voxels (a boolean). "
            "Mutually exclusive with mask"
        )
        interval_desc = (
            "Write out the median + 3.5 MAD of outlier "
            "count with each timepoint. "
        )
        mask_file_desc = (
            "Mask image. Compute correlation only across "
            "masked voxels.(a pathlike object or string "
            "representing a file). Mutually exclusive "
            "with automask and autoclip"
        )
        fraction_desc = (
            "Combine the final measurements " "along each axis (a bool)."
        )
        interval_desc = (
            "Write out the median + 3.5 MAD of outlier "
            "count with each timepoint."
        )
        legendre_desc = "Use Legendre polynomials."
        polort_desc = "Detrend each voxel timeseries with polynomials. "
        qthr_desc = "Indicate a value for q to compute alpha."
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the output image file "
            "(a string)."
        )

        # Outputs description
        out_file_desc = (
            "The output text file (a pathlike object or a "
            "string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "mask_file",
            File(Undefined, output=False, optional=True, desc=mask_file_desc),
        )

        self.add_trait(
            "autoclip",
            Bool(False, output=False, optional=True, desc=autoclip_desc),
        )

        self.add_trait(
            "automask",
            Bool(False, output=False, optional=True, desc=automask_desc),
        )

        self.add_trait(
            "fraction",
            Bool(True, optional=True, output=False, desc=fraction_desc),
        )

        self.add_trait(
            "interval",
            Bool(False, optional=True, output=False, desc=interval_desc),
        )

        self.add_trait(
            "legendre",
            Bool(False, optional=True, output=False, desc=legendre_desc),
        )

        self.add_trait(
            "out_prefix",
            String(
                "outliers_", output=False, optional=True, desc=out_prefix_desc
            ),
        )

        self.add_trait(
            "polort",
            Either(
                Int(),
                Undefined,
                default=Undefined,
                optional=True,
                output=False,
                desc=polort_desc,
            ),
        )

        self.add_trait(
            "qthr",
            Int(
                0.001,
                min=0.0,
                max=1.0,
                optional=True,
                output=False,
                desc=qthr_desc,
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.OutlierCount")

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
            if (self.automask and self.mask_file) or (
                self.autoclip and self.mask_file
            ):
                print(
                    "Initisalisation failed.."
                    "Automask and mask_file parameters are "
                    "mutually esclusive"
                )
                return

            if self.out_prefix == Undefined:
                self.out_prefix = "outliers_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "outliers" ...'
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
        super(OutlierCount, self).run_process_mia()

        self.process.in_file = self.in_file
        if self.automask:
            self.process.automask = self.automask
        if self.autoclip:
            self.process.autoclip = self.autoclip
        if self.mask_file:
            self.process.autoclip = Undefined
            self.process.automask = Undefined
            self.process.mask = self.mask_file
        self.process.fraction = self.fraction
        self.process.interval = self.interval
        self.process.legendre = self.legendre
        self.process.polort = self.polort
        self.process.qthr = self.qthr
        self.process.out_file = self.out_file

        return self.process.run(configuration_dict={})


class QualityIndex(ProcessMIA):
    """
    *Computes a quality index for each sub-brick in a 3D+time dataset*

    Please, see the complete documentation for the `QualityIndex brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/QualityIndex.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(QualityIndex, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Inputs description
        in_file_desc = (
            "A bold file (a pathlike object or string " "representing a file)."
        )
        automask_desc = (
            "Clip off small voxels (a boolean). "
            "Mutually exclusive with mask"
        )
        autoclip_desc = (
            "Clip off small voxels (a boolean). "
            "Mutually exclusive with mask"
        )
        interval_desc = (
            "Write out the median + 3.5 MAD of outlier "
            "count with each timepoint. "
        )
        mask_file_desc = (
            "Mask image. Compute correlation only across "
            "masked voxels.(a pathlike object or string "
            "representing a file). Mutually exclusive "
            "with automask and autoclip"
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the output image file "
            "(a string)."
        )
        quadrant_desc = (
            "Similar to -spearman, but using 1 minus "
            "the quadrant correlation coefficient as "
            "the quality index"
        )
        spearman_desc = (
            "Quality index is 1 minus the Spearman (rank) "
            "correlation coefficient of each sub-brick with "
            "the median sub-brick. "
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
            "automask",
            Bool(True, output=False, optional=True, desc=automask_desc),
        )

        self.add_trait(
            "autoclip",
            Bool(False, output=False, optional=True, desc=autoclip_desc),
        )

        self.add_trait(
            "interval",
            Bool(False, optional=True, output=False, desc=interval_desc),
        )

        self.add_trait(
            "mask_file",
            File(Undefined, output=False, optional=True, desc=mask_file_desc),
        )

        self.add_trait(
            "out_prefix",
            String("QI_", output=False, optional=True, desc=out_prefix_desc),
        )

        self.add_trait(
            "quadrant",
            Bool(False, optional=True, output=False, desc=quadrant_desc),
        )

        self.add_trait(
            "spearman",
            Bool(False, optional=True, output=False, desc=spearman_desc),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.QualityIndex")

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

        if (self.automask and self.mask_file) or (
            self.autoclip and self.mask_file
        ):
            print(
                "Initisalisation failed.."
                "Automask and mask_file parameters are mutually esclusive"
            )
            return

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.out_prefix == Undefined:
                self.out_prefix = "QI_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "QI" ...'
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
        super(QualityIndex, self).run_process_mia()

        self.process.in_file = self.in_file
        if self.autoclip:
            self.process.autoclip = self.autoclip
        if self.automask:
            self.process.automask = self.automask
        if self.mask_file:
            self.process.autoclip = Undefined
            self.process.automask = Undefined
            self.process.mask = self.mask_file
        self.process.interval = self.interval
        self.process.quadrant = self.quadrant
        self.process.spearman = self.spearman
        self.process.out_file = self.out_file

        return self.process.run(configuration_dict={})


class RefitDeoblique(ProcessMIA):
    """
    *Deoblique dataset*

    Please, see the complete documentation for the `RefitDeoblique brick in
    the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/RefitDeoblique.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(RefitDeoblique, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Inputs description
        in_file_desc = (
            "Input file (a pathlike object or string " "representing a file)."
        )
        deoblique_desc = "Deoblique dataset only if true (boolean)."

        # Outputs description
        out_file_desc = (
            "The deobliqued file (a pathlike object or a "
            "string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "deoblique",
            Bool(True, output=False, optional=True, desc=deoblique_desc),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.Refit")

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
        super(RefitDeoblique, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            self.outputs["out_file"] = self.in_file

            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(RefitDeoblique, self).run_process_mia()

        self.process.in_file = self.in_file
        self.process.deoblique = self.deoblique

        return self.process.run(configuration_dict={})


class SkullStrip(ProcessMIA):
    """
    *From MRI T1-weighted images, extract the brain from surrounding tissue*

    Please, see the complete documentation for the `SkullStrip brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/SkullStrip.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(SkullStrip, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Inputs description
        in_file_desc = (
            "A 3D-T1 file to skull-strip (a pathlike object"
            "or string representing a file)."
        )
        output_type_desc = (
            "Typecodes of the output image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the skull-stripped image file(s) "
            "(a string)."
        )

        # Outputs description
        out_file_desc = (
            "The skull-stripped files (a pathlike object or a "
            "string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "output_type",
            Enum(
                "NIFTI",
                "NIFTI_GZ",
                output=False,
                optional=True,
                desc=output_type_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String("ss_", output=False, optional=True, desc=out_prefix_desc),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.SkullStrip")

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
        super(SkullStrip, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.out_prefix == Undefined:
                self.out_prefix = "ss_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "ss" ...'
                )

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is" " not recognized...!")
                    return
                else:
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        self.out_prefix
                        + fileName
                        + "."
                        + EXT[self.output_type],
                    )

                    self.inheritance_dict[
                        self.outputs["out_file"]
                    ] = self.in_file

            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(SkullStrip, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.outputtype = self.output_type
        self.process.out_file = self.out_file

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


class TShift(ProcessMIA):
    """
    *Slice-time correction of bold images*

    Please, see the complete documentation for the `TShift brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/TShift.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(TShift, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Mandatory inputs description
        in_file_desc = (
            "A bold file to be time-shifted (a pathlike object"
            "or string representing a file)."
        )
        # Optional inputs with default value description
        interpolation_desc = (
            "Interpolation methods (one of ‘Fourier’ or"
            "‘linear’ or ‘cubic’ or ‘quintic’ or ‘heptic’)"
        )
        output_type_desc = (
            "Typecodes of the output image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "output filename (a string)."
        )
        rlt_desc = (
            "Before shifting, remove the mean and linear trend" "(a boolean)"
        )
        rltplus_desc = (
            "Before shifting, remove the mean and linear trend "
            "and later put back the mean. (a boolean)"
        )
        slice_encoding_dir_desc = (
            "Direction in which slice_timing is"
            "specified (default: k)."
            "If negative,slice_timing is defined"
            "in reverse order, that is, the first"
            "entry corresponds to the slice with "
            "the largest index, and the final entry"
            "corresponds to slice index zero."
        )
        # Optional inputs description
        ignore_desc = (
            "Ignore the first set of points specified." "(an integer)"
        )
        slice_timing_desc = (
            "Time offsets from the volume acquisition onset"
            "for each slice. (a string representing an"
            "existing file or a list of floats)."
            'Mutually exclusive with "tpattern" parameters'
        )
        tpattern_desc = (
            "Use specified slice time pattern rather than one in"
            "header . One of (‘alt+z’ or ‘altplus’ or ‘alt+z2’"
            "or ‘alt-z’ or ‘altminus’ or ‘alt-z2’ or ‘seq+z’"
            "or ‘seqplus’ or‘seq-z’ or ‘seqminus’"
            'Mutually exclusive with "slice_timming" parameters'
        )
        tr_desc = (
            "Manually set the TR. You can add suffix “s” for seconds or"
            "“ms” for milliseconds (a string)"
        )
        tslice_desc = (
            "Align each slice to time offset of given slice."
            "(an integer). Mutually exclusive with tzero parameter."
        )
        tzero_desc = (
            " Align each slice to given time offset"
            "(a float). Mutually exclusive with tslice parameter."
        )
        # Outputs description
        out_file_desc = (
            "The time shifted file (a pathlike object or a "
            "string representing a file)."
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits
        self.add_trait(
            "interpolation",
            Enum(
                "Fourier",
                "linear",
                "cubic",
                "quintic",
                "heptic",
                output=False,
                optional=True,
                desc=interpolation_desc,
            ),
        )

        self.add_trait(
            "output_type",
            Enum(
                "NIFTI",
                "NIFTI_GZ",
                output=False,
                optional=True,
                desc=output_type_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String(
                "st_corr_", output=False, optional=True, desc=out_prefix_desc
            ),
        )

        self.add_trait(
            "rlt", Bool(False, output=False, optional=True, desc=rlt_desc)
        )

        self.add_trait(
            "rltplus",
            Bool(False, output=False, optional=True, desc=rltplus_desc),
        )

        self.add_trait(
            "slice_encoding_dir",
            Enum(
                "k",
                "k-",
                default="k",
                output=False,
                optional=True,
                desc=slice_encoding_dir_desc,
            ),
        )
        # Optional inputs traits
        self.add_trait(
            "ignore",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=ignore_desc,
            ),
        )

        self.add_trait(
            "slice_timing",
            Either(
                List(Float()),
                File(),
                output=False,
                optional=True,
                desc=slice_timing_desc,
            ),
        )

        self.add_trait(
            "tpattern",
            Either(
                Undefined,
                Enum(
                    "alt+z",
                    "altplus",
                    "alt+z2",
                    "alt-z",
                    "altminus",
                    "alt-z2",
                    "seq+z",
                    "seqplus",
                    "seq-z",
                    "seqminus",
                ),
                default=Undefined,
                output=False,
                optional=True,
                desc=tpattern_desc,
            ),
        )

        self.add_trait(
            "tr",
            Either(
                Undefined,
                String(),
                default=Undefined,
                output=False,
                optional=True,
                desc=tr_desc,
            ),
        )

        self.add_trait(
            "tslice",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=tslice_desc,
            ),
        )

        self.add_trait(
            "tzero",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=tzero_desc,
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.TShift")

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
        super(TShift, self).list_outputs()

        if self.slice_timing and self.tpattern != Undefined:
            print(
                "\nInitialisation failed. "
                'Please, define only "slice_timing" paremeters or '
                '"tpattern" parameters (set the other as Undefined) ...!'
            )
            return

        if self.tslice != Undefined and self.tzero != Undefined:
            print(
                "\nInitialisation failed. "
                'Please, define only "tslice" paremeters or '
                '"tzero" parameters (set the other as Undefined) ...!'
            )
            return

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.slice_timing or self.tpattern != Undefined:
                if self.out_prefix == Undefined:
                    self.out_prefix = "st_corr_"
                    print(
                        "The out_prefix parameter is undefined. "
                        'Automatically set to "st_corr" ...'
                    )

                if self.output_directory:
                    valid_ext, in_ext, fileName = checkFileExt(
                        self.in_file, EXT
                    )

                    if not valid_ext:
                        print(
                            "\nThe input image format is" " not recognized...!"
                        )
                        return
                    else:
                        self.outputs["out_file"] = os.path.join(
                            self.output_directory,
                            self.out_prefix
                            + fileName
                            + "."
                            + EXT[self.output_type],
                        )

                else:
                    print("No output_directory was found...!\n")
                    return
            else:
                self.outputs["out_file"] = self.in_file

            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(TShift, self).run_process_mia()
        if not self.slice_timing and self.tpattern == Undefined:
            return

        self.process.in_file = self.in_file
        self.process.outputtype = self.output_type
        self.process.out_file = self.out_file
        if self.tpattern:
            self.process.tpattern = self.tpattern
        if self.slice_timing:
            self.process.slice_timing = self.slice_timing
        self.process.slice_encoding_direction = self.slice_encoding_dir
        self.process.rlt = self.rlt
        self.process.rltplus = self.rltplus
        self.process.interp = self.interpolation
        if self.out_prefix:
            self.process.out_prefix = self.out_prefix
        if self.tslice:
            self.process.tslice = self.tslice
        if self.tzero:
            self.process.tzero = self.tzero
        if self.ignore:
            self.process.ignore = self.ignore
        if self.tr:
            self.process.tr = self.tr

        return self.process.run(configuration_dict={})


class TStatMean(ProcessMIA):
    """
    *Mean of bold images (using AFNI 3dTstat)*

    Please, see the complete documentation for the `TStatMean brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/TStatMean.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(TStatMean, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Inputs description
        in_file_desc = (
            "A bold file to be averaged (a pathlike object"
            "or string  representing a file)."
        )
        output_type_desc = (
            "Typecodes of the output image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the skull-stripped image file(s) "
            "(a string)."
        )

        # Outputs description
        out_file_desc = (
            "The time shifted file (a pathlike object or a "
            "string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "output_type",
            Enum(
                "NIFTI",
                "NIFTI_GZ",
                output=False,
                optional=True,
                desc=output_type_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String("mean_", output=False, optional=True, desc=out_prefix_desc),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.TStat")

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
        super(TStatMean, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.out_prefix == Undefined:
                self.out_prefix = "mean_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "mean" ...'
                )

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is" " not recognized...!")
                    return
                else:
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        self.out_prefix
                        + fileName
                        + "."
                        + EXT[self.output_type],
                    )

                    self.inheritance_dict[
                        self.outputs["out_file"]
                    ] = self.in_file

            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(TStatMean, self).run_process_mia()

        self.process.in_file = self.in_file
        self.process.args = "-mean"
        self.process.outputtype = self.output_type
        self.process.out_file = self.out_file

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


class Volreg(ProcessMIA):
    """
    *Register an input volume to a base volume using AFNI 3dvolreg*

    Please, see the complete documentation for the `Volreg brick in
    the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/Volreg.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Volreg, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["afni", "nipype"]

        # Mandatory inputs description
        in_file_desc = (
            "Input (a pathlike object or string representing a file)."
        )
        # Optional inputs with default value description
        copyorigin_desc = "Copy base file origin coords to output (a boolean)"
        interpolation_desc = (
            "Spatial interpolation methods (Either "
            "‘Fourier’ or ‘cubic’ or ‘heptic’ "
            "or ‘quintic’ or ‘linear’. Default is heptic)"
        )
        output_type_desc = (
            "Typecodes of the output image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the registered image file(s) "
            "(a string)."
        )
        save_md1d_file_desc = "Save md1d file (a boolean)"
        save_oned_matrix_desc = "Save oned matrix (a boolean)"
        timeshift_desc = "Time shift to mean slice time offset."
        twopass_desc = (
            "Do two passes of the registration algorithm:"
            "(1) with smoothed base and data bricks, with linear"
            " interpolation, to get a crude alignment, then"
            "(2) with the input base and data bricks, to get a "
            "fine alignment."
            "This method is useful when aligning high-resolution"
            "datasets that may need to be moved more than"
            " a few voxels to be aligned. (a boolean)."
        )
        zpad_desc = (
            "Zeropad around the edges by ‘n’ voxels during rotations "
            "(an integer)."
        )
        # Optional inputs description
        in_weight_volume_desc = (
            "Weights for each voxel specified by a file"
            "with an optional volume number "
            "(defaults to 0). A a tuple of the form:"
            "(a pathlike object or string representing"
            " an existing file, an integer)"
        )
        # Outputs description
        md1d_file_desc = (
            "Max displacement outputfile (a pathlike object"
            "or a  string representing a file)."
        )
        oned_file_desc = (
            "The movement parameters file (a pathlike object"
            "or a  string representing a file)."
        )
        oned_matrix_desc = (
            "The transformation matrix (a pathlike object"
            "or a  string representing a file)."
        )
        out_file_desc = (
            "The registered file (a pathlike object or a "
            "string representing a file)."
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits
        self.add_trait(
            "copyorigin",
            Bool(False, output=False, optional=True, desc=copyorigin_desc),
        )

        self.add_trait(
            "interpolation",
            Enum(
                "heptic",
                "Fourier",
                "cubic",
                "quintic",
                "linear",
                default="heptic",
                output=False,
                optional=True,
                desc=interpolation_desc,
            ),
        )

        self.add_trait(
            "output_type",
            Enum(
                "NIFTI",
                "NIFTI_GZ",
                output=False,
                optional=True,
                desc=output_type_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String("reg_", output=False, optional=True, desc=out_prefix_desc),
        )

        self.add_trait(
            "save_oned_matrix",
            Bool(
                False, output=False, optional=True, desc=save_oned_matrix_desc
            ),
        )

        self.add_trait(
            "save_md1d_file",
            Bool(False, output=False, optional=True, desc=save_md1d_file_desc),
        )

        self.add_trait(
            "timeshift",
            Bool(False, output=False, optional=True, desc=timeshift_desc),
        )

        self.add_trait(
            "twopass",
            Bool(False, output=False, optional=True, desc=twopass_desc),
        )

        self.add_trait(
            "zpad", Int(4, output=False, optional=True, desc=zpad_desc)
        )

        # Optional inputs traits
        self.add_trait(
            "in_weight_volume",
            Either(
                Undefined,
                Tuple(String, Int(0)),
                default=Undefined,
                output=False,
                optional=True,
                desc=in_weight_volume_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "md1d_file", File(output=True, optional=True, desc=md1d_file_desc)
        )

        self.add_trait("oned_file", File(output=True, desc=oned_file_desc))

        self.add_trait(
            "oned_matrix",
            File(output=True, optional=True, desc=oned_matrix_desc),
        )

        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.afni.Volreg")

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
        super(Volreg, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.out_prefix == Undefined:
                self.out_prefix = "reg_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "reg" ...'
                )

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is" " not recognized...!")
                    return
                else:
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        self.out_prefix
                        + fileName
                        + "."
                        + EXT[self.output_type],
                    )

                    self.outputs["oned_file"] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName + "_oned.txt",
                    )

                    self.inheritance_dict[
                        self.outputs["out_file"]
                    ] = self.in_file

                    self.inheritance_dict[
                        self.outputs["oned_file"]
                    ] = self.in_file

                    if self.save_oned_matrix:
                        self.outputs["oned_matrix"] = os.path.join(
                            self.output_directory,
                            self.out_prefix
                            + fileName
                            + "_oned_matrix.aff12.1D",
                        )

                        self.inheritance_dict[
                            self.outputs["oned_matrix"]
                        ] = self.in_file

                    if self.save_md1d_file:
                        self.outputs["md1d_file"] = os.path.join(
                            self.output_directory,
                            self.out_prefix + fileName + "_md.1D",
                        )

                        self.inheritance_dict[
                            self.outputs["md1d_file"]
                        ] = self.in_file

            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Volreg, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.outputtype = self.output_type
        self.process.interp = self.interpolation
        if self.twopass:
            self.process.args = "-twopass"
        if self.zpad:
            self.process.zpad = self.zpad
        if self.copyorigin:
            self.process.copyorigin = self.copyorigin
        if self.in_weight_volume:
            self.process.in_weight_volume = self.in_weight_volume
        self.process.out_file = self.out_file
        self.process.oned_file = self.oned_file
        if self.save_oned_matrix:
            self.process.oned_matrix_save = self.oned_matrix
        if self.save_md1d_file:
            self.process.md1d_file = self.md1d_file
        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})
