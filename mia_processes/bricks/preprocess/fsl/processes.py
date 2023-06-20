# -*- coding: utf-8 -*-

"""The fsl preprocess library of the mia_processes package.

The purpose of this module is to customise the main fsl preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - BetSurfacesExtraction
        - FastSegment
        - Smooth

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
from traits.api import Bool, Either, Enum, Float, Int, List, String, Undefined

from mia_processes.utils import checkFileExt

EXT = {"NIFTI_GZ": "nii.gz", "NIFTI": "nii"}


class BetSurfacesExtraction(ProcessMIA):
    """
    *Surfaces (skull, inskull, outskull, outskin) extraction using BET (FSL)*

    Please, see the complete documentation for the `SurfacesExtraction brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/fsl/SurfacesExtraction.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(BetSurfacesExtraction, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["fsl", "nipype"]

        # Mandatory inputs description
        in_file_desc = (
            "Input file to skull strip (a pathlike object"
            "string representing an existing file)"
        )
        # Optional inputs with default value description
        output_type_desc = (
            "Typecodes of the output NIfTI image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        # Outputs description
        inskull_mask_file_desc = (
            "Inskull mask file (a pathlike object "
            "or string representing a file)"
        )
        inskull_mesh_file_desc = (
            "Inskull mesh file (a pathlike object "
            "or string representing a file)"
        )
        out_file_desc = (
            "Skull stripped image (a pathlike object "
            "or string representing a file)"
        )
        outskin_mask_file_desc = (
            "Outskin mask file (a pathlike object "
            "or string representing a file)"
        )
        outskin_mesh_file_desc = (
            "Outskin mesh file (a pathlike object "
            "or string representing a file)"
        )
        outskull_mask_file_desc = (
            "Outskull mask file (a pathlike object "
            "or string representing a file)"
        )
        outskull_mesh_file_desc = (
            "Outskull mesh file (a pathlike object "
            "or string representing a file)"
        )
        skull_mask_file_desc = (
            "Skull mask file (a pathlike object "
            "or string representing a file)"
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

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=True, desc=out_file_desc)
        )

        self.add_trait(
            "inskull_mask_file",
            File(output=True, optional=True, desc=inskull_mask_file_desc),
        )

        self.add_trait(
            "inskull_mesh_file",
            File(output=True, optional=True, desc=inskull_mesh_file_desc),
        )

        self.add_trait(
            "outskin_mask_file",
            File(output=True, optional=True, desc=outskin_mask_file_desc),
        )

        self.add_trait(
            "outskin_mesh_file",
            File(output=True, optional=True, desc=outskin_mesh_file_desc),
        )

        self.add_trait(
            "outskull_mask_file",
            File(output=True, optional=True, desc=outskull_mask_file_desc),
        )

        self.add_trait(
            "outskull_mesh_file",
            File(output=True, optional=True, desc=outskull_mesh_file_desc),
        )

        self.add_trait(
            "skull_mask_file",
            File(output=True, optional=True, desc=skull_mask_file_desc),
        )

        self.init_default_traits()

        # To suppress the "FSLOUTPUTTYPE environment
        # variable is not set" nipype warning:
        if "FSLOUTPUTTYPE" not in os.environ:
            os.environ["FSLOUTPUTTYPE"] = self.output_type

        self.init_process("nipype.interfaces.fsl.BET")

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
        super(BetSurfacesExtraction, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is" " not recognized...!")
                return

            self.process.output_type = self.output_type
            self.process.in_file = self.in_file
            self.process.surfaces = True

            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory,
                    os.path.split(self.process._out_file)[1],
                )
                self.outputs["outskin_mask_file"] = os.path.join(
                    self.output_directory,
                    os.path.split(self.process._outskin_mask_file)[1],
                )
                self.outputs["outskin_mesh_file"] = os.path.join(
                    self.output_directory,
                    os.path.split(self.process._outskin_mesh_file)[1],
                )
                self.outputs["outskull_mask_file"] = os.path.join(
                    self.output_directory,
                    os.path.split(self.process._outskull_mask_file)[1],
                )
                self.outputs["outskull_mesh_file"] = os.path.join(
                    self.output_directory,
                    os.path.split(self.process._outskull_mesh_file)[1],
                )
                self.outputs["inskull_mask_file"] = os.path.join(
                    self.output_directory,
                    os.path.split(self.process._inskull_mask_file)[1],
                )
                self.outputs["inskull_mesh_file"] = os.path.join(
                    self.output_directory,
                    os.path.split(self.process._inskull_mesh_file)[1],
                )
                self.outputs["skull_mask_file"] = os.path.join(
                    self.output_directory,
                    os.path.split(self.process._skull_mask_file)[1],
                )
        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file
            self.inheritance_dict[
                self.outputs["outskin_mask_file"]
            ] = self.in_file
            self.inheritance_dict[
                self.outputs["outskin_mesh_file"]
            ] = self.in_file
            self.inheritance_dict[
                self.outputs["outskull_mask_file"]
            ] = self.in_file
            self.inheritance_dict[
                self.outputs["outskull_mesh_file"]
            ] = self.in_file
            self.inheritance_dict[
                self.outputs["inskull_mask_file"]
            ] = self.in_file
            self.inheritance_dict[
                self.outputs["inskull_mesh_file"]
            ] = self.in_file
            self.inheritance_dict[
                self.outputs["skull_mask_file"]
            ] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(BetSurfacesExtraction, self).run_process_mia()
        self.process.output_type = self.output_type
        self.process.in_file = self.in_file

        # default inputs
        self.process.surfaces = True

        return self.process.run(configuration_dict={})


class FastSegment(ProcessMIA):
    """
    *Brain tissue segmentation using fsl.FAST*

    Please, see the complete documentation for the `FastSegment brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/htmldocumentation/bricks/preprocess/fsl/FastSegment.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(FastSegment, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["fsl", "nipype"]

        # Mandatory inputs description
        in_file_desc = "File to Segment"
        # Optional inputs with default value description
        img_type_desc = (
            "Int specifying type of image:" "(1 = T1, 2 = T2, 3 = PD). "
        )
        output_type_desc = (
            "Typecodes of the output NIfTI image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        segments_desc = (
            "Outputs a separate binary image for each tissue type "
            "(a boolean)."
        )
        # Outputs description
        mixeltype_desc = (
            "Path/name of mixeltype volume file _mixeltype"
            "(a pathlike object or string representing"
            "a file)"
        )
        partial_volume_map_desc = (
            "Partial volume map (_pveseg)"
            "(a pathlike object or string representing"
            "a file)"
        )
        partial_volume_files_desc = (
            "Partial volume files (a list of items " "which are file names."
        )
        tissue_class_files_desc = (
            "Binary segmented volume files, one image"
            "per class (a list of items which are"
            "file names.)"
        )
        tissue_class_map_desc = (
            "Binary segmented volume file "
            "(a pathlike object or string representing"
            "a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits
        self.add_trait(
            "img_type", Int(1, output=False, optional=True, desc=img_type_desc)
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
            "segments",
            Bool(True, output=False, optional=True, desc=segments_desc),
        )

        # Outputs traits
        self.add_trait(
            "mixeltype", File(output=True, optional=True, desc=mixeltype_desc)
        )

        self.add_trait(
            "partial_volume_map",
            File(output=True, optional=True, desc=partial_volume_map_desc),
        )

        self.add_trait(
            "partial_volume_files",
            List(File(), output=True, desc=partial_volume_files_desc),
        )

        self.add_trait(
            "tissue_class_files",
            List(
                File(),
                output=True,
                optional=True,
                desc=tissue_class_files_desc,
            ),
        )

        self.add_trait(
            "tissue_class_map", File(output=True, desc=tissue_class_map_desc)
        )

        self.init_default_traits()

        # To suppress the "FSLOUTPUTTYPE environment
        # variable is not set" nipype warning:
        if "FSLOUTPUTTYPE" not in os.environ:
            os.environ["FSLOUTPUTTYPE"] = self.output_type

        self.init_process("nipype.interfaces.fsl.FAST")

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
        super(FastSegment, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is" " not recognized...!")
                return

            self.process.output_type = self.output_type
            self.process.segments = self.segments

            if self.output_directory:
                _, fileIval = os.path.split(self.in_file)
                self.process.out_basename = os.path.join(
                    self.output_directory, fileIval
                )

                self.outputs["tissue_class_map"] = os.path.join(
                    self.output_directory,
                    os.path.split(self.process._tissue_class_map)[1],
                )

                self.outputs["partial_volume_map"] = os.path.join(
                    self.output_directory,
                    os.path.split(self.process._partial_volume_map)[1],
                )

                self.outputs["mixeltype"] = os.path.join(
                    self.output_directory,
                    os.path.split(self.process._mixeltype)[1],
                )

                self.outputs["partial_volume_files"] = []
                for out_val in self.process._partial_volume_files:
                    self.outputs["partial_volume_files"].append(
                        os.path.join(
                            self.output_directory, os.path.split(out_val)[1]
                        )
                    )

                self.outputs["tissue_class_files"] = []
                if self.segments:
                    for out_val in self.process._tissue_class_files:
                        self.outputs["tissue_class_files"].append(
                            os.path.join(
                                self.output_directory,
                                os.path.split(out_val)[1],
                            )
                        )
            else:
                print("No output_directory was found...!\n")
                return

        if self.outputs:
            self.inheritance_dict[
                self.outputs["tissue_class_map"]
            ] = self.in_file
            self.inheritance_dict[
                self.outputs["partial_volume_map"]
            ] = self.in_file
            self.inheritance_dict[self.outputs["mixeltype"]] = self.in_file
            if self.outputs["partial_volume_files"][0]:
                for out_val in self.outputs["partial_volume_files"]:
                    self.inheritance_dict[out_val] = self.in_file
            if self.outputs["tissue_class_files"][0]:
                for out_val in self.outputs["tissue_class_files"]:
                    self.inheritance_dict[out_val] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(FastSegment, self).run_process_mia()
        self.process.img_type = self.img_type
        self.process.in_files = self.in_file
        _, fileIval = os.path.split(self.in_file)
        self.process.out_basename = os.path.join(
            self.output_directory, fileIval
        )
        self.process.output_type = self.output_type
        self.process.segments = self.segments

        return self.process.run(configuration_dict={})


class Smooth(ProcessMIA):
    """
    *3D Gaussian smoothing of image volumes*

    Please, see the complete documentation for the `Smooth' brick in
    the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/fsl/Smooth.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Smooth, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["fsl", "nipype"]

        # Mandatory inputs description
        in_file_desc = (
            "A file to smooth (a pathlike object or string "
            "representing a file)."
        )
        # Optional inputs with default value description
        fwhm_desc = (
            "Gaussian kernel fwhm in mm (a float). "
            "Mutually exclusive with sigma. "
            "Basically, 2.3548 * sigma = fwhm. "
            "Default is 6"
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the smoothed image file(s) "
            "(a string)."
        )
        output_type_desc = (
            "Typecodes of the output NIfTI image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        # Optional inputs description
        sigma_desc = (
            "Gaussian kernel sigma in mm (a float). Mutually "
            "exclusive with fwhm. Basically, 2.3548 * sigma = fwhm."
        )
        # Outputs description
        out_file_desc = (
            "The smoothed files (a pathlike object or a "
            "string representing a file)."
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits
        self.add_trait(
            "fwhm",
            Either(
                Undefined,
                Float(),
                default=6.0,
                output=False,
                optional=True,
                desc=fwhm_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String("s_", output=False, optional=True, desc=out_prefix_desc),
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

        # Optional inputs traits
        self.add_trait(
            "sigma",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=sigma_desc,
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        # To suppress the "FSLOUTPUTTYPE environment
        # variable is not set" nipype warning:
        if "FSLOUTPUTTYPE" not in os.environ:
            os.environ["FSLOUTPUTTYPE"] = self.output_type

        self.init_process("nipype.interfaces.fsl.Smooth")

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
        super(Smooth, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.sigma == Undefined and self.fwhm == Undefined:
            print(
                "\nInitialisation failed. Please, set one of the two input "
                "parameters sigma or fwhm ...!"
            )
            return

        elif self.sigma != Undefined and self.fwhm != Undefined:
            print(
                '\nInitialisation failed. Both input parameters "sigma" and '
                '"fwhm" are mutually exclusive. Please, define only one of '
                "these two parameters (set the other as Undefined) ...!"
            )
            return

        if self.in_file:
            if self.out_prefix == Undefined:
                self.out_prefix = "s_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "s" ...'
                )

            if self.output_directory:
                output_type = self.output_type
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is" " not recognized...!")
                    return
                else:
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName + "." + EXT[output_type],
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
        super(Smooth, self).run_process_mia()
        self.process.in_file = self.in_file

        if self.fwhm == Undefined:
            self.process.sigma = self.sigma
            self.process.trait("fwhm").optional = True

        else:
            self.process.fwhm = self.fwhm
            self.process.trait("sigma").optional = True

        self.process.output_type = self.output_type
        self.process.smoothed_file = self.out_file
        return self.process.run(configuration_dict={})
