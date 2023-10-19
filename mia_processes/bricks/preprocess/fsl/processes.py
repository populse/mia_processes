# -*- coding: utf-8 -*-

"""The fsl preprocess library of the mia_processes package.

The purpose of this module is to customise the main fsl preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - BetSurfacesExtraction
        - ExtractROI
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
    in the mia_processes website
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
                print("\nThe input image format is not recognized...!")
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
            for k in (
                "out_file",
                "outskin_mask_file",
                "outskin_mesh_file",
                "outskull_mask_file",
                "outskull_mesh_file",
                "inskull_mask_file",
                "inskull_mesh_file",
                "skull_mask_file",
            ):
                self.tags_inheritance(
                    in_file=self.in_file,
                    out_file=self.outputs[k],
                )

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


class ExtractROI(ProcessMIA):
    """
    *Extract region of interest (ROI) from an image using FslROI*

    Please, see the complete documentation for the `ExtractROI brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/htmldocumentation/bricks/preprocess/fsl/ExtractROI.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ExtractROI, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["fsl", "nipype"]

        # Mandatory inputs description
        in_file_desc = "Input file"
        # Optional inputs with default value description
        suffix_desc = "Output suffix (a string, default value is roi)"
        output_type_desc = (
            "Typecodes of the output NIfTI image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        t_min_desc = "t minimum index (an integer)."
        t_size_desc = "t size (an integer)."
        x_min_desc = "x minimum index (an integer)."
        x_size_desc = "x size (an integer)."
        y_min_desc = "y minimum index (an integer)."
        y_size_desc = "y size (an integer)."
        z_min_desc = "z minimum index (an integer)."
        z_size_desc = "z size (an integer)."

        # Outputs description
        roi_file_desc = (
            "Output file (a pathlike object or string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits

        self.add_trait(
            "suffix",
            String("roi", output=False, optional=True, desc=suffix_desc),
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
            "t_min",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=t_min_desc,
            ),
        )

        self.add_trait(
            "t_size",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=t_size_desc,
            ),
        )

        self.add_trait(
            "x_min",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=x_min_desc,
            ),
        )

        self.add_trait(
            "x_size",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=x_size_desc,
            ),
        )

        self.add_trait(
            "y_min",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=y_min_desc,
            ),
        )

        self.add_trait(
            "y_size",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=y_size_desc,
            ),
        )

        self.add_trait(
            "z_min",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=z_min_desc,
            ),
        )

        self.add_trait(
            "z_size",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=z_size_desc,
            ),
        )

        # Outputs traits
        self.add_trait("roi_file", File(output=True, desc=roi_file_desc))

        self.init_default_traits()

        # To suppress the "FSLOUTPUTTYPE environment
        # variable is not set" nipype warning:
        if "FSLOUTPUTTYPE" not in os.environ:
            os.environ["FSLOUTPUTTYPE"] = self.output_type

        self.init_process("nipype.interfaces.fsl.ExtractROI")

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
        super(ExtractROI, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return

            if self.output_directory:
                self.outputs["roi_file"] = os.path.join(
                    self.output_directory,
                    fileName + "_" + self.suffix + "." + EXT[self.output_type],
                )
            else:
                print("No output_directory was found...!\n")
                return

        if self.outputs:
            if self.roi_file:
                self.tags_inheritance(
                    in_file=self.in_file,
                    out_file=self.outputs["roi_file"],
                )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ExtractROI, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.output_type = self.output_type
        self.process.roi_file = self.roi_file
        if self.t_min is not Undefined:
            self.process.t_min = self.t_min
        if self.t_size is not Undefined:
            self.process.t_size = self.t_size
        if self.x_min is not Undefined:
            self.process.x_min = self.x_min
        if self.x_size is not Undefined:
            self.process.x_size = self.x_size
        if self.y_min is not Undefined:
            self.process.y_min = self.y_min
        if self.y_size is not Undefined:
            self.process.y_size = self.y_size
        if self.z_min is not Undefined:
            self.process.z_min = self.z_min
        if self.z_size is not Undefined:
            self.process.z_size = self.z_size

        return self.process.run(configuration_dict={})


class FastSegment(ProcessMIA):
    """
    *Brain tissue segmentation using fsl.FAST*

    Please, see the complete documentation for the `FastSegment brick
    in the mia_processes website
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
            "Int specifying type of image: (1 = T1, 2 = T2, 3 = PD). "
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
            "Partial volume files (a list of items which are file names."
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
                print("\nThe input image format is not recognized...!")
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
            for k in ("tissue_class_map", "partial_volume_map", "mixeltype"):
                self.tags_inheritance(
                    in_file=self.in_file,
                    out_file=self.outputs[k],
                )

            if self.outputs["partial_volume_files"][0]:
                for out_val in self.outputs["partial_volume_files"]:
                    self.tags_inheritance(
                        in_file=self.in_file,
                        out_file=out_val,
                    )

            if self.outputs["tissue_class_files"][0]:
                for out_val in self.outputs["tissue_class_files"]:
                    self.tags_inheritance(
                        in_file=self.in_file,
                        out_file=out_val,
                    )

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


class Flirt(ProcessMIA):
    """
    *Linear (affine) intra and inter-modal brain image registration
    using fsl FLIRT*

    Please, see the complete documentation for the `Flirt brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/htmldocumentation/bricks/preprocess/fsl/Flirt.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Flirt, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["fsl", "nipype"]

        # Mandatory inputs description
        in_file_desc = (
            "Input file (a pathlike object"
            "string representing an existing file)"
        )
        in_reference_file_desc = (
            "Reference file (a pathlike object"
            "string representing an existing file)"
        )
        # Optional inputs with default value description
        apply_xfm_desc = (
            "Apply transformation supplied by in_matrix_file "
            "(a boolean). Mutually exclusive with apply_isoxfm. "
            "Required in_matrix_file"
        )
        apply_isoxfm_desc = (
            "Apply transformation supplied by in_matrix_file but "
            "forces isotropic resampling (a float). "
            "Mutually exclusive with apply_xfm."
            "Required in_matrix_file"
        )
        in_matrix_file_desc = (
            "Input 4x4 affine matrix. (a pathlike object or string "
            "representing a file)"
        )
        get_registered_file_desc = "Get registered file (a boolean)"
        angle_rep_desc = (
            "Representation of rotation angles (quaternion or euler)"
            "Default is euler"
        )
        bbrslope_desc = "Value of bbr slope. (a float)"
        bbrtype_desc = (
            "Type of bbr cost function (signed, global_abs orlocal_abs)"
            "Default is signed"
        )
        bgvalue_desc = (
            "Use specified background value for points outside FOV. (a float)"
        )
        bins_desc = "Number of histogram bins (an integer). Default is 256"
        coarse_search_desc = (
            "Coarse search delta angle (an integer). Default is 60"
        )
        cost_desc = (
            "Cost function (mutualinfo or corratio or normcorr or normmi or "
            "leastsq or labeldiff or bbr). Default is corratio"
        )
        cost_func_desc = (
            "Cost function (searchcost) (mutualinfo or corratio or normcorr "
            "or normmi or leastsq or labeldiff or bbr). Default is corratio"
        )
        datatype_desc = (
            "Force output data type (char, short, int, float, double)"
        )
        dof_desc = "Number of transform degrees of freedom (an integer) "
        echospacing_desc = (
            "Value of EPI echo spacing - units of seconds. (a float)"
        )
        fine_search_desc = (
            "Fine search delta angle (an integer, default is 18)"
        )
        filedmap_desc = (
            "Fieldmap image in rads/s - must be already registered "
            "to the reference image (a pathlike object or string "
            "representing a file)"
        )
        filedmapmask_desc = (
            "Mask for fieldmap image (a pathlike object or string "
            "representing a file)"
        )
        force_scaling_desc = (
            "Force rescaling even for low-res images. (a boolean)"
        )
        interp_desc = (
            "Final interpolation method used in reslicing. "
            "(trilinear or nearestneighbour or sinc or spline)"
        )
        in_weight_desc = (
            "File for input weighting volume (a pathlike object or"
            "string representing an existing file)"
        )
        min_sampling_desc = (
            "Set minimum voxel dimension for sampling (a float)"
        )
        no_clamp_desc = "Do not use intensity clampinp (a boolean)"
        no_resample_desc = "Do not change input sampling (a boolean)"
        no_resample_blur_desc = (
            "Do not use blurring on downsampling (a boolean)"
        )
        no_search_desc = (
            "Set all angular searches to ranges 0 to 0 (a boolean)"
        )
        output_type_desc = (
            "Typecodes of the output NIfTI image formats (one "
            "of NIFTI, NIFTI_GZ)."
        )
        padding_size_desc = (
            "For applyxfm: interpolates outside image by size (an integer)"
        )
        pedir_desc = (
            "Phase encode direction of EPI - 1/2/3=x/y/z & -1/-2/-3=-x/-y/-z "
            "(an integer)"
        )
        rigid2D_desc = "Use 2D rigid body mode ie ignore dof (a boolean)"
        ref_weight_desc = (
            "File for reference weighting volume (a pathlike "
            "object or string representing an existing file) "
        )
        save_log_desc = "Save log (a boolean)"
        schedule_desc = (
            "Replaces default schedule (a pathlike object or "
            "string representing an existing file) "
        )
        searchr_x_desc = (
            "Search angles along x-axis, in degrees (a list of two items "
            "which are integer [min angle, max angle], default is [-90, 90])"
        )
        searchr_y_desc = (
            "Search angles along y-axis, in degrees (a list of two items "
            "which are integer [min angle, max angle], default is [-90, 90])"
        )
        searchr_z_desc = (
            "Search angles along z-axis, in degrees (a list of two items "
            "which are integer [min angle, max angle], default is [-90, 90])"
        )
        sinc_width_desc = "Full-width in voxels (an integer). Default is 7."
        sinc_window_desc = "Sinc window (rectangular, hanning, blackman)"
        uses_qform_desc = "Initialize using sform or qform (a boolean)"
        wm_seg_desc = (
            "White matter segmentation volume needed by BBR cost function"
            "(a pathlike object or string representing a file)"
        )
        wmcoords_desc = (
            "White matter boundary coordinates for BBR cost function"
            "(a pathlike object or string representing a file)"
        )
        wmnorms_desc = (
            "White matter boundary normals for BBR cost function"
            "(a pathlike object or string representing a file)"
        )

        # Outputs description
        out_file_desc = (
            "Path/name of registered file (if generated)"
            "(a pathlike object or string representing"
            "a file)"
        )
        out_log_desc = (
            "Path/name of output log (if generated)"
            "(a pathlike object or string representing"
            "a file)"
        )
        out_matrix_file_desc = (
            "Path/name of calculated affine transform (if generated)"
            "(a pathlike object or string representing"
            "a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "in_reference_file",
            File(output=False, optional=False, desc=in_reference_file_desc),
        )

        # Optional inputs with default value traits
        self.add_trait(
            "get_registered_file",
            Bool(
                True,
                output=False,
                optional=True,
                desc=get_registered_file_desc,
            ),
        )

        self.add_trait(
            "apply_xfm",
            Bool(
                False,
                output=False,
                optional=True,
                desc=apply_xfm_desc,
            ),
        )

        self.add_trait(
            "apply_isoxfm",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=apply_isoxfm_desc,
            ),
        )

        self.add_trait(
            "in_matrix_file",
            File(output=False, optional=True, desc=in_matrix_file_desc),
        )

        self.add_trait(
            "angle_rep",
            Enum(
                "euler",
                "quaternion",
                output=False,
                optional=True,
                desc=angle_rep_desc,
            ),
        )

        self.add_trait(
            "bbrslope",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=bbrslope_desc,
            ),
        )

        self.add_trait(
            "bbrtype",
            Enum(
                "signed",
                "global_abs",
                "local_abs",
                output=False,
                optional=True,
                desc=bbrtype_desc,
            ),
        )

        self.add_trait(
            "bgvalue",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=bgvalue_desc,
            ),
        )

        self.add_trait(
            "bins",
            Int(
                256,
                output=False,
                optional=True,
                desc=bins_desc,
            ),
        )

        self.add_trait(
            "coarse_search",
            Int(
                60,
                output=False,
                optional=True,
                desc=coarse_search_desc,
            ),
        )

        self.add_trait(
            "cost",
            Enum(
                "corratio",
                "mutualinfo",
                "normcorr",
                "normmi",
                "leastsq",
                "labeldiff",
                "bbr",
                output=False,
                optional=True,
                desc=cost_desc,
            ),
        )

        self.add_trait(
            "cost_func",
            Enum(
                "corratio",
                "mutualinfo",
                "normcorr",
                "normmi",
                "leastsq",
                "labeldiff",
                "bbr",
                output=False,
                optional=True,
                desc=cost_func_desc,
            ),
        )

        self.add_trait(
            "datatype",
            Either(
                Enum("char", "short", "int", "float", "double"),
                Undefined,
                default=Undefined,
                output=False,
                optional=True,
                desc=datatype_desc,
            ),
        )

        self.add_trait(
            "dof",
            Int(
                12,
                output=False,
                optional=True,
                desc=dof_desc,
            ),
        )

        self.add_trait(
            "echospacing",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=echospacing_desc,
            ),
        )

        self.add_trait(
            "fieldmap",
            File(output=False, optional=True, desc=filedmap_desc),
        )

        self.add_trait(
            "fieldmapmask",
            File(output=False, optional=True, desc=filedmapmask_desc),
        )

        self.add_trait(
            "fine_search",
            Int(
                18,
                output=False,
                optional=True,
                desc=fine_search_desc,
            ),
        )

        self.add_trait(
            "force_scaling",
            Bool(
                False,
                output=False,
                optional=True,
                desc=force_scaling_desc,
            ),
        )

        self.add_trait(
            "interp",
            Enum(
                "trilinear",
                "nearestneighbour",
                "sinc",
                "spline",
                output=False,
                optional=True,
                desc=interp_desc,
            ),
        )

        self.add_trait(
            "in_weight",
            File(output=False, optional=True, desc=in_weight_desc),
        )

        self.add_trait(
            "min_sampling",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=min_sampling_desc,
            ),
        )

        self.add_trait(
            "no_clamp",
            Bool(
                False,
                output=False,
                optional=True,
                desc=no_clamp_desc,
            ),
        )

        self.add_trait(
            "no_resample",
            Bool(
                False,
                output=False,
                optional=True,
                desc=no_resample_desc,
            ),
        )

        self.add_trait(
            "no_resample_blur",
            Bool(
                False,
                output=False,
                optional=True,
                desc=no_resample_blur_desc,
            ),
        )

        self.add_trait(
            "no_search",
            Bool(
                False,
                output=False,
                optional=True,
                desc=no_search_desc,
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
            "padding_size",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=padding_size_desc,
            ),
        )

        self.add_trait(
            "pedir",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=pedir_desc,
            ),
        )

        self.add_trait(
            "ref_weight",
            File(output=False, optional=True, desc=ref_weight_desc),
        )

        self.add_trait(
            "rigid2D",
            Bool(
                False,
                output=False,
                optional=True,
                desc=rigid2D_desc,
            ),
        )

        self.add_trait(
            "save_log",
            Bool(
                False,
                output=False,
                optional=True,
                desc=save_log_desc,
            ),
        )

        self.add_trait(
            "searchr_x",
            List(
                [-90, 90],
                output=False,
                optional=True,
                desc=searchr_x_desc,
            ),
        )

        self.add_trait(
            "searchr_y",
            List(
                [-90, 90],
                output=False,
                optional=True,
                desc=searchr_y_desc,
            ),
        )

        self.add_trait(
            "searchr_z",
            List(
                [-90, 90],
                output=False,
                optional=True,
                desc=searchr_z_desc,
            ),
        )

        self.add_trait(
            "schedule",
            File(output=False, optional=True, desc=schedule_desc),
        )

        self.add_trait(
            "sinc_width",
            Either(
                Int(),
                Undefined,
                default=Undefined,
                output=False,
                optional=True,
                desc=sinc_width_desc,
            ),
        )

        self.add_trait(
            "sinc_window",
            Either(
                Enum("rectangular", "hanning", "blackman"),
                Undefined,
                default=Undefined,
                output=False,
                optional=True,
                desc=sinc_window_desc,
            ),
        )

        self.add_trait(
            "uses_qform",
            Bool(
                False,
                output=False,
                optional=True,
                desc=uses_qform_desc,
            ),
        )

        self.add_trait(
            "wm_seg",
            File(output=False, optional=True, desc=wm_seg_desc),
        )

        self.add_trait(
            "wmcoords",
            File(output=False, optional=True, desc=wmcoords_desc),
        )

        self.add_trait(
            "wmnorms",
            File(output=False, optional=True, desc=wmnorms_desc),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=True, desc=out_file_desc)
        )

        self.add_trait(
            "out_log",
            File(output=True, optional=True, desc=out_log_desc),
        )

        self.add_trait(
            "out_matrix_file",
            File(output=True, optional=True, desc=out_matrix_file_desc),
        )

        self.init_default_traits()

        # To suppress the "FSLOUTPUTTYPE environment
        # variable is not set" nipype warning:
        if "FSLOUTPUTTYPE" not in os.environ:
            os.environ["FSLOUTPUTTYPE"] = self.output_type

        self.init_process("nipype.interfaces.fsl.FLIRT")

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
        super(Flirt, self).list_outputs()

        if self.apply_isoxfm and self.apply_xfm:
            print(
                "\nInitialisation failed. Both input parameters apply_isoxfm "
                "and apply_xfm are mutually exclusive. Please, define only "
                "one of these two parameters...!"
            )
            return self.make_initResult()
        if (self.apply_isoxfm or self.apply_xfm) and not self.in_matrix_file:
            print(
                "\nInitialisation failed. Parameters apply_isoxfm or "
                "apply_xfm required in_matrix_file. Please, define "
                "in_matrix_file "
            )
            return self.make_initResult()
        if self.padding_size and not self.apply_xfm:
            print(
                "\nInitialisation failed. Parameters padding_size "
                "required apply_xfm. "
            )
            return self.make_initResult()

        # Outputs definition and tags inheritance (optional)
        if self.in_file and self.in_reference_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)
            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return

            valid_ext_ref, in_ext_ref, fileName_ref = checkFileExt(
                self.in_reference_file, EXT
            )
            if not valid_ext_ref:
                print(
                    "\nThe input reference image format is "
                    "not recognized...!"
                )
                return

            if self.output_directory:
                if self.get_registered_file:
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        fileName
                        + "_registered_with_"
                        + fileName_ref
                        + "."
                        + in_ext,
                    )
                if not self.apply_xfm or not self.apply_isoxfm:
                    self.outputs["out_matrix_file"] = os.path.join(
                        self.output_directory,
                        fileName + "_" + fileName_ref + "_flirt.mat",
                    )
                if self.save_log:
                    self.outputs["out_log"] = os.path.join(
                        self.output_directory,
                        fileName + "_flirt_log.txt",
                    )

            else:
                print("No output_directory was found...!\n")
                return

        if self.outputs:
            if self.save_log:
                self.tags_inheritance(
                    in_file=self.in_file,
                    out_file=self.outputs["out_log"],
                )
            if self.get_registered_file:
                self.tags_inheritance(
                    in_file=self.in_file,
                    out_file=self.outputs["out_file"],
                )

            self.tags_inheritance(
                in_file=self.in_file,
                out_file=self.outputs["out_matrix_file"],
            )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Flirt, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.reference = self.in_reference_file
        if self.out_file:
            self.process.out_file = self.out_file
        if self.out_matrix_file:
            self.process.out_matrix_file = self.out_matrix_file
        if self.apply_xfm:
            self.process.apply_xfm = self.apply_xfm
        if self.apply_isoxfm:
            self.process.apply_isoxfm = self.apply_isoxfm
        if self.in_matrix_file:
            self.process.in_matrix_file = self.in_matrix_file
        self.process.angle_rep = self.angle_rep
        if self.bbrslope:
            self.process.bbrslope = self.bbrslope
        self.process.bbrtype = self.bbrtype
        if self.bgvalue:
            self.process.bbrslope = self.bgvalue
        self.process.bins = self.bins
        self.process.coarse_search = self.coarse_search
        self.process.cost = self.cost
        self.process.cost_func = self.cost_func
        if self.datatype:
            self.process.datatype = self.datatype
        self.process.dof = self.dof
        if self.echospacing:
            self.process.echospacing = self.echospacing
        if self.fieldmap:
            self.process.fieldmap = self.fieldmap
        if self.fieldmapmask:
            self.process.fieldmapmask = self.fieldmapmask
        self.process.fine_search = self.fine_search
        self.process.force_scaling = self.force_scaling
        if self.in_weight:
            self.process.in_weight = self.in_weight
        self.process.interp = self.interp
        if self.min_sampling:
            self.process.min_sampling = self.min_sampling
        self.process.no_clamp = self.no_clamp
        self.process.no_resample = self.no_resample
        self.process.no_resample_blur = self.no_resample_blur
        self.process.no_search = self.no_search
        self.process.output_type = self.output_type
        if self.padding_size:
            self.process.padding_size = self.padding_size
        if self.pedir:
            self.process.pedir = self.pedir
        if self.ref_weight:
            self.process.ref_weight = self.ref_weight
        self.process.rigid2D = self.rigid2D
        if self.save_log:
            self.process.save_log = self.save_log
            self.process.out_log = self.out_log
        self.process.searchr_x = self.searchr_x
        self.process.searchr_y = self.searchr_y
        self.process.searchr_z = self.searchr_z
        if self.schedule:
            self.process.schedule = self.schedule
        self.process.sinc_width = self.sinc_width
        if self.sinc_window:
            self.process.sinc_window = self.sinc_window
        self.process.uses_qform = self.uses_qform
        if self.wm_seg:
            self.process.wm_seg = self.wm_seg
        if self.wmcoords:
            self.process.wmcoords = self.wmcoords
        if self.wmnorms:
            self.process.wmnorms = self.wmnorms

        return self.process.run(configuration_dict={})


class Smooth(ProcessMIA):
    """
    *3D Gaussian smoothing of image volumes*

    Please, see the complete documentation for the `Smooth' brick in
    the mia_processes website
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
                    print("\nThe input image format is not recognized...!")
                    return
                else:
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName + "." + EXT[output_type],
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
