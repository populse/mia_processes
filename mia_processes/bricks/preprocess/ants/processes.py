# -*- coding: utf-8 -*-

"""The ants preprocess library of the mia_processes package.

The purpose of this module is to customise the main ants preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - AffineInitializer
        - ApplyTransforms
        - N4BiasFieldCorrection
        - Registration

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
from nipype.interfaces.base import File, InputMultiPath

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA
from traits.api import (
    Bool,
    Either,
    Enum,
    Float,
    Int,
    List,
    Range,
    String,
    Tuple,
    Undefined,
)

from mia_processes.utils import checkFileExt

EXT = {"NIFTI_GZ": "nii.gz", "NIFTI": "nii"}


class AffineInitializer(ProcessMIA):
    """
    *A multi-start optimizer for affine registration*

    Please, see the complete documentation for the `AffineInitializer brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/ants/AffineInitializer.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(AffineInitializer, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["ants", "nipype"]

        # Mandatory inputs description
        moving_image_desc = "The moving image to be mapped to the fixed space."
        fixed_image_desc = (
            "The fixed reference image (a pathlike object or "
            "string representing a file)."
        )
        # Optional inputs with default value description
        dimension_desc = "Image dimension (3 or 2). Default is 3."
        local_search_desc = (
            "Determines if a local optimization is run at"
            "each search point for the set number of "
            "iterations (an integer). Default is 10."
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the corrected image file "
            "(a string)."
        )
        principal_axes_desc = (
            "Whether the rotation is searched around"
            "an initial principal axis alignment. "
            "Default is False."
        )
        radian_fraction_desc = (
            "Search this arc +/- principal axes "
            "(a float between 0.0 and 1.0)."
            "Default is 0.1"
        )
        search_factor_desc = (
            "Increments (degrees) for affine search"
            "(a float). Default is 15.0"
        )

        # Outputs description
        out_file_desc = (
            "Output transform file (a pathlike object or a "
            "string representing a file)."
        )

        # Mandatory inputs traits
        self.add_trait(
            "moving_image",
            File(output=False, optional=False, desc=moving_image_desc),
        )

        self.add_trait(
            "fixed_image",
            File(output=False, optional=False, desc=fixed_image_desc),
        )

        # Optional inputs with default value traits
        self.add_trait(
            "dimension",
            Enum(3, 2, output=False, optional=True, desc=dimension_desc),
        )

        self.add_trait(
            "local_search",
            Int(10, output=False, optional=True, desc=local_search_desc),
        )

        self.add_trait(
            "out_prefix",
            String(
                "AffineTransform_",
                output=False,
                optional=True,
                desc=out_prefix_desc,
            ),
        )

        self.add_trait(
            "principal_axes",
            Bool(False, output=False, optional=True, desc=principal_axes_desc),
        )

        self.add_trait(
            "radian_fraction",
            Range(
                value=0.1,
                low=0.0,
                high=1.0,
                output=False,
                optional=True,
                desc=radian_fraction_desc,
            ),
        )

        self.add_trait(
            "search_factor",
            Float(15.0, output=False, optional=True, desc=search_factor_desc),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.ants.AffineInitializer")

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
        super(AffineInitializer, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.moving_image:
            if self.out_prefix == Undefined:
                self.out_prefix = "AffineTransform_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "AffineTransform_" ...'
                )

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(
                    self.moving_image, EXT
                )

                if not valid_ext:
                    print("\nThe input image format is not recognized ...!")
                    return
                else:
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName + ".mat",
                    )

                    self.inheritance_dict[
                        self.outputs["out_file"]
                    ] = self.moving_image

            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(AffineInitializer, self).run_process_mia()

        # TODO: We see in soma_workflow: nipype.interface INFO:
        #       bad det -1 v 1 u -1

        # Mandatory inputs
        self.process.fixed_image = self.fixed_image
        self.process.moving_image = self.moving_image
        # Others inputs
        self.process.dimension = self.dimension
        self.process.local_search = self.local_search
        self.process.out_file = self.out_file
        self.process.principal_axes = self.principal_axes
        self.process.radian_fraction = self.radian_fraction
        self.process.search_factor = self.search_factor

        return self.process.run(configuration_dict={})


class ApplyTransforms(ProcessMIA):
    """
    *Transforms an image according to a reference image and transformation \
    (or set of transformation)*

    Please, see the complete documentation for the `ApplyTransforms brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/ants/ApplyTransforms.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ApplyTransforms, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["ants", "nipype"]

        # Mandatory inputs description
        input_image_desc = (
            "Image to apply transformation to. (a pathlike "
            "object or string representing an existing file) "
        )
        reference_image_desc = (
            "Reference image space that you wish to warp "
            "into (a pathlike object or string "
            " representing an existing file)"
        )
        transforms_desc = (
            "Transform files that will be applied in "
            "reverse order. (ie the last specified "
            "transform will be applied first)"
        )
        # Optional inputs with default value description
        default_value_desc = "Default value. Default is 0.0."
        dimension_desc = (
            "This option forces the image to be treated as a "
            "specified-dimensional image (2 or 3 or 4). "
            "Default is 3"
        )
        float_desc = (
            "Use float instead of double for computations. "
            "Default is False."
        )
        input_image_type_desc = (
            "Option specifying the input image type of "
            "scalar (0), vector (1), tensor(2), or "
            "time series(3). Default is 0 (scalar)."
        )
        interpolation_desc = (
            "Choice of interpolator. "
            "(‘Linear’ or ‘NearestNeighbor’ "
            "or ‘CosineWindowedSinc’ or ‘WelchWindowedSinc’ "
            "or ‘HammingWindowedSinc’ or "
            "‘LanczosWindowedSinc’ or ‘MultiLabel’ or "
            "‘Gaussian’ or ‘BSpline’)"
        )
        print_out_composite_warp_file_desc = (
            "Output a composite warp file "
            "instead of a transformed "
            "image. (a boolean)"
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the corrected image file "
            "(a string)."
        )
        # Optional inputs description
        interpolation_parameters_desc = (
            "Interpolation parameters for "
            "example for BSpline order or for "
            "sigma/alphaparameters for "
            "Gaussian/Multilabel (a tuple "
            "of the form: (an integer) "
            "or a tuple of the form: "
            "(a float, a float))"
        )
        invert_transform_flags_desc = (
            "Invert transform flags " "(a list of boolean)"
        )

        # Outputs description
        output_image_desc = (
            "Warped image (a pathlike object or string"
            "representing an existing file)."
        )

        # Mandatory inputs traits
        self.add_trait(
            "input_image",
            File(output=False, optional=False, desc=input_image_desc),
        )

        self.add_trait(
            "reference_image",
            File(output=False, optional=False, desc=reference_image_desc),
        )

        self.add_trait(
            "transforms",
            InputMultiPath(
                File(), output=False, optional=False, desc=transforms_desc
            ),
        )
        # Optional inputs with default value traits
        self.add_trait(
            "default_value",
            Float(0.0, output=False, optional=True, desc=default_value_desc),
        )

        self.add_trait(
            "dimension",
            Enum(3, 2, 4, output=False, optional=True, desc=dimension_desc),
        )

        self.add_trait(
            "float", Bool(True, output=False, optional=True, desc=float_desc)
        )

        self.add_trait(
            "input_image_type",
            Enum(
                0,
                1,
                2,
                3,
                output=False,
                optional=True,
                desc=input_image_type_desc,
            ),
        )

        self.add_trait(
            "interpolation",
            Enum(
                "Linear",
                "NearestNeighbor",
                "CosineWindowedSinc",
                "WelchWindowedSinc",
                "HammingWindowedSinc",
                "LanczosWindowedSinc",
                "MultiLabel",
                "Gaussian" "BSpline",
                default="Linear",
                output=False,
                optional=True,
                desc=interpolation_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String("t_", output=False, optional=True, desc=out_prefix_desc),
        )

        self.add_trait(
            "print_out_composite_warp_file",
            Bool(
                False,
                output=False,
                optional=True,
                desc=print_out_composite_warp_file_desc,
            ),
        )

        # Optional inputs with default value traits
        self.add_trait(
            "interpolation_parameters",
            Either(
                Undefined,
                Tuple(Int()),
                Tuple(Float(), Float()),
                default=Undefined,
                output=False,
                optional=True,
                desc=interpolation_parameters_desc,
            ),
        )

        self.add_trait(
            "invert_transform_flags",
            Either(
                Undefined,
                List(Bool()),
                default=Undefined,
                output=False,
                optional=True,
                desc=invert_transform_flags_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "output_image", File(output=True, desc=output_image_desc)
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.ants.ApplyTransforms")

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
        super(ApplyTransforms, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.input_image:
            if self.out_prefix == Undefined:
                self.out_prefix = "t_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "t_" ...'
                )

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(
                    self.input_image, EXT
                )

                if not valid_ext:
                    print("\nThe input image format is not recognized ...!")
                    return
                else:
                    self.outputs["output_image"] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName + "." + in_ext,
                    )

                    self.inheritance_dict[
                        self.outputs["output_image"]
                    ] = self.input_image

            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ApplyTransforms, self).run_process_mia()

        self.process.input_image = self.input_image
        self.process.reference_image = self.reference_image
        self.process.transforms = self.transforms

        self.process.dimension = self.dimension
        self.process.default_value = self.default_value
        self.process.float = self.float
        self.process.input_image_type = self.input_image_type
        self.process.interpolation = self.interpolation
        if self.interpolation_parameters:
            self.process.interpolation_parameters = (
                self.interpolation_parameters
            )
        if self.invert_transform_flags:
            self.process.invert_transform_flags = self.invert_transform_flags
        if self.process.print_out_composite_warp_file:
            self.process.print_out_composite_warp_file = (
                self.print_out_composite_warp_file
            )

        self.process.output_image = self.output_image

        return self.process.run(configuration_dict={})


class N4BiasFieldCorrection(ProcessMIA):
    """
    *N4 Bias field correction*

    Please, see the complete documentation for the `N4BiasFieldCorrection
    brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/ants/N4BiasFieldCorrection.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(N4BiasFieldCorrection, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["ants", "nipype"]

        # Mandatory inputs description
        in_file_desc = (
            "A image file (2D, 3D or 4D) to correct (a pathlike"
            "object or string representing a file)."
        )
        # Optional inputs with default value description
        copy_header_desc = (
            " Copy headers of the original image into the "
            "output (corrected) file (a boolean). "
            "Default is True"
        )
        dimension_desc = "Image dimension (2, 3 or 4). Default is 3"
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filenames of the corrected image file "
            "(a string)."
        )
        rescale_intensities_desc = (
            "Rescales to the [min,max] range of the"
            "original image intensities within the "
            "user-specified mask (a boolean)."
        )
        save_bias_desc = (
            "Save the estimated bias in a file (a boolean). " "Default is True"
        )
        # Optional inputs description
        bspline_fitting_distance_desc = (
            "Set bspline fitting distance " "(a float)"
        )
        bspline_order_desc = (
            "bspline order (an integer). "
            "Requires bspline_fitting_distance parameter"
        )
        convergence_threshold_desc = (
            "Convergence threshodl (a float). "
            "Require n_iterations parameter"
        )
        histogram_sharpening_desc = (
            "Three-values tuple of histogram "
            "sharpening parameters (FWHM, "
            "wienerNose, numberOfHistogramBins). "
            "(a tuple of the form: "
            "(a float, a float, an integer))"
        )
        mask_image_desc = (
            "Image to specify region to perform final "
            "bias correction in. (a pathlike"
            "object or string representing a file)"
        )
        n_iterations_desc = (
            "Number of iterations (convergence option) " "(a list of integer)"
        )
        shrink_factor_desc = "Shrink factor (an integer)."
        weight_image_desc = (
            "Image for relative weighting (e.g. "
            "probability map of the white matter) of voxels "
            "during the B-spline fitting (a pathlike"
            "object or string representing a file)"
        )
        # Outputs description
        bias_image_desc = (
            "Estimated bias (a pathlike object or a "
            "string representing a file)."
        )
        negative_values_desc = (
            "True if negative values are present in" "False otherwise."
        )
        out_file_desc = (
            "The corrected file (a pathlike object or a "
            "string representing a file)."
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs whith default value traits
        self.add_trait(
            "copy_header",
            Bool(True, output=False, optional=True, desc=copy_header_desc),
        )

        self.add_trait(
            "dimension",
            Enum(3, 2, 4, output=False, optional=True, desc=dimension_desc),
        )

        self.add_trait(
            "out_prefix",
            String("n4c_", output=False, optional=True, desc=out_prefix_desc),
        )

        self.add_trait(
            "rescale_intensities",
            Bool(
                False,
                output=False,
                optional=True,
                desc=rescale_intensities_desc,
            ),
        )

        self.add_trait(
            "save_bias",
            Bool(True, output=False, optional=True, desc=save_bias_desc),
        )

        # Optional inputs traits
        self.add_trait(
            "bspline_fitting_distance",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=bspline_fitting_distance_desc,
            ),
        )

        self.add_trait(
            "bspline_order",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=bspline_order_desc,
            ),
        )

        self.add_trait(
            "convergence_threshold",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=convergence_threshold_desc,
            ),
        )

        self.add_trait(
            "histogram_sharpening",
            Either(
                Undefined,
                Tuple(Float(), Float(), Int()),
                default=Undefined,
                output=False,
                optional=True,
                desc=histogram_sharpening_desc,
            ),
        )

        self.add_trait(
            "mask_image",
            Either(
                Undefined,
                String(),
                default=Undefined,
                output=False,
                optional=True,
                desc=mask_image_desc,
            ),
        )

        self.add_trait(
            "n_iterations",
            Either(
                Undefined,
                List(Int()),
                default=Undefined,
                output=False,
                optional=True,
                desc=n_iterations_desc,
            ),
        )

        self.add_trait(
            "shrink_factor",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=shrink_factor_desc,
            ),
        )

        self.add_trait(
            "weight_image",
            Either(
                Undefined,
                String(),
                default=Undefined,
                output=False,
                optional=True,
                desc=weight_image_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "bias_image",
            File(output=True, optional=True, desc=bias_image_desc),
        )

        self.add_trait(
            "negative_values",
            Bool(output=True, optional=True, desc=negative_values_desc),
        )

        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

        self.init_process("nipype.interfaces.ants.N4BiasFieldCorrection")

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
        super(N4BiasFieldCorrection, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.out_prefix == Undefined:
                self.out_prefix = "n4c_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "n4c" ...'
                )
            if (
                self.bspline_order
                and self.bspline_fitting_distance == Undefined
            ):
                print(
                    'Initialization failed.. "bspline_order" parameter '
                    'required "bspline_fitting_distance" parameter'
                )
                return

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is not recognized ...!")
                    return
                else:
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName + "." + in_ext,
                    )

                    if self.save_bias:
                        self.outputs["bias_image"] = os.path.join(
                            self.output_directory,
                            self.out_prefix + fileName + "_bias." + in_ext,
                        )

                        self.inheritance_dict[
                            self.outputs["out_file"]
                        ] = self.in_file

                    self.inheritance_dict[
                        self.outputs["bias_image"]
                    ] = self.in_file

            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(N4BiasFieldCorrection, self).run_process_mia()

        # If negative values, scale image
        input_nii = nib.load(self.in_file)
        data = input_nii.get_fdata()
        data_min = data.min()

        if data_min < 0:
            data_new = data - data_min
            new_nii = input_nii.__class__(
                data_new, input_nii.affine, input_nii.header
            )
            file_name, trail = self.in_file.rsplit(".", 1)
            nib.save(new_nii, file_name + "_scaled." + trail)
            self.in_file = file_name + "_scaled." + trail
            self.negative_values = True

            print(
                "\nThe input image contains negative values "
                "A scaled copy was created and set as input"
            )

        self.process.input_image = self.in_file
        self.process.copy_header = self.copy_header
        self.process.dimension = self.dimension
        if self.rescale_intensities:
            self.process.rescale_intensities = self.rescale_intensities
        if self.save_bias:
            self.process.bias_image = self.bias_image
        else:
            self.process.save_bias = self.save_bias
        if self.bspline_fitting_distance:
            self.process.bspline_fitting_distance = (
                self.bspline_fitting_distance
            )
        if self.bspline_order:
            self.process.bspline_order = self.bspline_order
        if self.convergence_threshold:
            self.process.convergence_threshold = self.convergence_threshold
        if self.histogram_sharpening:
            self.process.histogram_sharpening = self.histogram_sharpening
        if self.mask_image:
            self.process.mask_image = self.mask_image
        if self.n_iterations:
            self.process.n_iterations = self.n_iterations
        if self.shrink_factor:
            self.process.shrink_factor = self.shrink_factor
        if self.weight_image not in [Undefined, "<undefined>"]:
            self.process.weight_image = self.weight_image
        self.process.output_image = self.out_file

        return self.process.run(configuration_dict={})


class Registration(ProcessMIA):
    """
    *Registers a moving image to a fixed image*

    Please, see the complete documentation for the `Registration brick in
    the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/ants/Registration.html>`

    Here can be found some of classical registration parameters
    (used in fMRIPrep and MRIQC)
    https://github.com/nipreps/niworkflows/tree/master/niworkflows/data

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Registration, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["ants", "nipype"]

        # Mandatory inputs description
        fixed_image_desc = (
            "Image to which the moving image should be "
            "transformed (a pathlike object or string "
            "representing an existing file)."
        )
        moving_image_desc = (
            "Image that will be registered to the space "
            "of the fixed image (a pathlike object"
            "or string representing an existing file)."
        )
        # Optional inputs  whith default value description
        collapse_output_transforms_desc = (
            "Collapse output transforms " "(a boolean). Default is True"
        )
        dimension_desc = "Image dimension (3 or 2). Default is 3."
        float_desc = (
            "Use float instead of double for computations "
            "(a boolean). Default is False"
        )
        initialize_transforms_per_stage_desc = (
            "Initialize linear transforms "
            "from the previous stag "
            "(a boolean). "
            "Default is False"
        )
        interpolation_desc = (
            "(‘Linear’ or ‘NearestNeighbor’ or"
            "‘CosineWindowedSinc’ or ‘WelchWindowedSinc’ "
            "or ‘HammingWindowedSinc’ or"
            "‘LanczosWindowedSinc’ or ‘BSpline’ or"
            "‘MultiLabel’ or ‘Gaussian’ or ‘GenericLabel’)."
            "Default is Linear"
        )
        metric_desc = (
            "The metric to use for each stage (a list of strings"
            "which are "
            "CC"
            " or "
            "MeanSquares"
            " or "
            "Demons"
            ""
            "or "
            "GC"
            " or "
            "MI"
            " or Mattes)"
        )
        metric_weight_desc = (
            "The metric weight(s) for each stage."
            "The weights must sum to 1 per stage"
            "(A list of floats)."
        )
        output_inverse_warped_image_desc = (
            "Get inverse_warped image " "(a boolean). Default is False"
        )
        out_prefix_desc = (
            "Specify the string to be prepended to the "
            "filename of the warped image file "
            "(a string)."
        )
        radius_bins_item_trait_desc = (
            "Radius bins item (an interger). " "Default is 5"
        )
        sigma_units_desc = (
            "Units for smoothing sigmas (mm or vox). "
            'Default is ["vox"] * len(metric)'
        )
        smoothing_sigmas_desc = (
            "Smoothing sigmas (A list of items which "
            "are a list of items (which are a float,)"
        )
        shrink_factors_desc = (
            "(Shrink factor (A list of items which are "
            "a list of items which are an integer)"
        )
        transforms_desc = (
            "(A list of items which are Rigid or Affine or "
            "CompositeAffine or Similarity or Translation or "
            "BSpline or GaussianDisplacementField or "
            "TimeVaryingVelocityField or "
            "TimeVaryingBSplineVelocityField or SyN or "
            "BSplineSyN or Exponential or BSplineExponential"
        )
        transform_parameters_desc = "(A list of tuples)"
        use_estimate_learning_rate_once_desc = (
            "Use estimate learning rate "
            "(a list of boolean). Default "
            "is [True] * len(metric)"
        )
        use_histogram_matching_desc = (
            "Use histogram matching "
            "(a list of boolean). Default "
            "is [True] * len(metric)"
        )
        winsorize_lower_quantile_desc = (
            "The Lower quantile to clip image "
            "ranges (a float between 0.0 and "
            "1.0). Default is 0.005."
        )
        winsorize_upper_quantile_desc = (
            "The upper quantile to clip image "
            "ranges (a float between 0.0 and "
            "1.0). Default is 0.995."
        )
        write_composite_transform_desc = (
            "Write composite transform (a " "boolean). Default is True"
        )
        # Optional inputs description
        convergence_threshold_desc = (
            "Convergence threshold (a list of "
            "at least 1 float)."
            "Requires inputs: number_of_iterations"
        )
        convergence_window_size_desc = (
            "Convergence window size " "(a list of integer)"
        )
        fixed_image_masks_desc = (
            "Mask used to limit metric sampling region "
            "of the fixed image in all stages (a list"
            "of items which are a pathlike object or"
            "string representing an existing file"
            "or ‘NULL’)."
        )
        initial_moving_transform_desc = (
            "A transform or a list of transform "
            "that should be applied before the "
            "registration begin. (a list of items"
            "which are a pathlike object or"
            "string representing an existing"
            "file). Mutually exclusive with "
            "initial_moving_transform_com"
        )
        initial_moving_transform_com_desc = (
            "Align the moving_image and "
            "fixed_image before registration "
            "using the geometric center of "
            "the images (=0), the image "
            "intensities (=1), or the "
            "origin of the images (=2) "
            "(0 or 1 or 2). "
            "Mutually exclusive with "
            "initial_moving_transform"
        )
        interpolation_parameters_desc = (
            "Interpolation parameters for "
            "example for BSpline order or for "
            "sigma/alphaparameters for "
            "Gaussian/Multilabel (a tuple "
            "of the form: (an integer) "
            "or a tuple of the form: "
            "(a float, a float) or "
            "a tuple of the form: (a string))"
        )
        invert_initial_moving_transform_desc = (
            "A list of booleans that "
            "indicatewhether the "
            "inverse(s) of the "
            "transform(s) defined in "
            "initial_moving_transform "
            "should be used. "
            "Mutually exclusive with "
            "initial_moving_transform_com."
            "Required "
            "initial_moving_transform."
        )
        moving_image_masks_desc = (
            "Masks used to limit metric sampling"
            "region of the moving image, defined"
            "per registration stage (Use “NULL” to omit"
            "a mask at a given stage - a list of items"
            "which are a pathlike object or string "
            "representing an existing file or ‘NULL’)."
        )
        number_of_iterations_desc = (
            "Number of iterations. (A list of items "
            "which are a list of items which are an "
            "integer)"
        )
        sampling_percentage_desc = (
            "The metric sampling percentages to use"
            "for each stage (a list of"
            " 0.0 <= floats <=1; requires sampling"
            "strategy)"
        )
        sampling_strategy_desc = (
            "The metric sampling strategies for each"
            "stage (A list of strings which are "
            "None"
            ""
            "or "
            "Regular"
            " or "
            "Random"
            ")."
        )
        radius_or_number_of_bins_desc = (
            "The number of bins in each stage for"
            " the MI  and Mattes metric,"
            "the radius for other metrics"
            "(a list of integers)"
        )
        random_seed_desc = (
            "Fixed seed for random number generation. " "(an integer)"
        )

        # Outputs description
        composite_transform_desc = (
            "Composite transform file"
            "(a pathlike object or a string"
            "representing a file)."
        )
        inverse_composite_transform_desc = (
            "Inverse composite"
            "transform file (a pathlike object"
            "or a string representing"
            "a file)."
        )
        inverse_warped_image_desc = (
            "Inverse warped image. (a pathlike object "
            "or string representing an existing file)."
        )
        warped_image_desc = (
            "Warped image. (a pathlike object "
            "or string representing an existing file)."
        )

        # Mandatory inputs traits
        self.add_trait(
            "fixed_image",
            File(
                value=Undefined,
                output=False,
                optional=False,
                desc=fixed_image_desc,
            ),
        )

        self.add_trait(
            "moving_image",
            File(output=False, optional=False, desc=moving_image_desc),
        )

        # Optional inputs with default value traits

        self.add_trait(
            "collapse_output_transforms",
            Bool(
                True,
                output=False,
                optional=True,
                desc=collapse_output_transforms_desc,
            ),
        )

        self.add_trait(
            "dimension",
            Enum(3, 2, output=False, optional=True, desc=dimension_desc),
        )

        self.add_trait(
            "float", Bool(False, output=False, optional=True, desc=float_desc)
        )

        self.add_trait(
            "initialize_transforms_per_stage",
            Bool(
                False,
                output=False,
                optional=True,
                desc=initialize_transforms_per_stage_desc,
            ),
        )

        self.add_trait(
            "interpolation",
            Enum(
                "Linear",
                "NearestNeighbor",
                "CosineWindowedSinc",
                "WelchWindowedSinc",
                "HammingWindowedSinc",
                "LanczosWindowedSinc",
                "BSpline",
                "MultiLabel",
                "Gaussian",
                "GenericLabel",
                output=False,
                optional=True,
                desc=interpolation_desc,
            ),
        )

        self.add_trait(
            "metric",
            List(
                Enum("CC", "MeanSquares", "Demon", "GC", "MI", "Mattes"),
                default=["Mattes", "Mattes", "Mattes"],
                output=False,
                optional=True,
                desc=metric_desc,
            ),
        )

        self.add_trait(
            "metric_weight",
            List(
                Float(),
                default=[1.0, 1.0, 1.0],
                output=False,
                optional=True,
                desc=metric_weight_desc,
            ),
        )

        self.add_trait(
            "output_inverse_warped_image",
            Bool(
                False,
                output=False,
                optional=True,
                desc=output_inverse_warped_image_desc,
            ),
        )

        self.add_trait(
            "out_prefix",
            String("w_", output=False, optional=True, desc=out_prefix_desc),
        )

        self.add_trait(
            "radius_bins_item_trait",
            Int(
                5,
                output=False,
                optional=True,
                desc=radius_bins_item_trait_desc,
            ),
        )

        self.add_trait(
            "sigma_units",
            List(
                Enum("vox", "mm"),
                default=Undefined,
                output=False,
                optional=True,
                desc=sigma_units_desc,
            ),
        )

        self.add_trait(
            "smoothing_sigmas",
            List(
                List(Float()),
                default=[[4.0], [4.0, 2.0, 0.0], [1.0, 0.0]],
                output=False,
                optional=True,
                desc=smoothing_sigmas_desc,
            ),
        )

        self.add_trait(
            "shrink_factors",
            List(
                List(Int()),
                default=[[4], [4, 2, 1], [2, 1]],
                output=False,
                optional=True,
                desc=shrink_factors_desc,
            ),
        )

        self.add_trait(
            "transforms",
            List(
                Enum(
                    "Rigid",
                    "Affine",
                    "CompositeAffine",
                    "Similarity",
                    "Translation",
                    "BSpline",
                    "GaussianDisplacementField",
                    "TimeVaryingVelocityField",
                    "TimeVaryingBSplineVelocityField",
                    "SyN",
                    "BSplineSyN",
                    "Exponential",
                    "BSplineExponential",
                ),
                default=["Rigid", "Affine", "SyN"],
                output=False,
                optional=True,
                desc=transforms_desc,
            ),
        )

        self.add_trait(
            "transform_parameters",
            List(
                Either(
                    Tuple(Float),
                    Tuple(Float, Float, Float),
                    Tuple(Float, Int, Int, Int),
                    Tuple(Float, Int, Float, Float, Float, Float),
                    Tuple(Float, Float, Float, Int),
                    Tuple(Float, Int, Int, Int, Int),
                ),
                default=[(0.01,), (0.08,), (0.1, 3.0, 0.0)],
                output=False,
                optional=True,
                desc=transform_parameters_desc,
            ),
        )

        self.add_trait(
            "use_estimate_learning_rate_once",
            List(
                Bool(),
                default=Undefined,
                output=False,
                optional=True,
                desc=use_estimate_learning_rate_once_desc,
            ),
        )

        self.add_trait(
            "use_histogram_matching",
            List(
                Bool(),
                default=Undefined,
                output=False,
                optional=True,
                desc=use_histogram_matching_desc,
            ),
        )

        self.add_trait(
            "winsorize_lower_quantile",
            Float(
                0.005,
                min=0.0,
                max=1.0,
                output=False,
                optional=True,
                desc=winsorize_lower_quantile_desc,
            ),
        )

        self.add_trait(
            "winsorize_upper_quantile",
            Float(
                0.995,
                min=0.0,
                max=1.0,
                output=False,
                optional=True,
                desc=winsorize_upper_quantile_desc,
            ),
        )

        self.add_trait(
            "write_composite_transform",
            Bool(
                True,
                output=False,
                optional=True,
                desc=write_composite_transform_desc,
            ),
        )

        # Optional input traits
        self.add_trait(
            "convergence_threshold",
            List(
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=convergence_threshold_desc,
            ),
        )

        self.add_trait(
            "convergence_window_size",
            List(
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=convergence_window_size_desc,
            ),
        )

        self.add_trait(
            "fixed_image_masks",
            Either(
                File(),
                List(File()),
                "NULL",
                default="NULL",
                output=False,
                optional=True,
                desc=fixed_image_masks_desc,
            ),
        )

        self.add_trait(
            "initial_moving_transform",
            InputMultiPath(
                File(),
                default=Undefined,
                output=False,
                optional=True,
                desc=initial_moving_transform_desc,
            ),
        )

        self.add_trait(
            "initial_moving_transform_com",
            Either(
                Undefined,
                Enum(0, 1, 2),
                default=Undefined,
                output=False,
                optional=True,
                desc=initial_moving_transform_com_desc,
            ),
        )

        self.add_trait(
            "interpolation_parameters",
            Either(
                Undefined,
                Tuple(Int()),
                Tuple(Float(), Float()),
                Tuple(String()),
                default=Undefined,
                output=False,
                optional=True,
                desc=interpolation_parameters_desc,
            ),
        )

        self.add_trait(
            "invert_initial_moving_transform",
            Either(
                Undefined,
                List(Bool()),
                default=Undefined,
                output=False,
                optional=True,
                desc=invert_initial_moving_transform_desc,
            ),
        )

        self.add_trait(
            "moving_image_masks",
            Either(
                File(),
                List(File()),
                "NULL",
                default="NULL",
                output=False,
                optional=True,
                desc=moving_image_masks_desc,
            ),
        )

        self.add_trait(
            "number_of_iterations",
            List(
                List(Int()),
                default=Undefined,
                output=False,
                optional=True,
                desc=number_of_iterations_desc,
            ),
        )

        self.add_trait(
            "radius_or_number_of_bins",
            List(
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=radius_or_number_of_bins_desc,
            ),
        )

        self.add_trait(
            "random_seed",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=random_seed_desc,
            ),
        )

        self.add_trait(
            "sampling_percentage",
            List(
                Float(min=0.0, max=1.0),
                default=Undefined,
                output=False,
                optional=True,
                desc=sampling_percentage_desc,
            ),
        )

        self.add_trait(
            "sampling_strategy",
            List(
                Enum("None", "Regular", "Random"),
                default=Undefined,
                output=False,
                optional=True,
                desc=sampling_strategy_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "composite_transform",
            File(
                Undefined,
                output=True,
                optional=True,
                desc=composite_transform_desc,
            ),
        )

        self.add_trait(
            "inverse_composite_transform",
            File(
                Undefined,
                output=True,
                optional=True,
                desc=inverse_composite_transform_desc,
            ),
        )

        self.add_trait(
            "inverse_warped_image",
            File(
                Undefined,
                output=True,
                optional=True,
                desc=inverse_warped_image_desc,
            ),
        )

        self.add_trait(
            "warped_image",
            File(Undefined, output=True, desc=warped_image_desc),
        )

        self.init_default_traits()

        self.init_process("nipype.interfaces.ants.Registration")

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
        super(Registration, self).list_outputs()

        if (
            self.sampling_percentage == Undefined
            and self.sampling_strategy != Undefined
        ) or (
            self.sampling_percentage != Undefined
            and self.sampling_strategy == Undefined
        ):
            print(
                "\nInitialisation failed. Please, set both (or none) of the "
                "two inputparameters sampling_percentage "
                "and sampling_strategy ...!"
            )
            return

        if (
            self.convergence_threshold != Undefined
            and self.number_of_iterations == Undefined
        ):
            print(
                "\nInitialisation failed. Please, set number_of_iterations "
                "input if convergence_threshold input is set...!"
            )
            return

        if (
            self.convergence_window_size != Undefined
            and self.convergence_threshold == Undefined
        ):
            print(
                "\nInitialisation failed. Please, set convergence_threshold "
                "input if convergence_window_size input is set...!"
            )
            return

        if (
            self.initial_moving_transform_com != Undefined
            and self.initial_moving_transform
        ):
            print(
                "\nInitialisation failed. initial_moving_transform_com and "
                "initial_moving_transform are mutually exclusif...!"
            )
            return

        if (
            self.invert_initial_moving_transform != Undefined
            and not self.initial_moving_transform
        ):
            print(
                "\nInitialisation failed. invert_initial_moving_transform "
                "require initial_moving_transform...!"
            )
            return

        # Outputs definition and tags inheritance (optional)
        if self.moving_image:
            if self.out_prefix == Undefined:
                self.out_prefix = "w_"
                print(
                    "The out_prefix parameter is undefined. Automatically "
                    'set to "w" ...'
                )

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(
                    self.moving_image, EXT
                )

                if not valid_ext:
                    print("\nThe input image format is not recognized ...!")
                    return
                else:
                    self.outputs["warped_image"] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName + "." + in_ext,
                    )

                    self.outputs["composite_transform"] = os.path.join(
                        self.output_directory, fileName + "_Composite.h5"
                    )

                    self.outputs["inverse_composite_transform"] = os.path.join(
                        self.output_directory,
                        fileName + "_InverseComposite.h5",
                    )

                    self.inheritance_dict[
                        self.outputs["warped_image"]
                    ] = self.moving_image

                    self.inheritance_dict[
                        self.outputs["composite_transform"]
                    ] = self.moving_image

                    self.inheritance_dict[
                        self.outputs["inverse_composite_transform"]
                    ] = self.moving_image

                    if self.output_inverse_warped_image:
                        self.outputs["inverse_warped_image"] = os.path.join(
                            self.output_directory,
                            self.out_prefix
                            + "_inverse_"
                            + fileName
                            + "."
                            + in_ext,
                        )
                        self.inheritance_dict[
                            self.outputs["inverse_warped_image"]
                        ] = self.moving_image

            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Registration, self).run_process_mia()

        # TODO: We see in soma_workflow: nipype.interface INFO: file NULL does
        #       not exist, with MRIQC anat

        # Mandatory inputs (in Nipype)
        self.process.fixed_image = self.fixed_image
        self.process.metric = self.metric
        self.process.metric_weight = self.metric_weight
        self.process.moving_image = self.moving_image
        self.process.transforms = self.transforms
        self.process.shrink_factors = self.shrink_factors
        self.process.smoothing_sigmas = self.smoothing_sigmas

        # Optional inputs
        self.process.collapse_output_transforms = (
            self.collapse_output_transforms
        )
        self.process.convergence_threshold = self.convergence_threshold
        self.process.convergence_window_size = self.convergence_window_size
        self.process.dimension = self.dimension
        self.process.fixed_image_masks = self.fixed_image_masks
        self.process.float = self.float
        if self.initial_moving_transform:
            self.process.initial_moving_transform = (
                self.initial_moving_transform
            )
        if self.initial_moving_transform_com:
            self.process.initial_moving_transform_com = (
                self.initial_moving_transform_com
            )
        self.process.initialize_transforms_per_stage = (
            self.initialize_transforms_per_stage
        )
        self.process.interpolation = self.interpolation
        if self.interpolation_parameters:
            self.process.interpolation_parameters = (
                self.interpolation_parameters
            )
        if self.invert_initial_moving_transform:
            self.process.invert_initial_moving_transform = (
                self.invert_initial_moving_transform
            )
        self.process.moving_image_masks = self.moving_image_masks
        self.process.number_of_iterations = self.number_of_iterations
        self.process.radius_bins_item_trait = self.radius_bins_item_trait
        self.process.radius_or_number_of_bins = self.radius_or_number_of_bins
        self.process.random_seed = self.random_seed
        self.process.sampling_percentage = self.sampling_percentage
        self.process.sampling_strategy = self.sampling_strategy
        if self.sigma_units:
            self.process.sigma_units = self.sigma_units
        else:
            self.process.sigma_units = ["vox"] * len(self.metric)
        if self.use_estimate_learning_rate_once:
            self.process.use_estimate_learning_rate_once = (
                self.use_estimate_learning_rate_once
            )
        else:
            self.process.use_estimate_learning_rate_once = [True] * len(
                self.metric
            )
        self.process.transform_parameters = self.transform_parameters
        if self.use_histogram_matching:
            self.process.use_histogram_matching = self.use_histogram_matching
        else:
            self.process.use_histogram_matching = [True] * len(self.metric)
        self.process.winsorize_lower_quantile = self.winsorize_lower_quantile
        self.process.winsorize_upper_quantile = self.winsorize_upper_quantile

        # Ouput flags
        if self.output_inverse_warped_image:
            self.process.output_inverse_warped_image = (
                self.inverse_warped_image
            )
        ifile = os.path.split(self.moving_image)[-1]
        fileName, _ = ifile.rsplit(".", 1)
        self.process.output_transform_prefix = fileName + "_"
        self.process.output_warped_image = self.warped_image
        self.process.write_composite_transform = self.write_composite_transform

        return self.process.run(configuration_dict={})
