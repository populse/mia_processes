# -*- coding: utf-8 -*- #

"""The afni preprocess library of the mia_processes package.

The purpose of this module is to customise the main afni preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - SkullStripping

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# nipype imports
from nipype.interfaces.base import (File, InputMultiPath)

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# Other import
import os
from traits.api import Enum, Bool, String, Undefined, Either, Float, List, Int, Tuple
import nibabel as nib


class N4BiasFieldCorrection(ProcessMIA):
    """
    * N4 retrospective bias correction algorithm *

    Please, see the complete documentation for the `N4BiasFieldCorrection' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/afni/N4BiasFieldCorrection.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(N4BiasFieldCorrection, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['ants', 'nipype']

        # Inputs description
        in_file_desc = ('A image file (2D, 3D or 4D) to correct (a pathlike object or string '
                        'representing a file).')
        dimension_desc = ('Image dimension (2, 3 or 4).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the corrected image file '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The corrected file (a pathlike object or a '
                         'string representing a file).')
        bias_image_desc = ('Estimated bias (a pathlike object or a '
                         'string representing a file).')
        negative_values_desc = 'True if negative values are present in input data file, ' \
                               'False otherwise.'

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("dimension",
                       Enum(2,
                            3,
                            4,
                            output=False,
                            optional=True,
                            desc=dimension_desc))

        self.add_trait("out_prefix",
                       String('n4c',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.add_trait("bias_image",
                       File(output=True,
                            desc=bias_image_desc))

        self.add_trait("negative_values",
                       Bool(output=True,
                            desc=negative_values_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.ants.N4BiasFieldCorrection')

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
                self.out_prefix = 'n4c'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "n4c" ...')

            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:

                    if trail in ['nii', 'nii.gz', 'img']:
                        print('\nThe input image format does not seem to be '
                              'nii or img. This can prevent the process '
                              'launch ...!')

                    self.outputs['out_file'] = os.path.join(
                        self.output_directory,
                        self.out_prefix + ifile)

                    self.inheritance_dict[self.outputs[
                        'out_file']] = self.in_file

            else:
                print('No output_directory was found...!\n')
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(N4BiasFieldCorrection, self).run_process_mia()

        input_nii = nib.load(self.in_file)
        datamin = input_nii.get_fdata().min()
        if datamin < 0:
            data = input_nii.get_fdata() - datamin
            newnii = input_nii.__class__(data, input_nii.affine, input_nii.header)

            fileName, trail = self.in_file.rsplit('.', 1)

            nib.save(newnii, fileName + '_scaled.' + trail)
            self.in_file = fileName + '_scaled.' + trail
            self.negative_values = True

            print('\nThe input image contains negative values '
                  'A scaled copy was created and set as input')

        self.process.input_image = self.in_file
        self.process.output_image = self.out_file
        self.process.save_bias = True
        self.process.rescale_intensities = True
        self.process.copy_header = True
        self.process.out_file = self.out_file
        self.process.bias_image = self.bias_image

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


class AffineInitializer(ProcessMIA):
    """
    * A multi-start optimizer for affine registration *

    Please, see the complete documentation for the `AffineInitializer' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/afni/N4BiasFieldCorrection.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(AffineInitializer, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['ants', 'nipype']

        # Inputs description
        fixed_image_desc = ('The fixed reference image (a pathlike object or string '
                        'representing a file).')
        moving_image_desc = 'The moving image to be mapped to the fixed space.'
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the corrected image file '
                           '(a string).')

        # Outputs description
        out_file_desc = ('Output transform file (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait("fixed_image",
                       File(value=Undefined,
                            output=False,
                            optional=True,
                            desc=fixed_image_desc))

        self.add_trait("moving_image",
                       File(output=False,
                            optional=False,
                            desc=moving_image_desc))

        self.add_trait("out_prefix",
                       String('AffineTransform_',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.ants.AffineInitializer')

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
                self.out_prefix = 'AffineTransform_'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "AffineTransform_" ...')

            if self.output_directory:
                ifile = os.path.split(self.moving_image)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:

                    if trail in ['nii', 'nii.gz', 'img']:
                        print('\nThe input image format does not seem to be '
                              'nii or img. This can prevent the process '
                              'launch ...!')

                    self.outputs['out_file'] = os.path.join(
                        self.output_directory,
                        self.out_prefix + ifile)

                    self.inheritance_dict[self.outputs[
                        'out_file']] = self.moving_image

            else:
                print('No output_directory was found...!\n')
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(AffineInitializer, self).run_process_mia()

        self.process.fixed_image = self.fixed_image
        self.process.moving_image = self.moving_image
        self.process.out_file = self.out_file

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


class Registration(ProcessMIA):
    """
    * Registers a moving image to a fixed image using a predefined
     (sequence of) cost function(s) and transformation operations.
     Uses a sequence of three transforms: Rigid, Affine and SyN*

    Please, see the complete documentation for the `Registration' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/afni/N4BiasFieldCorrection.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Registration, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['ants', 'nipype']

        # Inputs description
        fixed_image_desc = ('Image to which the moving image should be '
                            'transformed.')

        moving_image_desc = 'Image that will be registered to the space ' \
                            'of the fixed image.'

        initial_moving_transform_desc = ('A transform or a list of transform '
                                         'that should be applied before the '
                                         'registration begin. (a list of items '
                                         'which are a pathlike object or string '
                                         'representing an existing file).')

        metric_desc = 'The metric to use for each stage (a list of strings which are ' \
                      ' ''CC'' or ''MeanSquares'' or ''Demons'' or ''GC'' or ''MI'' ' \
                      'or Mattes)'

        shrink_factor_desc = '(A list of integers)'

        smoothing_sigmas_desc = '(A list of floats)'

        number_of_iterations_desc = '(A list of integers)'

        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filename of the warped image file '
                           '(a string).')

        output_transform_prefix_desc = ('Specify the string to be prepended to the '
                                        'filename of the transform image file '
                                        '(a string).')

        radius_or_number_of_bins_desc = 'The number of bins in each stage for the MI ' \
                                        'and Mattes metric, the radius for other metrics ' \
                                        '(a list of integers)'

        sampling_percentage_desc = 'The metric sampling percentages to use for each stage ' \
                                   '(a list of 0.0 <= floats <= 1; requires sampling strategy)'

        sampling_strategy_desc = 'The metric sampling strategies for each stage (A list of ' \
                                 'strings which are ''None'' or ''Regular'' or ''Random'').'

        transform_parameters_desc = '(A list of tuples)'

        # Outputs description
        composite_transform_desc = ('Output composite transform file (a pathlike '
                                    'object or a string representing a file).')

        inverse_composite_transform_desc = ('Output inverse composite transform'
                                            ' file (a pathlike object or a '
                                            'string representing a file).')

        warped_image_desc = ('Outputs the warped image. (a pathlike object '
                             'or string representing an existing file).')

        # Inputs traits
        self.add_trait("fixed_image",
                       File(value=Undefined,
                            output=False,
                            optional=True,
                            desc=fixed_image_desc))

        self.add_trait("moving_image",
                       File(output=False,
                            optional=False,
                            desc=moving_image_desc))

        self.add_trait("initial_moving_transform",
                       InputMultiPath(File(),
                                      default=Undefined,
                                      output=False,
                                      optional=True,
                                      desc=initial_moving_transform_desc))

        self.add_trait("metric",
                       List(Enum('CC',
                                 'MeanSquares',
                                 'Demon',
                                 'GC',
                                 'MI',
                                 'Mattes'),
                            minlen=3, maxlen=3,
                            output=False,
                            optional=False,
                            desc=metric_desc))

        self.add_trait("shrink_factor",
                       List(Int(),
                            minlen=3, maxlen=3,
                            output=False,
                            optional=False,
                            desc=shrink_factor_desc))

        self.add_trait("smoothing_sigmas",
                       List(Float(),
                            minlen=3, maxlen=3,
                            output=False,
                            optional=False,
                            desc=smoothing_sigmas_desc))

        self.add_trait("number_of_iterations",
                       List(List(Int()),
                            minlen=3, maxlen=3,
                            output=False,
                            optional=False,
                            desc=number_of_iterations_desc))

        self.add_trait("out_prefix",
                       String('w_',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        self.add_trait("output_transform_prefix",
                       String('transform_',
                              output=False,
                              optional=True,
                              desc=output_transform_prefix_desc))

        self.add_trait("radius_or_number_of_bins",
                       List(Int(),
                            minlen=3, maxlen=3,
                            output=False,
                            optional=False,
                            desc=radius_or_number_of_bins_desc))

        self.add_trait("sampling_percentage",
                       List(Float(min=0.0, max=1.0),
                            default=Undefined,
                            minlen=3, maxlen=3,
                            output=False,
                            optional=True,
                            desc=sampling_percentage_desc))

        self.add_trait("sampling_strategy",
                       List(Enum('None',
                                 'Regular',
                                 'Random'),
                            default=Undefined,
                            minlen=3, maxlen=3,
                            output=False,
                            optional=True,
                            desc=sampling_strategy_desc))

        self.add_trait("transform_parameters",
                       List(Tuple(),
                            minlen=3, maxlen=3,
                            output=False,
                            optional=False,
                            desc=transform_parameters_desc))

        # Outputs traits
        self.add_trait("composite_transform",
                       File(output=True,
                            desc=composite_transform_desc))

        self.add_trait("inverse_composite_transform",
                       File(output=True,
                            desc=inverse_composite_transform_desc))

        self.add_trait("warped_image",
                       File(output=True,
                            desc=warped_image_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.ants.Registration')

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

        if (self.sampling_percentage == Undefined and self.sampling_strategy != Undefined) or \
           (self.sampling_percentage != Undefined and self.sampling_strategy == Undefined):
            print('\nInitialisation failed. Please, set both (or none) of the two input '
                  'parameters sampling_percentage and sampling_strategy ...!')
            return

        # Outputs definition and tags inheritance (optional)
        if self.moving_image:

            if self.out_prefix == Undefined:
                self.out_prefix = 'w_'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "AffineTransform_" ...')

            if self.output_directory:
                ifile = os.path.split(self.moving_image)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:

                    if trail in ['nii', 'nii.gz', 'img']:
                        print('\nThe input image format does not seem to be '
                              'nii or img. This can prevent the process '
                              'launch ...!')

                    self.outputs['warped_image'] = os.path.join(
                        self.output_directory,
                        self.out_prefix + ifile)

                    self.outputs['composite_transform'] = os.path.join(
                        self.output_directory,
                        'comp_transf' + ifile)

                    self.outputs['inverse_composite_transform'] = os.path.join(
                        self.output_directory,
                        'inv_comp_transf' + ifile)

                    self.inheritance_dict[self.outputs[
                        'warped_image']] = self.moving_image

                    self.inheritance_dict[self.outputs[
                        'composite_transform']] = self.moving_image

                    self.inheritance_dict[self.outputs[
                        'inverse_composite_transform']] = self.moving_image

            else:
                print('No output_directory was found...!\n')
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Registration, self).run_process_mia()

        self.process.fixed_image = self.fixed_image
        self.process.moving_image = self.moving_image
        self.process.initial_moving_transform = self.initial_moving_transform

        self.process.metric = self.metric
        self.process.shrink_factor = self.shrink_factor
        self.process.smoothing_sigmas = self.smoothing_sigmas
        self.process.number_of_iterations = self.number_of_iterations
        self.process.output_transform_prefix = self.output_transform_prefix
        self.process.radius_or_number_of_bins = self.radius_or_number_of_bins
        self.process.sampling_percentage = self.sampling_percentage
        self.process.sampling_strategy = self.sampling_strategy
        self.process.transform_parameters = self.transform_parameters

        # default inputs
        self.process.transforms = ['Rigid', 'Affine', 'SyN']
        self.process.metric_weight = [1, 1, 1]
        self.process.collapse_output_transforms = True
        self.process.convergence_threshold = [1e-6, 1e-6, 1e-6]
        self.process.convergence_window_size = [20, 20, 10]
        self.process.interpolation = 'LanczosWindowedSinc'
        self.process.dimension = 3
        self.process.output_warped_image = True
        self.process.sigma_units = ['vox', 'vox', 'vox']
        self.process.use_estimate_learning_rate_once = [True, True, True]
        self.process.use_histogram_matching = [True, True, True]
        self.process.winsorize_lower_quantile = 0.005
        self.process.winsorize_upper_quantile = 0.995
        self.process.write_composite_transform = True

        self.process.warped_image = self.warped_image
        self.process.composite_transform = self.composite_transform
        self.process.inverse_composite_transform = self.inverse_composite_transform

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


class ApplyTransforms(ProcessMIA):
    """
    * Registers a moving image to a fixed image using a predefined
     (sequence of) cost function(s) and transformation operations.
     Uses a sequence of three transforms: Rigid, Affine and SyN*

    Please, see the complete documentation for the `Registration' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/afni/N4BiasFieldCorrection.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ApplyTransforms, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['ants', 'nipype']

        # Inputs description
        input_image_desc = ('Image to apply transformation to. (a pathlike object '
                            'or string representing an existing file)')

        reference_image_desc = 'Reference image space that you wish to warp into ' \
                               '(a pathlike object or string representing an existing file)'

        transforms_desc = 'Transform files that will be applied in reverse order.'

        interpolation_desc = 'Choice of interpolator. (‘Linear’ or ‘NearestNeighbor’ ' \
                             'or ‘CosineWindowedSinc’ or ‘WelchWindowedSinc’ or ' \
                             '‘HammingWindowedSinc’ or ‘LanczosWindowedSinc’ or ' \
                             '‘MultiLabel’ or ‘Gaussian’ or ‘BSpline’)'

        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the corrected image file '
                           '(a string).')

        # Outputs description
        output_image_desc = ('Warped image (a pathlike object or string representing '
                        'an existing file).')

        # Inputs traits
        self.add_trait("input_image",
                       File(output=False,
                            optional=False,
                            desc=input_image_desc))

        self.add_trait("reference_image",
                       File(output=False,
                            optional=False,
                            desc=reference_image_desc))

        self.add_trait("transforms",
                       Either(InputMultiPath(File()),
                              'Identity',
                              default='identity',
                              output=False,
                              optional=False,
                              desc=transforms_desc))

        self.add_trait("interpolation",
                       Enum('Linear',
                            'NearestNeighbor',
                            'CosineWindowedSinc',
                            'WelchWindowedSinc',
                            'HammingWindowedSinc',
                            'LanczosWindowedSinc',
                            'MultiLabel',
                            'Gaussian'
                            'BSpline',
                            default='Linear',
                            output=False,
                            optional=True,
                            desc=interpolation_desc))

        self.add_trait("out_prefix",
                       String('t_',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        # Outputs traits
        self.add_trait("output_image",
                       File(output=True,
                            desc=output_image_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.ants.ApplyTransforms')

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

        if (self.sampling_percentage == Undefined and self.sampling_strategy != Undefined) or \
           (self.sampling_percentage != Undefined and self.sampling_strategy == Undefined):
            print('\nInitialisation failed. Please, set both (or none) of the two input '
                  'parameters sampling_percentage and sampling_strategy ...!')
            return

        # Outputs definition and tags inheritance (optional)
        if self.input_image:

            if self.out_prefix == Undefined:
                self.out_prefix = 't_'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "AffineTransform_" ...')

            if self.output_directory:
                ifile = os.path.split(self.moving_image)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:

                    if trail in ['nii', 'nii.gz', 'img']:
                        print('\nThe input image format does not seem to be '
                              'nii or img. This can prevent the process '
                              'launch ...!')

                    self.outputs['output_image'] = os.path.join(
                        self.output_directory,
                        self.out_prefix + ifile)

                    self.inheritance_dict[self.outputs[
                        'output_image']] = self.input_image

            else:
                print('No output_directory was found...!\n')
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ApplyTransforms, self).run_process_mia()

        self.process.input_image = self.input_image
        self.process.reference_image = self.reference_image
        self.process.transforms = self.transforms

        self.process.dimension = 3
        self.process.default_value = 0.0
        self.process.interpolation = self.interpolation

        self.process.output_image = self.output_image
        self.process.composite_transform = self.composite_transform
        self.process.inverse_composite_transform = self.inverse_composite_transform

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})
