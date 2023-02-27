# -*- coding: utf-8 -*- #

"""The fsl preprocess library of the mia_processes package.

The purpose of this module is to customise the main fsl preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - Segment
        - Smooth
        - SurfacesExtraction

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from nipype.interfaces.base import File

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# Other import
import os
from traits.api import Either, Enum, Float, String, Undefined, List, Bool


class Segment(ProcessMIA):
    """
    * Brain tissue segmentation using fsl.FAST *

    Please, see the complete documention for the `Segment' brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/fsl/Smooth.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Segment, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['fsl', 'nipype']

        # Inputs description
        in_file_desc = ('Files to Segment (a list of items which are existing '
                         'file names).')
        segments_desc = ('Outputs a separate binary image for each tissue type '
                         '(a boolean).')
        output_type_desc = ('Typecodes of the output NIfTI image formats (one '
                            'of NIFTI, NIFTI_PAIR, NIFTI_GZ, NIFTI_PAIR_GZ).')

        # Outputs description
        tissue_class_map_desc = 'path/name of binary segmented volume file '
        partial_volume_files_desc = ('Partial volume files (a list of items '
                                     'which are file names.')

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("segments",
                       Bool(True,
                            output=False,
                            optional=True,
                            desc=segments_desc))

        self.add_trait("output_type",
                       Enum('NIFTI',
                            'NIFTI_PAIR',
                            'NIFTI_GZ',
                            'NIFTI_PAIR_GZ',
                            output=False,
                            optional=True,
                            desc=output_type_desc))

        # Outputs traits
        self.add_trait("tissue_class_map",
                       File(output=True,
                            desc=tissue_class_map_desc))

        self.add_trait("partial_volume_files",
                       List(File(),
                            output=True,
                            desc=partial_volume_files_desc))

        self.init_default_traits()

        # To suppress the "FSLOUTPUTTYPE environment
        # variable is not set" nipype warning:
        if not 'FSLOUTPUTTYPE' in os.environ:
            os.environ['FSLOUTPUTTYPE'] = self.output_type

        self.init_process('nipype.interfaces.fsl.FAST')

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
        super(Segment, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            self.process.output_type = self.output_type

            if self.output_directory:
                _, fileIval = os.path.split(self.in_file)
                self.process.out_basename = os.path.join(self.output_directory, fileIval)

                self.outputs['tissue_class_map'] = \
                    os.path.join(self.output_directory,
                                 os.path.split(
                                     self.process._tissue_class_map)[1])

                self.outputs['partial_volume_files'] = []
                for out_val in self.process._partial_volume_files:
                    self.outputs['partial_volume_files'].append(
                        os.path.join(self.output_directory,
                                     os.path.split(out_val)[1]))

            else:
                print('No output_directory was found...!\n')
                return

        if self.outputs:
            self.inheritance_dict[self.outputs[
                'tissue_class_map']] = self.in_file
            if self.outputs['partial_volume_files'][0]:
                for out_val in self.outputs['partial_volume_files']:
                    self.inheritance_dict[out_val] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Segment, self).run_process_mia()
        _, fileIval = os.path.split(self.in_file)
        self.process.out_basename = os.path.join(self.output_directory,
                                                 fileIval)
        self.process.output_type = self.output_type
        self.process.in_files = self.in_file
        self.process.segments = self.segments

        return self.process.run(configuration_dict={})


class Smooth(ProcessMIA):
    """
    *3D Gaussian smoothing of image volumes*

    Please, see the complete documention for the `Smooth' brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/fsl/Smooth.html>`_

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
        self.requirement = ['fsl', 'nipype']

        # Inputs description
        in_file_desc = ('A file to smooth (a pathlike object or string '
                        'representing a file).')
        fwhm_desc = ('Gaussian kernel fwhm in mm (a float). Mutually exclusive '
                     'with sigma. Basically, 2.3548 * sigma = fwhm.')
        sigma_desc = ('Gaussian kernel sigma in mm (a float). Mutually '
                      'exclusive with fwhm. Basically, 2.3548 * sigma = fwhm.')
        output_type_desc = ('Typecodes of the output NIfTI image formats (one '
                            'of NIFTI, NIFTI_PAIR, NIFTI_GZ, NIFTI_PAIR_GZ).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the smoothed image file(s) '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The smoothed files (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits 
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("fwhm",
                       Either(Undefined,
                              Float(),
                              default=6.0,
                              output=False,
                              optional=True,
                              desc=fwhm_desc))

        self.add_trait("sigma",
                       Either(Undefined,
                              Float(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=sigma_desc))

        self.add_trait("output_type",
                       Enum('NIFTI',
                            'NIFTI_PAIR',
                            'NIFTI_GZ',
                            'NIFTI_PAIR_GZ',
                            output=False,
                            optional=True,
                            desc=output_type_desc))

        self.add_trait("out_prefix",
                       String('s_',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        self.add_trait("outf",
                       String(output=False,
                              optional=True,
                              userlevel=1))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()

        # To suppress the "FSLOUTPUTTYPE environment
        # variable is not set" nipype warning:
        if 'FSLOUTPUTTYPE' not in os.environ:
            os.environ['FSLOUTPUTTYPE'] = self.output_type

        self.init_process('nipype.interfaces.fsl.Smooth')

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
            print('\nInitialisation failed. Please, set one of the two input '
                  'parameters sigma or fwhm ...!')
            return

        elif self.sigma != Undefined and self.fwhm != Undefined:
            print('\nInitialisation failed. Both input parameters "sigma" and '
                  '"fwhm" are mutually exclusive. Please, define only one of '
                  'these two parameters (set the other as Undefined) ...!')
            return

        if self.in_file:

            if self.out_prefix == Undefined:
                self.out_prefix = 's_'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "s" ...')

            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)

                except ValueError:
                    print('\nThe input image format is not recognised ...!')
                    return

                else:

                    if trail in ['nii', 'img']:

                        if self.output_type == 'NIFTI':
                            trail = 'nii'
                        elif self.output_type == 'NIFTI_PAIR':
                            trail = 'img'
                        elif self.output_type == 'NIFTI_GZ':
                            trail = 'nii.gz'
                        elif self.output_type == 'NIFTI_PAIR_GZ':
                            trail = 'img.gz'

                        self.outputs['out_file'] = os.path.join(
                            self.output_directory,
                            self.out_prefix + fileName + '.' + trail)
                        self.outf = os.path.join(self.output_directory,
                                                 self.out_prefix + fileName)

                    else:
                        self.outputs['out_file'] = os.path.join(
                            self.output_directory,
                            self.out_prefix + ifile)
                        print('\nThe input image format does not seem to be '
                              'nii or img. This can prevent the process '
                              'launch ...!')

                    self.inheritance_dict[self.outputs[
                        'out_file']] = self.in_file

            else:
                print('No output_directory was found...!\n')
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Smooth, self).run_process_mia()
        self.process.in_file = self.in_file

        if self.fwhm == Undefined:
            self.process.sigma = self.sigma
            self.process.trait('fwhm').optional = True

        else:
            self.process.fwhm = self.fwhm
            self.process.trait('sigma').optional = True

        self.process.output_type = self.output_type
        self.process.smoothed_file = self.outf
        return self.process.run(configuration_dict={})


class SurfacesExtraction(ProcessMIA):
    """
    * Surfaces (inskull, outskull, outskin) extraction using
    fsl.BET with option -A (bet2 and betsurf) and -n
    (no default brain image output)*

    Runs both bet2 and betsurf programs in order to get skull and scalp
    surfaces created by betsurf.
    This involves registering to standard space in order to allow betsurf
    to find the standard space masks it needs.

    Please, see the complete documention for the `SurfacesExtraction' brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/fsl/SurfacesExtraction.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(SurfacesExtraction, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['fsl', 'nipype']

        # Inputs description
        in_file_desc = ('File to delete non-brain tissue (a pathlike object'
                        'string representing a file)')

        # Outputs description
        outskin_mask_file_desc = 'outskin mask file'
        outskin_mesh_file_desc = 'outskin mesh file'
        outskull_mask_file_desc = 'outskull mask file'
        outskull_mesh_file_desc = 'outskull mesh file'
        inskull_mask_file_desc = 'inskull mask file'
        inskull_mesh_file_desc = 'inskull mesh file'
        skull_mask_file_desc = 'skull mask file'

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        # Outputs traits
        self.add_trait("outskin_mask_file",
                       File(output=True,
                            optional=True,
                            desc=outskin_mask_file_desc))

        self.add_trait("outskin_mesh_file",
                       File(output=True,
                            optional=True,
                            desc=outskin_mesh_file_desc))

        self.add_trait("outskull_mask_file",
                       File(output=True,
                            optional=True,
                            desc=outskull_mask_file_desc))

        self.add_trait("outskull_mesh_file",
                       File(output=True,
                            optional=True,
                            desc=outskull_mesh_file_desc))

        self.add_trait("inskull_mask_file",
                       File(output=True,
                            optional=True,
                            desc=inskull_mask_file_desc))

        self.add_trait("inskull_mesh_file",
                       File(output=True,
                            optional=True,
                            desc=inskull_mesh_file_desc))

        self.add_trait("skull_mask_file",
                       File(output=True,
                            optional=True,
                            desc=skull_mask_file_desc))

        self.init_default_traits()

        # To suppress the "FSLOUTPUTTYPE environment
        # variable is not set" nipype warning:
        if 'FSLOUTPUTTYPE' not in os.environ:
            os.environ['FSLOUTPUTTYPE'] = self.output_type

        self.init_process('nipype.interfaces.fsl.BET')

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
        super(SurfacesExtraction, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            self.process.in_file = self.in_file
            self.process.surfaces = True

            if self.output_directory:
                self.outputs['outskin_mask_file'] = os.path.join(
                    self.output_directory, os.path.split(
                        self.process._outskin_mask_file)[1])
                self.outputs['outskin_mesh_file'] = os.path.join(
                    self.output_directory, os.path.split(
                        self.process._outskin_mesh_file)[1])
                self.outputs['outskull_mask_file'] = os.path.join(
                    self.output_directory, os.path.split(
                        self.process._outskull_mask_file)[1])
                self.outputs['outskull_mesh_file'] = os.path.join(
                    self.output_directory, os.path.split(
                        self.process._outskull_mesh_file)[1])
                self.outputs['inskull_mask_file'] = os.path.join(
                    self.output_directory, os.path.split(
                        self.process._inskull_mask_file)[1])
                self.outputs['inskull_mesh_file'] = os.path.join(
                    self.output_directory, os.path.split(
                        self.process._inskull_mesh_file)[1])
                self.outputs['skull_mask_file'] = os.path.join(
                    self.output_directory, os.path.split(
                        self.process._skull_mask_file)[1])
        if self.outputs:
            self.inheritance_dict[self.outputs[
                'outskin_mask_file']] = self.in_file
            self.inheritance_dict[self.outputs[
                'outskin_mesh_file']] = self.in_file
            self.inheritance_dict[self.outputs[
                'outskull_mask_file']] = self.in_file
            self.inheritance_dict[self.outputs[
                'outskull_mesh_file']] = self.in_file
            self.inheritance_dict[self.outputs[
                'inskull_mask_file']] = self.in_file
            self.inheritance_dict[self.outputs[
                'inskull_mesh_file']] = self.in_file
            self.inheritance_dict[self.outputs[
                'skull_mask_file']] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(SurfacesExtraction, self).run_process_mia()
        self.process.in_file = self.in_file

        # default inputs
        self.process.surfaces = True
        self.process.no_output = True

        return self.process.run(configuration_dict={})
