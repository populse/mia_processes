# -*- coding: utf-8 -*- #

"""The afni preprocess library of the mia_processes package.

The purpose of this module is to customise the main afni preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - SkullStripping
        - Calc

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# nipype imports
from nipype.interfaces.base import File

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# Other import
import os
from traits.api import Either, Enum, Float, String, Undefined


class Calc(ProcessMIA):
    """
    * Voxel-by-voxel arithmetic on 3D datasets *

    Please, see the complete documentation for the `Calc' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/afni/Calc.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Calc, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_a_desc = ('First input 3D file (a pathlike object or string '
                          'representing a file).')
        in_file_b_desc = ('Second input 3D file (a pathlike object or string '
                          'representing a file).')
        in_file_c_desc = ('Third input 3D file (a pathlike object or string '
                          'representing a file).')
        expr_desc = ('Arithmetic expression to apply between a, b and c '
                     '(a string).')
        output_type_desc = ('Typecodes of the output image formats (one '
                            'of NIFTI, AFNI, NIFTI_GZ).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the output image file '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The calculated files (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait("in_file_a",
                       File(output=False,
                            optional=False,
                            desc=in_file_a_desc))

        self.add_trait("in_file_b",
                       File(value=Undefined,
                            output=False,
                            optional=True,
                            desc=in_file_b_desc))

        self.add_trait("in_file_c",
                       File(value=Undefined,
                            output=False,
                            optional=True,
                            desc=in_file_c_desc))

        self.add_trait("expr",
                       String(Undefined,
                              output=False,
                              optional=False,
                              desc=expr_desc))

        self.add_trait("output_type",
                       Enum('NIFTI',
                            'AFNI',
                            'NIFTI_GZ',
                            output=False,
                            optional=True,
                            desc=output_type_desc))

        self.add_trait("out_prefix",
                       String('c',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))


        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.afni.Calc')

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
        if self.in_file:

            if self.out_prefix == Undefined:
                self.out_prefix = 'c'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "c" ...')

            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:

                    if trail in ['nii', '3D']:

                        if self.output_type == 'NIFTI':
                            trail = 'nii'
                        elif self.output_type == 'AFNI':
                            trail = '3D'
                        elif self.output_type == 'NIFTI_GZ':
                            trail = 'nii.gz'

                        self.outputs['out_file'] = os.path.join(
                            self.output_directory,
                            self.out_prefix + fileName + '.' + trail)

                    else:
                        self.outputs['out_file'] = os.path.join(
                            self.output_directory,
                            self.out_prefix + ifile)
                        print('\nThe input image format does not seem to be '
                              'nii or 3D. This can prevent the process '
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
        super(Calc, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.outputtype = self.output_type
        self.process.out_file = self.out_file

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


class SkullStripping(ProcessMIA):
    """
    * Extract the brain from surrounding tissue from MRI T1-weighted images *

    Please, see the complete documentation for the `SkullStripping' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/afni/SkullStripping.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(SkullStripping, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_desc = ('A 3D-T1 file to skull-strip (a pathlike object or string '
                        'representing a file).')
        output_type_desc = ('Typecodes of the output image formats (one '
                            'of NIFTI, AFNI, NIFTI_GZ).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the skull-stripped image file(s) '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The skull-stripped files (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))


        self.add_trait("output_type",
                       Enum('NIFTI',
                            'AFNI',
                            'NIFTI_GZ',
                            output=False,
                            optional=True,
                            desc=output_type_desc))

        self.add_trait("out_prefix",
                       String('ss_',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))


        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.afni.SkullStrip')

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
        super(SkullStripping, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:

            if self.out_prefix == Undefined:
                self.out_prefix = 'ss'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "ss" ...')

            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:

                    if trail in ['nii', '3D']:

                        if self.output_type == 'NIFTI':
                            trail = 'nii'
                        elif self.output_type == 'AFNI':
                            trail = '3D'
                        elif self.output_type == 'NIFTI_GZ':
                            trail = 'nii.gz'

                        self.outputs['out_file'] = os.path.join(
                            self.output_directory,
                            self.out_prefix + fileName + '.' + trail)

                    else:
                        self.outputs['out_file'] = os.path.join(
                            self.output_directory,
                            self.out_prefix + ifile)
                        print('\nThe input image format does not seem to be '
                              'nii or 3D. This can prevent the process '
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
        super(SkullStripping, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.outputtype = self.output_type
        self.process.out_file = self.out_file

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})
