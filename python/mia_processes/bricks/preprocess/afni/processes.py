# -*- coding: utf-8 -*- #

"""The afni preprocess library of the mia_processes package.

The purpose of this module is to customise the main afni preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - Automask
        - Calc
        - Despike
        - DropTRs
        - Mean
        - RefitDeoblique
        - SkullStripping
        - TShift
        - Volreg

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
from mia_processes.utils import checkFileExt

# Other import
import os
from traits.api import Bool, Either, Enum, Float, Int, List, String, Undefined
import nibabel as nib

EXT = {'NIFTI_GZ': 'nii.gz',
       'NIFTI': 'nii'
       }


class Automask(ProcessMIA):
    """
    * Create a brain-only mask of the image using AFNI 3dAutomask command *

    Please, see the complete documentation for the `Automask' brick in
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
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_desc = ('Input file (a pathlike object or string '
                        'representing a file).')
        output_type_desc = ('Typecodes of the output image formats (one '
                            'of NIFTI, NIFTI_GZ).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the output image file '
                           '(a string).')
        out_brain_suffix_desc = ('Suffix of the brain masked image (a string)')
        clfrac_desc = ('Sets the clip level fraction (must be 0.1-0.9).'
                       'A small value will tend to make the mask larger '
                       '(a float)')
        dilate_desc = ('Dilate the mask outwards (an integer)')
        erode_desc = ('Erode the mask inwards (an integer)')

        # Outputs description
        out_file_desc = ('The brain mask file (a pathlike object or a '
                         'string representing a file).')

        brain_file_desc = ('The masked brain file (a pathlike object or a '
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
                       String('automask_',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        self.add_trait("out_brain_suffix",
                       String('_masked',
                              output=False,
                              optional=True,
                              desc=out_brain_suffix_desc))

        self.add_trait("clfrac",
                       Float(0.5,
                             output=False,
                             optional=True,
                             desc=clfrac_desc))

        self.add_trait("erode",
                       Either(Undefined,
                              Int(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=erode_desc))

        self.add_trait("dilate",
                       Either(Undefined,
                              Int(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=dilate_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.add_trait("brain_file",
                       File(output=True,
                            optional=True,
                            desc=brain_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.afni.Automask')

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
                self.out_prefix = 'automask_'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "automask" ...')

            if self.out_brain_suffix == Undefined:
                self.out_brain_suffix = '_masked'
                print('The out_brain_suffix parameter is undefined.'
                      'Automatically set to "_masked" ...')

            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print('\nThe input image format is not recognized ...!')
                    return
                else:
                    self.outputs['out_file'] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName + '.' + EXT[
                            self.output_type])

                    self.outputs['brain_file'] = os.path.join(
                        self.output_directory,
                        fileName + self.out_brain_suffix + '.' + EXT[
                            self.output_type])

                    self.inheritance_dict[self.outputs[
                        'out_file']] = self.in_file

                    self.inheritance_dict[self.outputs[
                        'brain_file']] = self.in_file

            else:
                print('No output_directory was found...!\n')
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
    * Voxel-by-voxel arithmetic on 3D datasets *

    Please, see the complete documentation for the `Calc' brick
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
                            'of NIFTI, NIFTI_GZ).')
        single_idx_desc = ('Volume index for in_file_a.'
                           '(an integer or Undefined)')
        start_idx_desc = ('Start index for in_file_a (an integer'
                          'or Undefined). Requires inputs: stop_idx')
        stop_idx_desc = ('Stop index for in_file_a (an integer or Undefined).'
                         'Requires inputs: start_idx.')
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
                              optional=True,
                              desc=expr_desc))
        self.expr = 'a*step(b)'

        self.add_trait("output_type",
                       Enum('NIFTI',
                            'NIFTI_GZ',
                            output=False,
                            optional=True,
                            desc=output_type_desc))

        self.add_trait("single_idx",
                       Either(Undefined,
                              Int(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=single_idx_desc))

        self.add_trait("start_idx",
                       Either(Undefined,
                              Int(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=start_idx_desc))
        self.add_trait("stop_idx",
                       Either(Undefined,
                              Int(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=stop_idx_desc))

        self.add_trait("out_prefix",
                       String('c_',
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
        if (self.start_idx != Undefined and
                self.stop_idx == Undefined) or (self.start_idx == Undefined and
                                                self.stop_idx != Undefined):
            print('\nInitialisation failed. "start_idx" parameter required'
                  '"stop_idx" parameters and vice versa')
            return
        if self.in_file_a:

            if self.out_prefix == Undefined:
                self.out_prefix = 'c_'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "c" ...')

            if self.output_directory:
                valid_ext, in_ext, fileName_a = checkFileExt(self.in_file_a,
                                                             EXT)

                if not valid_ext:
                    print('\nThe input image format is not recognized ...!')
                    return
                else:
                    self.outputs['out_file'] = os.path.join(
                        self.output_directory,
                        self.out_prefix + fileName_a + '.' + EXT[
                            self.output_type])

                    self.inheritance_dict[self.outputs[
                        'out_file']] = self.in_file_a

            else:
                print('No output_directory was found...!\n')
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


class Despike(ProcessMIA):
    """
    * Removes ‘spikes’ from the 3D+time input dataset *

    Please, see the complete documentation for the `Despike' brick
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
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_desc = ('Input 3D file (a pathlike object or string '
                        'representing a file).')
        despike_desc = 'Despike dataset only if true (boolean).'
        output_type_desc = ('Typecodes of the output image formats (one '
                            'of NIFTI, NIFTI_GZ).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the output image file '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The despiked file (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("despike",
                       Bool(True,
                            output=False,
                            optional=True,
                            desc=despike_desc))

        self.add_trait("output_type",
                       Enum('NIFTI',
                            'NIFTI_GZ',
                            output=False,
                            optional=True,
                            desc=output_type_desc))

        self.add_trait("out_prefix",
                       String('d_',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.afni.Despike')

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
                    self.out_prefix = 'd_'
                    print('The out_prefix parameter is undefined.'
                          'Automatically set to "d" ...')

                if self.output_directory:
                    valid_ext, in_ext, fileName = checkFileExt(self.in_file,
                                                               EXT)

                    if not valid_ext:
                        print('\nThe input image format is'
                              ' not recognized...!')
                        return
                    else:
                        self.outputs['out_file'] = os.path.join(
                            self.output_directory,
                            self.out_prefix + fileName + '.' + EXT[
                                self.output_type])

                else:
                    print('No output_directory was found...!\n')
                    return
            else:
                self.outputs['out_file'] = self.in_file

            self.inheritance_dict[self.outputs[
                'out_file']] = self.in_file

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


class DropTRs(ProcessMIA):
    """
    * DropTRs of bold datasets *

    Please, see the complete documentation for the `DropTRs' brick in
    the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/DropTRs.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(DropTRs, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_desc = ('First input 3D file (a pathlike object or string '
                        'representing a file).')
        start_idx_desc = 'start index (inclusive) for in_file (an Int).'
        stop_idx_desc = 'stop index (exclusive) for in_file (an Int).'
        output_type_desc = ('Typecodes of the output image formats (one '
                            'of NIFTI, AFNI, NIFTI_GZ).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the output image file '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The TR cropped file (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("start_idx",
                       Int(0,
                           output=False,
                           optional=True,
                           desc=start_idx_desc))

        self.add_trait("stop_idx",
                       Int(-1,
                           output=False,
                           optional=True,
                           desc=stop_idx_desc))

        self.add_trait("output_type",
                       Enum('NIFTI',
                            'AFNI',
                            'NIFTI_GZ',
                            output=False,
                            optional=True,
                            desc=output_type_desc))

        self.add_trait("out_prefix",
                       String('cropped_',
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
        super(DropTRs, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:

            if not self.stop_idx or self.stop_idx == -1:
                print("\nWarning: stop_idx will be automatically set to"
                      "the length of input file")
            elif self.stop_idx <= self.start_idx:
                print("\nError: stop_idx cannot be lower than or equal to"
                      "start_idx")
                return

            if self.start_idx == 0 and self.stop_idx == -1:
                self.outputs['out_file'] = self.in_file
            else:
                if self.out_prefix == Undefined:
                    self.out_prefix = 'cropped_'
                    print('The out_prefix parameter is undefined.'
                          'Automatically set to "cropped" ...')

                if self.output_directory:
                    valid_ext, in_ext, fileName = checkFileExt(self.in_file,
                                                               EXT)

                    if not valid_ext:
                        print('\nThe input image format is'
                              ' not recognized...!')
                        return
                    else:
                        self.outputs['out_file'] = os.path.join(
                            self.output_directory,
                            self.out_prefix + fileName + '.' + EXT[
                                self.output_type])
                else:
                    print('No output_directory was found...!\n')
                    return

            self.inheritance_dict[self.outputs[
                'out_file']] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(DropTRs, self).run_process_mia()

        if self.in_file == self.out_file:
            return

        try:
            imnii = nib.load(self.in_file)
            nb_volumes = imnii.dataobj.shape[3]

        except (nib.filebasedimages.ImageFileError,
                FileNotFoundError, TypeError) as e:
            print("\nError while opening input file"
                  ": ", e)
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


class Mean(ProcessMIA):
    """
    * Mean of bold images *

    Please, see the complete documentation for the `Mean' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/Mean.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Mean, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_desc = ('A bold file to be averaged (a pathlike object or string '
                        'representing a file).')
        output_type_desc = ('Typecodes of the output image formats (one '
                            'of NIFTI, AFNI, NIFTI_GZ).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the skull-stripped image file(s) '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The time shifted file (a pathlike object or a '
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
                       String('mean_',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.afni.TStat')

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
        super(Mean, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:

            if self.out_prefix == Undefined:
                self.out_prefix = 'mean_'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "mean" ...')

            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)
                    if trail == 'gz':
                        (fileName_2,
                         trail_2) = os.path.splitext(fileName)
                        if trail_2 == 'nii':
                            trail = 'nii.gz'

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:

                    if trail in ['nii', '3D', 'nii.gz']:

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
        super(Mean, self).run_process_mia()

        self.process.in_file = self.in_file
        self.process.args = '-mean'
        self.process.outputtype = self.output_type
        self.process.out_file = self.out_file

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


class RefitDeoblique(ProcessMIA):
    """
    * Deoblique dataset *

    Please, see the complete documentation for the `Deoblique' brick in
    the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/Deoblique.html>`_

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
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_desc = ('First input 3D file (a pathlike object or string '
                        'representing a file).')
        deoblique_desc = 'Deoblique dataset only if true (boolean).'

        # Outputs description
        out_file_desc = ('The deobliqued file (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("deoblique",
                       Bool(True,
                            output=False,
                            optional=True,
                            desc=deoblique_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.afni.Refit')

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

            self.outputs['out_file'] = self.in_file

            self.inheritance_dict[self.outputs[
                'out_file']] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(RefitDeoblique, self).run_process_mia()

        self.process.in_file = self.in_file
        self.process.deoblique = self.deoblique

        return self.process.run(configuration_dict={})


class SkullStripping(ProcessMIA):
    """
    * Extract the brain from surrounding tissue from MRI T1-weighted images *

    Please, see the complete documentation for the `SkullStripping' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/afni/SkullStripping.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

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
                self.out_prefix = 'ss_'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "ss" ...')

            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)
                    if trail == 'gz':
                        (fileName_2,
                         trail_2) = os.path.splitext(fileName)
                        if trail_2 == 'nii':
                            trail = 'nii.gz'

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:

                    if trail in ['nii', '3D', 'nii.gz']:

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


class TShift(ProcessMIA):
    """
    * Slice-time correction of bold images *

    Please, see the complete documentation for the `TShift' brick in the populse.mia_processes website
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
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_desc = ('A bold file to be time-shifted (a pathlike object or string '
                        'representing a file).')
        slice_encoding_dir_desc = ('Direction in which slice_timing is specified '
                                   '(default: k). If negative,slice_timing is defined '
                                   'in reverse order, that is, the first entry '
                                   'corresponds to the slice with the largest index, '
                                   'and the final entry corresponds to slice index zero.')
        slice_timing_desc = ('Time offsets from the volume acquisition onset for each '
                             'slice. (a list of floats).')
        output_type_desc = ('Typecodes of the output image formats (one '
                            'of NIFTI, AFNI, NIFTI_GZ).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the skull-stripped image file(s) '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The time shifted file (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("slice_encoding_dir",
                       Enum('k',
                            'k-',
                            default='k',
                            output=False,
                            optional=True,
                            desc=slice_encoding_dir_desc))

        self.add_trait("slice_timing",
                       List(Float(),
                            output=False,
                            optional=True,
                            desc=slice_timing_desc))

        self.add_trait("output_type",
                       Enum('NIFTI',
                            'AFNI',
                            'NIFTI_GZ',
                            output=False,
                            optional=True,
                            desc=output_type_desc))

        self.add_trait("out_prefix",
                       String('st_corr_',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.afni.TShift')

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

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.slice_timing:
                if self.out_prefix == Undefined:
                    self.out_prefix = 'st_corr_'
                    print('The out_prefix parameter is undefined. Automatically '
                          'set to "st_corr" ...')

                if self.output_directory:
                    ifile = os.path.split(self.in_file)[-1]

                    try:
                        fileName, trail = ifile.rsplit('.', 1)
                        if trail == 'gz':
                            (fileName_2,
                             trail_2) = os.path.splitext(fileName)
                            if trail_2 == 'nii':
                                trail = 'nii.gz'

                    except ValueError:
                        print('\nThe input image format is not recognized ...!')
                        return

                    else:

                        if trail in ['nii', '3D', 'nii.gz']:

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

                else:
                    print('No output_directory was found...!\n')
                    return
            else:
                self.outputs['out_file'] = self.in_file

            self.inheritance_dict[self.outputs[
                'out_file']] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(TShift, self).run_process_mia()
        if not self.slice_timing:
            return

        self.process.in_file = self.in_file
        self.process.outputtype = self.output_type
        self.process.slice_encoding_direction = self.slice_encoding_direction
        self.process.slice_timing = self.slice_timing
        self.process.out_file = self.out_file

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


class Volreg(ProcessMIA):
    """
    * Register input volumes to a base volume using AFNI 3dvolreg command *

    Please, see the complete documentation for the `Volreg' brick in the populse.mia_processes website
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
        self.requirement = ['afni', 'nipype']

        # Inputs description
        in_file_desc = ('A bold file to be time-shifted (a pathlike object or string '
                        'representing a file).')
        interp_desc = (' Spatial interpolation methods (Either ‘Fourier’ or ‘cubic’ '
                       'or ‘heptic’ or ‘quintic’ or ‘linear’ - default = heptic)')
        two_pass_desc = ('Do two passes of the registration algorithm:'
                         '(1) with smoothed base and data bricks, with linear'
                         ' interpolation, to get a crude alignment, then'
                         '(2) with the input base and data bricks, to get a fine '
                         'alignment.'
                         'This method is useful when aligning high-resolution datasets '
                         'that may need to be moved more than a few voxels to be aligned.'
                         '(a boolean).')
        zpad_desc = ('Zeropad around the edges by ‘n’ voxels during rotations '
                     '(an integer).')
        output_type_desc = ('Typecodes of the output image formats (one '
                            'of NIFTI, AFNI, NIFTI_GZ).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the skull-stripped image file(s) '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The registered file (a pathlike object or a '
                         'string representing a file).')
        oned_file_desc = ('The movement parameters file (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("interp",
                       Enum('Fourier',
                            'cubic',
                            'heptic',
                            'quintic',
                            'linear',
                            default='heptic',
                            output=False,
                            optional=True,
                            desc=interp_desc))

        self.add_trait("twopass",
                       Bool(False,
                            output=False,
                            optional=True,
                            desc=two_pass_desc))

        self.add_trait("zpad",
                       Int(4,
                           output=False,
                           optional=True,
                           desc=zpad_desc))

        self.add_trait("output_type",
                       Enum('NIFTI',
                            'AFNI',
                            'NIFTI_GZ',
                            output=False,
                            optional=True,
                            desc=output_type_desc))

        self.add_trait("out_prefix",
                       String('reg_',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.add_trait("oned_file",
                       File(output=True,
                            desc=oned_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.afni.Volreg')

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
                self.out_prefix = 'reg_'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "reg" ...')

            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)
                    if trail == 'gz':
                        (fileName_2,
                         trail_2) = os.path.splitext(fileName)
                        if trail_2 == 'nii':
                            trail = 'nii.gz'

                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                else:

                    if trail in ['nii', '3D', 'nii.gz']:

                        if self.output_type == 'NIFTI':
                            trail = 'nii'
                        elif self.output_type == 'AFNI':
                            trail = '3D'
                        elif self.output_type == 'NIFTI_GZ':
                            trail = 'nii.gz'

                        self.outputs['out_file'] = os.path.join(
                            self.output_directory,
                            self.out_prefix + fileName + '.' + trail)

                        self.outputs['oned_file'] = os.path.join(
                            self.output_directory,
                            self.out_prefix + fileName + '_oned.txt')

                    else:
                        self.outputs['out_file'] = os.path.join(
                            self.output_directory,
                            self.out_prefix + ifile)

                        self.outputs['oned_file'] = os.path.join(
                            self.output_directory,
                            self.out_prefix + fileName + '_oned.txt')
                        print('\nThe input image format does not seem to be '
                              'nii or 3D. This can prevent the process '
                              'launch ...!')

                    self.inheritance_dict[self.outputs[
                        'out_file']] = self.in_file

                    self.inheritance_dict[self.outputs[
                        'oned_file']] = self.in_file

            else:
                print('No output_directory was found...!\n')
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Volreg, self).run_process_mia()
        self.process.in_file = self.in_file
        self.process.outputtype = self.output_type
        self.process.interp = self.interp
        if self.twopass:
            self.process.args = '-twopass'
        if self.zpad:
            self.process.zpad = self.zpad
        self.process.out_file = self.out_file
        self.process.oned_file = self.oned_file

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})
