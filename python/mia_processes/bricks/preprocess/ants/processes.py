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
from nipype.interfaces.base import File

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# Other import
import os
from traits.api import Enum, Bool, String, Undefined
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

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})
