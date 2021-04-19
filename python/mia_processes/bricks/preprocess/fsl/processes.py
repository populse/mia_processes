# -*- coding: utf-8 -*- #

"""The fsl preprocess library of the mia_processes package.

The purpose of this module is to customise the main fsl preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - Smooth

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# Capsul import
#from capsul.api import StudyConfig, get_process_instance

# nibabel import
#import nibabel as nib

# nipype imports
#from nipype.interfaces import fsl
from nipype.interfaces.base import File

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# Other import
import os
from traits.api import Either, Enum, Float, String, Undefined, Union


class Smooth(ProcessMIA):
    """
    *3D Gaussian smoothing of image volumes*

    Please, see the complete documention for the `Smooth brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/fsl/Smooth.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.

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
                       String('s',
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
        if not 'FSLOUTPUTTYPE' in os.environ:
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
                self.out_prefix = 's'
                print('The out_prefix parameter is undefined. Automatically '
                      'set to "s" ...')
            
            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName ,trail = ifile.rsplit('.', 1)

                except ValueError:
                    print('\nThe input image format is not recognised ...!')
                    return

                else:

                    if trail in ['nii', 'img']:

                        if self.output_type == 'NIFTI': trail = 'nii'
                        elif self.output_type == 'NIFTI_PAIR': trail = 'img'
                        elif self.output_type == 'NIFTI_GZ': trail = 'nii.gz'
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
