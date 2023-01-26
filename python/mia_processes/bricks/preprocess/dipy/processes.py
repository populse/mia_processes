# -*- coding: utf-8 -*- #

"""The dipy preprocess library of the mia_processes package.

The purpose of this module is to customise the main dipy preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - Denoise

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# nibabel import
import nibabel as nib
import nibabel.processing as nibp

# nipype imports
from nipype.interfaces.base import File

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# Other import
import os
from traits.api import Either, Enum, Float, String, Undefined
import numpy as np

class Denoise(ProcessMIA):
    """
    * Non-local means for denoising 3D images *

    Please, see the complete documentation for the `Denoise' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/dipy/Denoise.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Denoise, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['nipype']

        # Inputs description
        in_file_desc = ('A file to denoise (a pathlike object or string '
                        'representing a file).')
        seg_file_desc = ('A segmentation file to calculate SNR (a pathlike '
                         'object or string representing a file). Mutually '
                         'exclusive with snr')
        snr_desc = ('Signal to noise ratio (a float). If undefined, SNR is '
                    'estimated from in_file and seg_file. If seg_file also '
                    'undefined, SNR is estimated from in_file only. Mutually '
                    'exclusive with seg_file')
        in_file_snr_desc = ('A input file to calculate SNR with seg_file. If'
                            'not specified, equal to in_file (a pathlike object'
                            ' or string representing a file).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the smoothed image file(s) '
                           '(a string).')

        # Outputs description
        out_file_desc = ('The denoised file (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait("in_file",
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait("in_file_snr",
                       File(output=False,
                            optional=True,
                            desc=in_file_snr_desc))

        self.add_trait("seg_file",
                       File(output=False,
                            optional=True,
                            desc=seg_file_desc))

        self.add_trait("snr",
                       Either(Undefined,
                              Float(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=snr_desc))

        # Outputs traits
        self.add_trait("out_file",
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()

        self.init_process('nipype.interfaces.dipy.Denoise')

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
        super(Denoise, self).list_outputs()

        if self.seg_file and self.snr is not Undefined:
                print('\nInitialisation failed. Please, do not set both seg_file '
                      'and snr ...!')
                return

        # Outputs definition and tags inheritance (optional)
        if self.in_file:

            if self.output_directory:
                ifile = os.path.split(self.in_file)[-1]

                try:
                    fileName, trail = ifile.rsplit('.', 1)

                except ValueError:
                    print('\nThe input image format is not recognised ...!')
                    return

                else:

                    if trail not in ['nii', 'img']:
                        print('\nThe input image format does not seem to be '
                              'nii or img. This can prevent the process '
                              'launch ...!')

                    self.outputs['out_file'] = os.path.join(
                        self.output_directory,
                        fileName + '_denoise.' + trail)

                    self.inheritance_dict[self.outputs[
                        'out_file']] = self.in_file

                # self.outputs['out_file'] = self.process._out_file
            else:
                print('No output_directory was found...!\n')
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Denoise, self).run_process_mia()
        self.process.in_file = self.in_file

        if self.seg_file:
            if self.in_file_snr:
                file_name = self.in_file_snr
            else:
                file_name = self.in_file

            seg_file_name = self.seg_file

            try:
                img = nib.load(file_name)
                seg_img = nib.load(seg_file_name)
            except (nib.filebasedimages.ImageFileError,
                    FileNotFoundError, TypeError) as e:
                print("\nError with files, during "
                      "initialisation: ", e)
                img = None
                seg_img = None

            if (img is not None) and (seg_img is not None):
                data = img.get_fdata()
                mask = seg_img.get_fdata() == 2  # WM label
                self.process.snr = float(np.mean(data[mask]) /
                                         (data[mask].std() * np.sqrt(mask.sum() / (mask.sum() - 1))))
        else:
            if self.snr is Undefined:
                self.process.snr = None
            else:
                self.process.snr = self.snr

        self.process._out_file = self.out_file

        return self.process.run(configuration_dict={})
