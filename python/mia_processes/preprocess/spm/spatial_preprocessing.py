# -*- coding: utf-8 -*- #

"""The spm preprocess library of the mia_processes package.

The purpose of this module is to customise the main spm preprocessing bricks
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

# mia_processes import
from mia_processes.process_mia import Process_Mia
#from .nipype_extension import NewSegmentMia

# nipype import
from nipype.interfaces import spm
from nipype.interfaces.base import OutputMultiPath, InputMultiPath, File, traits
from nipype.interfaces.spm.base import ImageFileSPM

# Other import
import os
#from traits.api import Float

# populse_mia import
#from populse_mia.software_properties import Config


class Smooth(Process_Mia):
    """Smooth (from mia_processes.preprocess.spm.spatial_preprocessing) <=>
                                                           Smooth (SPM12 names).
*** Smooth: 3D Gaussian smoothing of image volumes ***
    * Inputs parameters:
        # in_files <=> data: List of files to smooth. A list of items which are
                             an existing, uncompressed file (valid extensions:
                             [.img, .nii, .hdr]).
            <ex. ['/home/ArthurBlair/data/raw_data/Func.nii']>
        # fwhm <=> fwhm: Specify  the  full-width at half maximum (FWHM) of the
                         Gaussian smoothing kernel in mm. Three values should be
                         entered, denoting the FWHM in the x, y and z
                         directions. A list of 3 items which are a float of fwhm
                          for each dimension.
            <ex. [6, 6, 6]>
        # data_type <=> dtype: Data type of the output images (an integer 
                               [int or long]).
                               0: same as the original images
                               2: UINT8 (unsigned char)
                               4: INT16 (signed short)
                               6: INT32 (signed int)
                               8: FLOAT32 (single prec. float)
                              10: FLOAT64 (double prec. float)
                              <undefined> = 0 ?
            <ex. 0, MIA_processes default value>
        # implicit_masking <=> im: A mask implied by a particular voxel value
                                   (a boolean). If set to True, the implicit
                                   masking of the input image is preserved in
                                   the smoothed image.
            <ex. False>
        # out_prefix <=> prefix: Specify the string to be prepended to the
                                 filenames of the smoothed image file(s)
                                 (a string).
            <ex. s, capsul/nipype default value>
    * Outputs parameters:
        # smoothed_files: Smoothed files (a list of items which are an existing
                          file name).
            <ex. /home/ArthurBlair/data/raw_data/sFunc.nii>

    """

    
    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Smooth, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['matlab', 'spm']


        # Inputs description
        in_files_desc = ('List of files to smooth. A list of items which are '
                         'an existing, uncompressed file '
                         '(valid extensions: [.img, .nii, .hdr]).')
        fwhm_desc = ('Full-width at half maximum (FWHM) of the Gaussian '
                     'smoothing kernel in mm. A list of 3 items which are a '
                     'float of fwhm for each dimension.')
        data_type_desc = ('Data type of the output images '
                          '(an integer [int or long]).')
        implicit_masking_desc = ('A mask implied by a particular voxel value '
                                 '(a boolean).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the smoothed image file(s) '
                           '(a string).')

        # Outputs description
        smoothed_files_desc = ('The smoothed files (a list of items which are '
                               'an existing file name).')

        # Input traits 
        self.add_trait("in_files",
                       InputMultiPath(ImageFileSPM(),
                                      copyfile=False,
                                      output=False,
                                      desc=in_files_desc))
        self.add_trait("fwhm",
                       traits.Either(traits.Float(),
                                     traits.List(traits.Float(),
                                                 minlen=3, maxlen=3),
                                     default=[6.0,6.0,6.0],
                                     output=False,
                                     optional=True,
                                     desc= fwhm_desc))

        self.add_trait("data_type",
                       traits.Int(output=False,
                                  optional=True,
                                  desc=data_type_desc))
        
        self.add_trait("implicit_masking",
                       traits.Bool(output=False,
                                   optional=True,
                                   desc=implicit_masking_desc))
        
        self.add_trait("out_prefix",
                       traits.String('s',
                                     usedefault=True,
                                     output=False,
                                     optional=True,
                                     desc=out_prefix_desc))

        # Output traits 
        self.add_trait("smoothed_files",
                       OutputMultiPath(File(),
                                       output=True,
                                       desc=smoothed_files_desc))

        self.process = spm.Smooth()
        self.change_dir = True

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return 
        self.make_initResult() object.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Smooth, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_files:
            self.process.inputs.in_files = self.in_files

            if self.out_prefix:
                self.process.inputs.out_prefix = self.out_prefix

            self.outputs['smoothed_files'] = self.process._list_outputs()[
                                                               'smoothed_files']

        """raw_data_folder = os.path.join("data", "raw_data")
        derived_data_folder = os.path.join("data", "derived_data")
        for out_name, out_value in outputs.items():
            if type(out_value) is list:
                for idx, element in enumerate(out_value):
                    # To change the raw_data_folder to the derived_data_folder
                    element = element.replace(raw_data_folder, derived_data_folder)
                    outputs[out_name][idx] = element
                    self.smoothed_files = outputs["smoothed_files"]
        """
        if self.outputs:
        
            for key, values in self.outputs.items():
            
                if key == "smoothed_files":
                
                    for fullname in values:
                        path, filename = os.path.split(fullname)
                    
                        if self.out_prefix:
                            filename_without_prefix = filename[
                                                          len(self.out_prefix):]
                        
                        else:
                            filename_without_prefix = filename[len('s'):]

                        if (os.path.join(path,
                                         filename_without_prefix)
                              in self.in_files):
                            self.inheritance_dict[fullname] = os.path.join(path,
                                                        filename_without_prefix)
        
        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()
    
    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Smooth, self).run_process_mia()
        
        for idx, element in enumerate(self.in_files):
            full_path = os.path.relpath(element)
            self.in_files[idx] = full_path

        self.process.inputs.in_files = self.in_files
        self.process.inputs.fwhm = self.fwhm
        self.process.inputs.data_type = self.data_type
        self.process.inputs.implicit_masking = self.implicit_masking
        self.process.inputs.out_prefix = self.out_prefix
        self.process.run()
