# mia_processes import
from mia_processes.process_mia import Process_Mia
#from .nipype_extension import NewSegmentMia

# nipype import
from nipype.interfaces import spm
from nipype.interfaces.base import OutputMultiPath, InputMultiPath, File, traits
from nipype.interfaces.spm.base import ImageFileSPM

# Other import
import os
from traits.api import Float

# populse_mia import
from populse_mia.software_properties import Config


                        ### Bricks/classes in this module: ###
                        # - Coregister    Vu                 #
                        # - NewSegment    Vu                 #
                        # - Normalize     Vu                 #
                        # - Realign       Vu                 #
                        # - Smooth        Vu                 #
                        ######################################
                        ###   No functions in this module   ##
                        ######################################


class Smooth(Process_Mia):
    """
- Smooth (User_processes.preprocess.spm.spatial_preprocessing.Smooth) <=> Smooth (SPM12 names).
*** Smooth: 3D Gaussian smoothing of image volumes ***
    * Dependencies: Use_Matlab: True; Use_SPM: True
    * Inputs parameters:
        # in_files <=> data: List of files to smooth. A list of items which are an existing, uncompressed file
                             (valid extensions: [.img, .nii, .hdr]).
            <ex. ['/home/ArthurBlair/data/raw_data/Func.nii']>
        # fwhm <=> fwhm: Specify  the  full-width at half maximum (FWHM) of the Gaussian smoothing kernel in mm.
                         Three values should be entered, denoting the FWHM in the x, y and z directions. A list
                         of 3 items which are a float of fwhm for each dimension.
            <ex. [6, 6, 6]>
        # data_type <=> dtype: Data type of the output images (an integer [int or long]).
                               0: same as the original images
                               2: UINT8 (unsigned char)
                               4: INT16 (signed short)
                               6: INT32 (signed int)
                               8: FLOAT32 (single prec. float)
                              10: FLOAT64 (double prec. float)
                              <undefined> = 0 ?
            <ex. 0, MIA_processes default value>
        # implicit_masking <=> im: A mask implied by a particular voxel value (a boolean). If set to True, the
                                   implicit masking of the input image is preserved in the smoothed image.
            <ex. False>
        # out_prefix <=> prefix: Specify the string to be prepended to the filenames of the smoothed
                                 image file(s) (a string).
            <ex. s, capsul/nipype default value>
    * Outputs parameters:
        # smoothed_files: Smoothed files (a list of items which are an existing file name).
            <ex. /home/ArthurBlair/data/raw_data/sFunc.nii>
    """
    
    def __init__(self):
        super(Smooth, self).__init__()

        # Inputs description
        in_files_desc = 'List of files to smooth. A list of items which are an existing, uncompressed file (valid extensions: [.img, .nii, .hdr]).'
        fwhm_desc = 'Full-width at half maximum (FWHM) of the Gaussian smoothing kernel in mm. A list of 3 items which are a float of fwhm for each dimension.'
        data_type_desc = 'Data type of the output images (an integer [int or long]).'
        implicit_masking_desc = 'A mask implied by a particular voxel value (a boolean).'
        out_prefix_desc = 'Specify  the string to be prepended to the filenames of the smoothed image file(s) (a string).'

        # Outputs description
        smoothed_files_desc = 'The smoothed files (a list of items which are an existing file name).'

        # Input traits 
        self.add_trait("in_files",
                       InputMultiPath(ImageFileSPM(),
                                      copyfile=False,
                                      output=False,
                                      desc=in_files_desc))
        '''self.add_trait("fwhm", traits.Either(traits.Float(),
                       traits.List(traits.Float()), default_value=[6, 6, 6], output=False, optional=True))'''
        self.add_trait("fwhm",
                       traits.List([6, 6, 6],
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
# test start
        self.add_trait("test_files",
                       OutputMultiPath(File(),
                                       optional=True,
                                       output=True,
                                       desc='Output just for a '
                                            'plug link detection test'))        
# test end

        self.process = spm.Smooth()
        self.change_dir = True

    def list_outputs(self, plugs=None):

        print('\ntest Passing through populse_mia/processes/User_processes/preprocess/spm/spatial_preprocessing.py, list_outputs method, Smooth(Process_Mia) class')

        super(Smooth, self).list_outputs()

        if not self.in_files:
            return {}
        else:
            self.process.inputs.in_files = self.in_files

        if self.out_prefix:
            self.process.inputs.out_prefix = self.out_prefix

        outputs = self.process._list_outputs()
        
        if plugs['test_files'].links_to:                         #test
            print('\ntest plug test_files is linked in output')  #test
            print('\ntest self.test_files: ', self.test_files)   #test
            outputs['test_files'] = ['test_files.nii']           #test test_files.nii doesn't mean anything, it's just an example.
            
        if not plugs['test_files'].links_to: #test
            print('\ntest plug test_files is not linked in output')      #test
            print('test ... So I do nothing or reset to <undefined>')     #test
            if self.test_files not in ["<undefined>", traits.Undefined]: #test
                self.test_files = "<undefined>"                          #test

        """raw_data_folder = os.path.join("data", "raw_data")
        derived_data_folder = os.path.join("data", "derived_data")
        for out_name, out_value in outputs.items():
            if type(out_value) is list:
                for idx, element in enumerate(out_value):
                    # To change the raw_data_folder to the derived_data_folder
                    element = element.replace(raw_data_folder, derived_data_folder)
                    outputs[out_name][idx] = element

        self.smoothed_files = outputs["smoothed_files"]"""

        inheritance_dict = {}
        for key, values in outputs.items():

            if key == "smoothed_files":
                for fullname in values:
                    path, filename = os.path.split(fullname)
                    
                    if self.out_prefix:
                        filename_without_prefix = filename[len(self.out_prefix):]
                        
                    else:
                        filename_without_prefix = filename[len('s'):]

                    if os.path.join(path, filename_without_prefix) in self.in_files:
                        inheritance_dict[fullname] = os.path.join(path, filename_without_prefix)

        return outputs, inheritance_dict

    def run_process_mia(self):

        print('\nprout Passing through populse_mia/processes/User_processes/preprocess/spm/spatial_preprocessing.py, run_process_mia method, Smooth(Process_Mia) class')
        
        super(Smooth, self).run_process_mia()

        for idx, element in enumerate(self.in_files):
            full_path = os.path.relpath(element)
            self.in_files[idx] = full_path
        self.process.inputs.in_files = self.in_files
        self.process.inputs.fwhm = self.fwhm
        self.process.inputs.data_type = self.data_type
        self.process.inputs.implicit_masking = self.implicit_masking
        self.process.inputs.out_prefix = self.out_prefix


        if self.test_files in ["<undefined>", traits.Undefined]: #test
            print('\ntest the plug test_files is not linked !')  #test

        else:                                                    #test
            print('\ntest the plug test_files is linked !')      #test
            
        self.process.run()
