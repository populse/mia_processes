# -*- coding: utf-8 -*- #

"""The spm preprocess library of the mia_processes package.

The purpose of this module is to customise the main spm preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - Coregister
        - GM_WM_Normalize
        - NewSegment
        - Normalize12
        - Realign
        - SliceTiming
        - Smooth

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# capsul import
from capsul.api import capsul_engine

# mia_processes import
#from .nipype_extension import NewSegmentMia

# nipype imports
from nipype.interfaces import spm
from nipype.interfaces.base import (File, InputMultiPath, InputMultiObject,
                                    OutputMultiPath, traits, TraitListObject,
                                    traits_extension)
from nipype.interfaces.spm.base import ImageFileSPM

# populse_mia imports
from populse_mia.data_manager.project import COLLECTION_CURRENT
from populse_mia.software_properties import Config
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# soma-base import
from soma.controller.trait_utils import relax_exists_constraint
from soma.qt_gui.qt_backend.Qt import QMessageBox

# Other imports
import itertools, math, os
import nibabel as nib
import numpy as np
import tempfile
from pathlib import Path
from traits.api import (Bool, Either, Enum, Float, Int, List, Range, String,
                        Tuple, Undefined)


class Coregister(ProcessMIA):
    """
    *Align together scans of different modalities*

    Please, see the complete documention for the `Coregister brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/spm/Coregister.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Coregister, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['spm', 'nipype']

        # Inputs description
        target_desc = ('The reference file to register to. An existing,'
                       ' uncompressed file (valid extensions:'
                       ' [.img, .nii, .hdr]).')
        source_desc = ('File to register to target image. A list of items which'
                       ' are an existing, uncompressed file (valid extensions:'
                       ' [.img, .nii, .hdr]).')
        apply_to_files_desc = ('Files to apply transformation to (a list of'
                               ' items which are an existing file name, (valid'
                               ' extensions: [.img, .nii, .hdr]).')
        jobtype_desc = "One of 'estwrite' or 'estimate' or 'write'."
        cost_function_desc = "One of 'mi' or 'nmi' or 'ecc' or 'ncc'."
        separation_desc = ('Sampling separation in mm (a list of items which'
                           ' are a float).')
        tolerance_desc = ('The acceptable tolerance for each of 12 params (a'
                          ' list of items which are a float).')
        fwhm_desc = ('Gaussian smoothing kernel width (mm) applied to the'
                     ' 256*256 joint histogram (a list of 2 items which are'
                     ' a float).')
        write_interp_desc = ('0 <= a long integer <= 7. The method by which the'
                             ' images are sampled when being written in a'
                             ' different space.')
        write_wrap_desc = ('A list of from 3 to 3 items which are an integer'
                           ' (int or long).')
        write_mask_desc = 'Mask output image (a boolean)'
        out_prefix_desc = 'Coregisterd output prefix (a string).'
        
        # Outputs description
        coregistered_source_desc = ("Coregistered source files, corresponding "
                                    "to 'source' images, (a pathlike object "
                                    "or string representing a file, or a list "
                                    "of pathlike objects or strings "
                                    "representing a file).")
        coregistered_files_desc = ("Coregistered other files, corresponding"
                                   " to 'apply_to_files', (a pathlike object "
                                    "or string representing a file, or a list "
                                    "of pathlike objects or strings "
                                    "representing a file).")

        # Inputs traits 
        self.add_trait("target",
                       ImageFileSPM(copyfile=False,
                                    output=False,
                                    optional=False,
                                    desc=target_desc))
        
        self.add_trait("source",
                       InputMultiPath(Either(ImageFileSPM(),
                                             Undefined),
                                      copyfile=True,
                                      output=False,
                                      optional=False,
                                      desc=source_desc))
        
        self.add_trait("apply_to_files",
                       InputMultiPath(Either(ImageFileSPM(),
                                             Undefined),
                                      copyfile=True,
                                      output=False,
                                      optional=True,
                                      desc=apply_to_files_desc))

        self.add_trait("jobtype",
                       Enum('estimate',
                            'estwrite',
                            'write',
                            output=False,
                            optional=True,
                            desc=jobtype_desc))

        self.add_trait("cost_function",
                       Enum('nmi',
                            'mi',
                            'ecc',
                            'ncc',
                            output=False,
                            optional=True,
                            desc=cost_function_desc))

        self.add_trait("separation",
                       List(value=[4.0, 2.0],
                            trait=Range(low=0.0, high=None),
                            minlen=1,
                            maxlen=32,
                            output=False,
                            optional=True,
                            desc=separation_desc))

        self.add_trait("tolerance",
                       List(value=[.02, .02, .02, 0.001, 0.001, 0.001,
                                   .01, .01, .01, 0.001, 0.001, 0.001],
                            trait=Range(low=0.0, high=None),
                            minlen=12,
                            maxlen=12,
                            output=False,
                            optional=True,
                            desc=tolerance_desc))

        self.add_trait("fwhm",
                       List(value=[7.0, 7.0],
                            trait=Range(low=0.0, high=None),
                            minlen=2,
                            maxlen=2,
                            output=False,
                            optional=True,
                            desc=fwhm_desc))

        self.add_trait("write_interp",
                       Range(value=4,
                             low=0,
                             high=7,
                             output=False,
                             optional=True,
                             desc=write_interp_desc))

        self.add_trait("write_wrap",
                       List(value=[0, 0, 0],
                            trait=Range(low=0, high=1),
                            minlen=3,
                            maxlen=3,
                            output=False,
                            optional=True,
                            desc=write_wrap_desc))

        self.add_trait("write_mask",
                       Bool(default_value=False,
                            output=False,
                            optional=True,
                            desc=write_mask_desc))

        self.add_trait("out_prefix",
                       String('r',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        # Output traits
        self.add_trait("coregistered_source",
                       OutputMultiPath(ImageFileSPM(),
                                       output=True,
                                       desc=coregistered_source_desc))

        self.add_trait("coregistered_files",
                       OutputMultiPath(ImageFileSPM(),
                                       output=True,
                                       optional=True,
                                       desc=coregistered_files_desc))

        self.init_default_traits()
        self.init_process('nipype.interfaces.spm.Coregister')

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
        super(Coregister, self).list_outputs()
        
        # Outputs definition and tags inheritance (optional)
        if (self.target and self.target != Undefined and
              self.source and self.source != [Undefined] and
              self.jobtype):
            self.process.source = self.source
            self.process.target = self.target
            self.process.jobtype = self.jobtype

            if self.apply_to_files and self.apply_to_files != [Undefined]:
                self.process.apply_to_files = self.apply_to_files

            if self.out_prefix:
                self.process.out_prefix = self.out_prefix

            # The management of self.process.output_directory could be delegated
            # to the populse_mia.user_interface.pipeline_manager.process_mia
            # module. We can't do it at the moment because the
            # sync_process_output_traits() of the capsul/process/nipype_process
            # module raises an exception in nipype if the mandatory parameters
            # are not yet defined!
            if self.output_directory:
                self.process.output_directory = self.output_directory

            else:
                print('No output_directory was found...!\n')

            self.outputs[
                'coregistered_source'
                        ] = self.process._coregistered_source
            self.outputs[
                'coregistered_files'
                        ] = self.process._coregistered_files

        if self.outputs:
            out_coregsource = []
            out_coregfiles = []

            for key, val in self.outputs.items():

                if key == "coregistered_source":
                    
                    if not isinstance(val, list):
                        val = [val]

                    for in_val, out_val in zip(self.source, val):
                        pathOval, fileOval = os.path.split(out_val)
                        _, fileIval = os.path.split(in_val)

                        if not self.jobtype == "estimate":

                            if self.out_prefix:
                                fileOvalNoPref = fileOval[len(self.out_prefix):]
                        
                            else:
                                fileOvalNoPref = fileOval[len('r'):]

                        else:
                            fileOvalNoPref = fileOval

                        if fileOvalNoPref == fileIval:
                            self.inheritance_dict[out_val] = in_val

                            if self.jobtype == "estwrite":
                                self.inheritance_dict[
                                    os.path.join(pathOval,
                                                 fileOvalNoPref)] = in_val
                                out_coregsource.append(out_val)
                                out_coregsource.append(
                                                os.path.join(pathOval,
                                                             fileOvalNoPref))

                if key == "coregistered_files" and val != Undefined:

                    if not isinstance(val, list):
                        val = [val]

                    for in_val, out_val in zip(self.apply_to_files, val):
                        pathOval, fileOval = os.path.split(out_val)
                        _, fileIval = os.path.split(in_val)

                        if not self.jobtype == "estimate":

                            if self.out_prefix:
                                fileOvalNoPref = fileOval[len(self.out_prefix):]
                        
                            else:
                                fileOvalNoPref = fileOval[len('r'):]

                        else:
                            fileOvalNoPref = fileOval

                        if fileOvalNoPref == fileIval:
                            self.inheritance_dict[out_val] = in_val

                            if self.jobtype == "estwrite":
                                self.inheritance_dict[
                                    os.path.join(pathOval,
                                                 fileOvalNoPref)] = in_val
                                out_coregfiles.append(out_val)
                                out_coregfiles.append(
                                               os.path.join(pathOval,
                                                            fileOvalNoPref))

            if out_coregsource:
                self.outputs["coregistered_source"] = out_coregsource

            if out_coregfiles:
                self.outputs["coregistered_files"] = out_coregfiles

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Coregister, self).run_process_mia()
        self.process.target = self.target
        self.process.source = self.source
        self.process.jobtype = self.jobtype

        if self.apply_to_files and self.apply_to_files != [Undefined]:
            self.process.apply_to_files = self.apply_to_files

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        self.process.cost_function = self.cost_function
        self.process.separation = self.separation
        self.process.tolerance = self.tolerance
        self.process.fwhm = self.fwhm
        self.process.write_interp = self.write_interp
        self.process.write_wrap = self.write_wrap
        self.process.write_mask = self.write_mask
        return self.process.run(configuration_dict={})


class GM_WM_Normalize(ProcessMIA):
    """
    *Normalisation of the grey and/or white matter map(s)*

    Please, see the complete documention for the `GM_WM_Normalize brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/spm/GM_WM_Normalize.html>`_

    """
    
    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(GM_WM_Normalize, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['spm', 'nipype']

        # Inputs description
        deformation_file_desc = ('File y_*.nii containing 3 deformation fields '
                                 'for the deformation in x, y and z dimension '
                                 '(an uncompressed file; valid extensions in '
                                 '[.img, .nii, .hdr]).')
        apply_to_files_desc = ('Files to apply transformation to. A list of '
                               'items which are an existing, uncompressed file '
                               '(valid extensions: [.img, .nii, .hdr]).')

        in_filter_desc = ('One of "GM" (gray matter) or "WM" (white matter) or '
                          '"GM & WM" (uses GM and WM; 2 files) or "GM + WM" '
                          '(makes GM + WM, then uses the result; 1 file).')

        write_bounding_box_desc = ('The bounding box (in mm) of the volume '
                                  'which is to be written (a list of 2 items, '
                                  'which are a list of items, which are a '
                                  'float.')
        write_voxel_sizes_desc = ('The voxel sizes (x, y & z, in mm) of the '
                                  'written normalised images (a list of 3 '
                                  'items, which are a float).')
        write_interp_desc = ('Degree of b-spline used for interpolation '
                             '(0 <= a long integer <= 7; 1 is OK for PET, '
                             'realigned fMRI, or segmentations).')
        
        # Outputs description
        normalized_files_desc = ('Normalised files (a pathlike object or '
                                 'string representing a file, or a list of '
                                 'pathlike objects or strings representing a '
                                 'file).')
        
        # Inputs traits 
        self.add_trait("deformation_file",
                       ImageFileSPM(output=False,
                                    optional=False,
                                    desc=deformation_file_desc))
 
        self.add_trait("apply_to_files",
                       InputMultiPath(traits.Either(ImageFileSPM(),
                                                    traits.List(ImageFileSPM()),
                                                    Undefined),
                                      output=False,
                                      optional=False,
                                      desc=apply_to_files_desc))

        self.add_trait("in_filter",
                       Enum('GM',
                            'WM',
                            'GM & WM',
                            'GM + WM',
                            output=False,
                            optional=True,
                            desc=in_filter_desc))

        self.add_trait("write_bounding_box",
                       traits.List(traits.List(traits.Float()),
                                   value=[[-78, -112, -50], [78, 76, 85]],
                                   output=False,
                                   optional=True,
                                   desc=write_bounding_box_desc))

        self.add_trait("write_voxel_sizes",
                       traits.List(traits.Float(),
                                   value=[1, 1, 1],
                                   output=False,
                                   optional=True,
                                   desc=write_voxel_sizes_desc))

        self.add_trait("write_interp",
                       traits.Range(value=1,
                                    low=0,
                                    high=7,
                                    output=False,
                                    optional=True,
                                    desc=write_interp_desc))

        self.add_trait("files",
                       traits.List(output=False,
                                   optional=True,
                                   userlevel=1,
                                   desc='internal parameter'))

        # Outputs traits
        self.add_trait("normalized_files",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=normalized_files_desc))

        self.init_default_traits()
        self.init_process('nipype.interfaces.spm.Normalize12')
        
    def _namesFilter(self):
        """Filtering of apply_to_files input parameter for GM, WM, etc.

        :returns: from the apply_to_files input parameter, returns the selected 
                  mask(s).
        """
        files = []

        if self.runFlag is False:

            for file_name in self.apply_to_files:
            
                if type(file_name) in [list, TraitListObject]:
                    file_name = file_name[0]

                if (self.in_filter == 'GM' and
                        os.path.basename(os.path.normpath(
                                                     file_name))[:2] in ['c1']):
                    files.append(file_name)

                elif (self.in_filter == 'WM' and
                        os.path.basename(os.path.normpath(
                                                     file_name))[:2] in ['c2']):
                    files.append(file_name)         

                elif ((self.in_filter == 'GM & WM' or
                        self.in_filter == 'GM + WM') and
                            os.path.basename(os.path.normpath(
                                               file_name))[:2] in ['c1', 'c2']):
                    files.append(file_name)

            if len(files) == 0:
                print('\nWarning: No GM or WM were detected ... Please, check '
                      'the input data and/or the in_filter parameter!')
                return []

            self.files = files.copy()
    
            if self.in_filter == 'GM + WM':
                j = []
                self.files = []
                
                while files:
                    i = files.pop()
                    i_path, i_file =  os.path.split(i)
                
                    if i_file[:2] in ['c1', 'c2']:

                        try:
                            ind = [os.path.split(k)[1] for k in files].index(
                                                              'c2' + i_file[2:])

                        except ValueError:
                            pass

                        else:
                            self.files.append(i)
                            self.files.append(files.pop(ind))
                            j.append(os.path.join(i_path, 'c1c2' + i_file[2:]))
                        
                        try:
                            ind = [os.path.split(k)[1] for k in files].index(
                                                              'c1' + i_file[2:])

                        except ValueError:
                            pass

                        else:
                            self.files.append(i)
                            self.files.append(files.pop(ind))
                            j.append(os.path.join(i_path, 'c1c2' + i_file[2:]))

                files = j
               
                if len(files) == 0:
                    print('\nWarning: no c1+c2 doublet was detected ... '
                          'Please, check the input data and/or the in_filter '
                          'parameter!')
                    return []

        if self.runFlag is True:
        
            if self.in_filter == 'GM + WM':

                for i, j in zip(self.files[0::2], self.files[1::2]):
                    img1 = nib.load(i)
                    img2 = nib.load(j)

                    if (img1.affine != img2.affine).all():
                        print('\nWarning! The affine arrays for {0} and {1} '
                              'are not the same. This can produce an erroneous '
                              'result!'.format(i, j))

                    img1_data = img1.get_fdata()
                    img2_data = img2.get_fdata()
            
                    if img1_data.shape != img2_data.shape:
                        print('\nWarning! The shapes for {0} and {1} are not '
                              'the same. This can produce an erroneous '
                              'result!'.format(i, j))

                    c1c2_img_data = np.add(img1_data, img2_data)
                    c1c2_img = nib.Nifti1Image(c1c2_img_data,
                                               header=img1.header,
                                               affine=img1.affine)
                    c1c2_img.header['descrip'] = np.array(b'Tissue class 1 + '
                                                          b'Tissue class 2',
                                                          dtype='|S80')
                    temp_dir = tempfile.mkdtemp()
                    c1c2_file = os.path.join(temp_dir,
                                             'c1c2' + os.path.split(i)[1][2:])
                    nib.save(c1c2_img, c1c2_file)
                    files.append(c1c2_file)

            else:
                files = self.files   

        return files
    
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
        super(GM_WM_Normalize, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        self.runFlag = False
        self.process.jobtype = 'write'
        self.process.apply_to_files = self._namesFilter()

        if self.process.apply_to_files != Undefined and self.deformation_file:
            self.process.deformation_file = self.deformation_file

            # The management of self.process.output_directory could be delegated
            # to the populse_mia.user_interface.pipeline_manager.process_mia
            # module. We can't do it at the moment because the
            # sync_process_output_traits() of the capsul/process/nipype_process
            # module raises an exception in nipype if the mandatory parameters
            # are not yet defined!
            if self.output_directory:
                self.process.output_directory = self.output_directory

            else:
                print('No output_directory was found...!\n')

            self.outputs['normalized_files'] = self.process._normalized_files

        if self.outputs:

            for key, val in self.outputs.items():

                if (key == "normalized_files"):

                    if self.in_filter == 'GM + WM':
                        self.inheritance_dict[val] = self.apply_to_files[0]
                        
                    else:

                        if not isinstance(val, list):
                            val = [val]
                        
                        for (in_val,
                              out_val) in zip(self.process.apply_to_files, val):
                            _, fileOval = os.path.split(out_val)
                            _, fileIval = os.path.split(in_val)
                            fileOvalNoPref = fileOval[1:]

                            if fileOvalNoPref == fileIval:
                                self.inheritance_dict[out_val] = in_val

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()
 
    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        self.runFlag = True
        self.process.jobtype = 'write'
        self.process.apply_to_files = self._namesFilter()
        self.process.deformation_file = self.deformation_file
        self.process.write_bounding_box = self.write_bounding_box
        self.process.write_voxel_sizes = self.write_voxel_sizes
        self.process.write_interp = self.write_interp
        self.process.trait('image_to_align').optional = True

        # because the sync_process_output_traits() of the capsul/process/nipype_process
        # module raises an exception in nipype if the mandatory parameters
        # are not yet defined, the next line can't be write before!
        super(GM_WM_Normalize, self).run_process_mia()
        return self.process.run(configuration_dict={})


class NewSegment(ProcessMIA):
    """
    *Segmentation: Segments, bias corrects and spatially normalises - all in the same model*

    Please, see the complete documention for the `NewSegment brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/spm/NewSegment.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(NewSegment, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['spm', 'nipype']

        # Inputs description
        
        channel_files_desc = ('Path of the scans for processing, valid '
                              'extensions: [.img, .nii, .hdr]. A list with one '
                              'string element corresponding to an existing '
                              'path file.')
        channel_info_desc = ('A tuple with the following values: bias '
                             'reguralisation -a float between 0 and 10-, bias '
                             'FWHM -a float between 20 and infinity-, which '
                             'maps to save -a tuple of two boolean values '
                             '(estimated bias field, bias corrected image)-.')
        tissues_desc = ('A list of tuples with the following values for each '
                        'tissue types (grey matter, white matter, etc.): '
                        '(tissue probability map (4D), 1-based index to frame),'
                        ' number of gaussians, (which maps to save: Native,'
                        'DARTEL), (which maps to save: Unmodulated, Modulated),'
                        ' ...')
        warping_regularization_desc = ('Define of the roughness of the '
                                       'deformations for registration using '
                                       '1 or 5 elements: a float or a list of '
                                       'floats, the latter is required '
                                       'by SPM12.')
        affine_regularization_desc = ('Standard space for affine registration: '
                                      'mni, eastern, subj or none.')

        sampling_distance_desc = ('Approximate distance between sampled points '
                                 'when estimating the model parameters: '
                                 'a float.')

        write_deformation_fields = ('Which deformation fields, that can be used'
                                    ' by the deformation utility, to save: a '
                                    'list of 2 booleans for Inverse, Forward.')

        # Outputs description
        bias_corrected_images_desc = ('Bias corrected images (a pathlike '
                                      'object or string representing a file, '
                                      'or a list of pathlike objects or '
                                      'strings representing a file).')
        bias_field_images_desc = ('Estimated bias field (a pathlike object or'
                                  'string representing a file, or a list of '
                                  'pathlike objects or strings representing a '
                                  'file).')
        native_class_images_desc = ('Native space probability maps (a list of '
                                    'items which are a list of items which are '
                                    'a pathlike object or string representing '
                                    'a file).')
        dartel_input_images_desc = ('“Imported” class images into a form that '
                                    'can be used with the Dartel toolbox (a '
                                    'list of items which are a list of items '
                                    'which are a pathlike object or string '
                                    'representing a file).')
        modulated_class_images_desc = ('Modulated and normalised class images '
                                       '(a list of items which are a list of '
                                       'items which are a pathlike object or '
                                       'string representing a file).')
        normalized_class_images_desc = ('Normalised class images, without '
                                        'modulation (a list of items which are '
                                        'a list of items which are a pathlike '
                                        'object or string representing a file')
        inverse_deformation_field_desc = ('Inverse deformation field (a '
                                          'pathlike object or string '
                                          'representing a file, or a list of '
                                          'pathlike objects or strings '
                                          'representing a file).')
        forward_deformation_field_desc = ('Forward deformation field (a '
                                          'pathlike object or string '
                                          'representing a file).')
        transformation_mat_desc = ('Normalisation transformation (a pathlike '
                                   'object or string representing a file).')     

        # Tissues parameter definition
        config = Config()
        resources_path = os.path.join(config.get_mia_path(), 'resources')
        tpm_path = os.path.join(resources_path, 'spm12', 'tpm', 'TPM.nii')

        if not Path(tpm_path).is_file():
            print('\nNewSegment brick warning!: The {} file seems to '
                  'not exists ...'.format(tpm_path))
            tissues_list = Undefined

        else:
            tissues_list = [((tpm_path, 1), 2, (True, False), (False, False)),
                            ((tpm_path, 2), 2, (True, False), (False, False)),
                            ((tpm_path, 3), 2, (True, False), (False, False)),
                            ((tpm_path, 4), 3, (True, False), (False, False)),
                            ((tpm_path, 5), 4, (True, False), (False, False)),
                            ((tpm_path, 6), 2, (True, False), (False, False))]

        # Inputs traits
        self.add_trait("channel_files",
                       InputMultiPath(ImageFileSPM(),
                                      output=False,
                                      optional=False,
                                      desc=channel_files_desc))

        self.add_trait("channel_info",
                       Tuple(Range(value=0.0001,
                                   low=0.,
                                   high=10.),
                             Range(value=60.,
                                   low=20.,
                                   high=None),
                             Tuple(Bool(False),
                                   Bool(True)),
                             output=False,
                             optional=True,
                             desc=channel_info_desc))

        self.add_trait("tissues",
                       Either(
                           List(
                               Tuple(
                                   Tuple(ImageFileSPM(exists=True),
                                         Int()),
                                   Int(),
                                   Tuple(Bool, Bool),
                                   Tuple(Bool, Bool))),
                           Undefined,
                           output=False,
                           optional=True,
                           desc=tissues_desc))
        self.tissues = tissues_list

        self.add_trait("warping_regularization",
                       Either(List(Float(),
                                   minlen=5,
                                   maxlen=5),
                              Float(),
                              default=[0, 0.001, 0.5, 0.05, 0.2],
                              output=False,
                              optional=True,
                              desc=warping_regularization_desc))

        self.add_trait("affine_regularization",
                       Enum('mni',
                            'eastern',
                            'subj',
                            'none',
                            output=False,
                            optional=True,
                            desc=affine_regularization_desc))
        
        self.add_trait("sampling_distance",
                       Float(3.0,
                             output=False,
                             optional=True,
                             desc=sampling_distance_desc))

        self.add_trait("write_deformation_fields",
                       List(Bool(),
                            minlen=2,
                            maxlen=2,
                            value=[False, True],
                            output=False,
                            optional=True,
                            desc=write_deformation_fields))

        # Output traits
        self.add_trait("bias_corrected_images",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=bias_corrected_images_desc))

        self.add_trait("bias_field_images",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=bias_field_images_desc))

        self.add_trait("native_class_images",
                       List(List(File()),
                            output=True,
                            optional=True,
                            desc=native_class_images_desc))

        self.add_trait("dartel_input_images",
                       List(List(File()),
                            output=True,
                            optional=True,
                            desc=dartel_input_images_desc))

        self.add_trait("modulated_class_images",
                       List(List(File()),
                            output=True,
                            optional=True,
                            desc=modulated_class_images_desc))

        self.add_trait("normalized_class_images",
                       List(List(File()),
                            output=True,
                            optional=True,
                            desc=normalized_class_images_desc))

        self.add_trait("inverse_deformation_field",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=inverse_deformation_field_desc))

        self.add_trait("forward_deformation_field",
                       OutputMultiPath(ImageFileSPM(),
                                       output=True,
                                       optional=True,
                                       desc=forward_deformation_field_desc))

        self.add_trait("transformation_mat",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=transformation_mat_desc))

        self.init_default_traits()
        self.init_process('nipype.interfaces.spm.NewSegment')

        """ To be removed for V2
        # process instanciation
        self.process = NewSegmentMia() # workaround to try to decrease the
                                       # instantiation time of the NewSegment
                                       # class
        #self.process = spm.NewSegment()
        """

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
        super(NewSegment, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if (self.channel_files and self.channel_info and
                       self.write_deformation_fields and self.tissues):
            self.process.channel_files = self.channel_files
            self.process.channel_info = self.channel_info
            (self.process.
                  write_deformation_fields) = self.write_deformation_fields
            self.process.tissues = self.tissues

            # The management of self.process.output_directory could be delegated
            # to the populse_mia.user_interface.pipeline_manager.process_mia
            # module. We can't do it at the moment because the
            # sync_process_output_traits() of the capsul/process/nipype_process
            # module raises an exception in nipype if the mandatory parameters
            # are not yet defined!
            if self.output_directory:
                self.process.output_directory = self.output_directory

            else:
                print('No output_directory was found...!\n')

            for k in ('bias_corrected_images', 'bias_field_images',
                      'native_class_images', 'dartel_input_images',
                      'modulated_class_images', 'normalized_class_images',
                      'inverse_deformation_field', 'forward_deformation_field',
                      'transformation_mat'):
                self.outputs[k] = getattr(self.process, '_' + k)
        
        """
         When there is only one image at the input, the tags inheritance is
         directly managed in PipelineManagerTab.add_plug_value_to_database().
         However the inheritance of 2 output parameters is defined below, as
         an example
        """
        if self.outputs:

            for key, val in self.outputs.items():

                if key == "native_class_images":

                    if val[0]:
                        inp = self.channel_files.copy() * len(val)
                        out = [x[0] for x in val]

                        for in_val, out_val in zip(inp, out):
                            _, fileOval = os.path.split(out_val)
                            _, fileIval = os.path.split(in_val)
                            fileOvalNoPref = fileOval[2:]

                            if fileOvalNoPref == fileIval:
                                self.inheritance_dict[out_val] = in_val

                if key == "forward_deformation_field":

                    if val != Undefined:
                        val = [val]

                        for in_val, out_val in zip(self.channel_files, val):
                            _, fileOval = os.path.split(out_val)
                            _, fileIval = os.path.split(in_val)
                            fileOvalNoPref = fileOval[2:]

                            if fileOvalNoPref == fileIval:
                                self.inheritance_dict[out_val] = in_val

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()
    
    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(NewSegment, self).run_process_mia()
        self.process.channel_files = self.channel_files
        self.process.channel_info = self.channel_info
        self.process.tissues = self.tissues
        self.process.warping_regularization = self.warping_regularization
        self.process.affine_regularization = self.affine_regularization
        self.process.sampling_distance = self.sampling_distance
        self.process.write_deformation_fields = self.write_deformation_fields
        return self.process.run(configuration_dict={})


class Normalize12(ProcessMIA):
    """    
    *Computes the warp that best aligns the template (atlas) to the individual’s image*

    Please, see the complete documention for the `Normalize12 brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/spm/Normalize12.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Normalize12, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['spm', 'nipype']

        # Inputs description
        image_to_align_desc = ('The image that the atlas data is warped into '
                               'alignment with (an existing file; valid '
                               'extensions in [.img, .nii, .hdr]). Mutually '
                               'exclusive with the deformation_file parameter.')
        deformation_file_desc = ('File y_*.nii containing 3 deformation fields '
                                 'for the deformation in x, y and z dimension '
                                 '(an uncompressed file; valid extensions in '
                                 '[.img, .nii, .hdr]). Mutually exclusive with '
                                 'the image_to_align and tpm parameters.')
        apply_to_files_desc = ('Files to apply transformation to. A list of '
                               'items which are an existing, uncompressed file '
                               '(valid extensions: [.img, .nii, .hdr]).')
        jobtype_desc = 'One of "estwrite" or "estimate" or "write".'
        bias_regularization_desc =('(For low-intensity non-uniformity artifacts'
                                   ', use a high bias regularization (a float '
                                   ' between 0 and 10).')
        bias_fwhm_desc = ('Full Width at Half Maximum of Gaussian smoothness '
                          'of bias (a value in [30, 40, 50, 60, 70, 80, 90, '
                          '100, 110, 120, 130, 140, 150, ‘Inf’).')
        tpm_desc = ('The template in form of tissue probability atlas (a '
                    'pathlike object or string representing an existing file). '
                    'Mutually exclusive with the deformation_file parameter.')
        affine_regularization_type_desc = ("Standard space for affine "
                                           "registration (‘mni’ or ‘size’ or "
                                           "‘none’).")
        warping_regularization_desc = ('The measure of the roughness of the '
                                       'deformations for registration. Involve '
                                       'the sum of 5 elements (list of '
                                       'floats).')
        smoothness_desc = ('Value to smooth the data before normalisation (a '
                           'float; in mm). 0 is a good value for MRI.')
        sampling_distance_desc = ('Approximate distance between sampled points '
                                  'when estimating the model parameters (a '
                                  'float).')
        write_bounding_box_desc = ('The bounding box (in mm) of the volume '
                                  'which is to be written (a list of 2 items, '
                                  'which are a list of items, which are a '
                                  'float.')
        write_voxel_sizes_desc = ('The voxel sizes (x, y & z, in mm) of the '
                                  'written normalised images (a list of 3 '
                                  'items, which are a float).')
        write_interp_desc = ('Degree of b-spline used for interpolation '
                             '(0 <= a long integer <= 7; 1 is OK for PET, '
                             'realigned fMRI, or segmentations).')
        
        # Outputs description
        deformation_field_desc  = ('File y_*.nii containing 3 deformation '
                                   'fields for the deformation in x, y and z '
                                   'dimension (a pathlike object or string '
                                   'representing a file, or a list of pathlike '
                                   'objects or strings representing a file)')
        normalized_image_desc = ('Normalised file that needed to be aligned (a '
                                 'pathlike object or string representing a '
                                 'file, or a list of pathlike objects or '
                                 'strings representing a file).')
        normalized_files_desc = ('Normalised other files (a pathlike object or '
                                 'string representing a file, or a list of '
                                 'pathlike objects or strings representing a '
                                 'file).')

        # Tpm parameter definition
        config = Config()
        resources_path = os.path.join(config.get_mia_path(), 'resources')
        tpm_path = os.path.join(resources_path, 'spm12', 'tpm', 'TPM.nii')
        
        if not Path(tpm_path).is_file():
            print('\nThe {} file seems to not exists ...'.format(tpm_path))
            tpm_path = Undefined

        # Inputs traits
        self.add_trait("image_to_align",
                       ImageFileSPM(output=False,
                                    optional=True,
                                    desc=image_to_align_desc))

        self.add_trait("deformation_file",
                       ImageFileSPM(output=False,
                                    optional=True,
                                    desc=deformation_file_desc))
 
        self.add_trait("apply_to_files",
                       InputMultiPath(Either(ImageFileSPM(),
                                             List(ImageFileSPM()),
                                             Undefined),
                                      output=False,
                                      optional=True,
                                      desc=apply_to_files_desc))

        self.add_trait("jobtype",
                       Enum("write",
                            "est",
                            "estwrite",
                            output=False,
                            optional=True,
                            desc=jobtype_desc))

        self.add_trait("bias_regularization",
                       Enum(0, 0.00001, 0.0001, 0.001,
                            0.01, 0.1, 1, 10,
                            output=False,
                            optional=True,
                            desc=bias_regularization_desc))

        self.add_trait("bias_fwhm",
                       Enum(60, 30, 40, 50, 70, 80, 90, 100,
                            110, 120, 130, 140,150, "Inf",
                            output=False,
                            optional=True,
                            desc=bias_fwhm_desc))

        self.add_trait("tpm",
                       File(exists=True,
                            value=tpm_path,
                            output=False,
                            optional=True,
                            desc=tpm_desc))

        self.add_trait("affine_regularization_type",
                       Enum("mni",
                            "size",
                            "none",
                            output=False,
                            optional=True,
                            desc=affine_regularization_type_desc))

        self.add_trait("warping_regularization",
                       List(value=[0, 0.001, 0.5, 0.05, 0.2],
                            trait=Float(),
                            minlen=5,
                            maxlen=5,
                            output=False,
                            optional=True,
                            desc=warping_regularization_desc))

        self.add_trait("smoothness",
                       Float(0.,
                             output=False,
                             optional=True,
                             desc=smoothness_desc))

        self.add_trait("sampling_distance",
                       Float(3.,
                             output=False,
                             optional=True,
                             desc=sampling_distance_desc))

        self.add_trait("write_bounding_box",
                       List(List(Float()),
                            value=[[-78, -112, -50], [78, 76, 85]],
                            output=False,
                            optional=True,
                            desc=write_bounding_box_desc))

        self.add_trait("write_voxel_sizes",
                       List(Float(),
                            value=[1, 1, 1],
                            output=False,
                            optional=True,
                            desc=write_voxel_sizes_desc))

        self.add_trait("write_interp",
                       Range(value=1,
                             low=0,
                             high=7,
                             output=False,
                             optional=True,
                             desc=write_interp_desc))

        # Outputs traits
        self.add_trait("deformation_field",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=deformation_field_desc))

        self.add_trait("normalized_files",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=normalized_files_desc))

        self.init_default_traits()
        self.init_process('nipype.interfaces.spm.Normalize12')
    

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
        super(Normalize12, self).list_outputs()

        # Outputs definition
        _flag = False
        
        if ((self.apply_to_files) and (self.apply_to_files != [Undefined]) and
                         (self.deformation_file) and (self.jobtype == 'write')):
            self.process.jobtype = self.jobtype
            self.process.apply_to_files = self.apply_to_files
            self.process.deformation_file = self.deformation_file
            _flag = True

            if self.image_to_align:
                self.image_to_align = Undefined

            if self.tpm:
                self.tpm = Undefined
                
        elif ((self.image_to_align) and (self.image_to_align != Undefined) and
                                        (self.tpm) and (self.jobtype == 'est')):
            self.process.image_to_align = self.image_to_align
            self.process.tpm = self.tpm
            self.process.jobtype = self.jobtype
            _flag = True

            if self.apply_to_files and (self.apply_to_files != [Undefined]):
                self.apply_to_files = [Undefined]

            if self.deformation_file:
                self.deformation_file = Undefined
                
        elif ((self.image_to_align) and (self.image_to_align != Undefined) and
                                   (self.tpm) and (self.jobtype == 'estwrite')):
            self.process.image_to_align = self.image_to_align
            self.process.tpm = self.tpm
            self.process.jobtype = self.jobtype
            _flag = True

            if (self.apply_to_files) and (self.apply_to_files != [Undefined]):
                self.process.apply_to_files = self.apply_to_files

            if self.deformation_file:
                self.deformation_file = Undefined

        if _flag:

            # The management of self.process.output_directory could be delegated
            # to the populse_mia.user_interface.pipeline_manager.process_mia
            # module. We can't do it at the moment because the
            # sync_process_output_traits() of the capsul/process/nipype_process
            # module raises an exception in nipype if the mandatory parameters
            # are not yet defined!
            if self.output_directory:
                self.process.output_directory = self.output_directory

            else:
                print('No output_directory was found...!\n')

            for k in ('deformation_field', 'normalized_files'):
                self.outputs[k] = getattr(self.process, '_' + k)

        #Tags inheritance (optional)
        if self.outputs:

            for key, val in self.outputs.items():

                if (key == "normalized_files") and (val != Undefined):

                    if not isinstance(val, list):
                        val = [val]

                    for in_val, out_val in zip(self.apply_to_files, val):
                        pathOval, fileOval = os.path.split(out_val)
                        _, fileIval = os.path.split(in_val)
                        fileOvalNoPref = fileOval[1:]

                        if fileOvalNoPref == fileIval:
                            self.inheritance_dict[out_val] = in_val

                if (key == "deformation_field") and (val != Undefined):
                    pathOval, fileOval = os.path.split(val)
                    _, fileIval = os.path.split(self.image_to_align)
                    fileOvalNoPref = fileOval[2:]
                    
                    if fileOvalNoPref == fileIval:
                        self.inheritance_dict[val] = self.image_to_align

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        self.process.trait('image_to_align').optional = True
        self.process.image_to_align = self.image_to_align
        self.process.jobtype = self.jobtype
        self.process.trait('deformation_file').optional = True
        self.process.deformation_file = self.deformation_file
        
        if self.apply_to_files == [Undefined]:
            self.process.apply_to_files = Undefined

        else:
            self.process.apply_to_files = self.apply_to_files

        self.process.bias_regularization = self.bias_regularization
        self.process.bias_fwhm = self.bias_fwhm

        if self.jobtype.startswith('est'):
            self.process.tpm = self.tpm

        (self.process.
                   affine_regularization_type) = self.affine_regularization_type
        self.process.warping_regularization = self.warping_regularization
        self.process.smoothness = self.smoothness
        self.process.sampling_distance = self.sampling_distance
        self.process.write_bounding_box = self.write_bounding_box
        self.process.write_voxel_sizes = self.write_voxel_sizes
        self.process.write_interp = self.write_interp

        # because the sync_process_output_traits() of the capsul/process/nipype_process
        # module raises an exception in nipype if the mandatory parameters
        # are not yet defined, the next line can't be write before!
        super(Normalize12, self).run_process_mia()
        return self.process.run(configuration_dict={})


class Realign(ProcessMIA):
    """
    *Realigns a time-series of images acquired from the same subject*

    Please, see the complete documention for the `Realign brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/spm/Realign.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Realign, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['spm', 'nipype']
        
        # Inputs description
        in_files_desc = ('The images to realign (a list of pathlike objects or '
                         'strings representing a file or of list of items '
                         'which are a pathlike object or strings representing '
                         'a file or a _Undefined).')
        jobtype_desc = 'One of "estwrite" or "estimate" or "write".'
        quality_desc = ('Quality versus speed trade-off (0.0 <= a floating'
                        'point number <= 1.0).')
        separation_desc = ('Sampling separation (in mm) in the reference image'
                           '(a floating point number >= 0.0).')
        fwhm_desc = ('The gaussian smoothing kernel width (mm, a floating point'
                     ' number >= 0.0) applied to the images before estimating'
                    ' the realignment parameters.')
        register_to_mean_desc = ('Indicate whether realignment is done to the'
                                 ' mean image (True) or to the first image'
                                 ' (False).')
        interp_desc = ('Degree of b-spline (1 <= a long integer <= 7) used for'
                       ' interpolation.')
        wrap_desc = ('Check if interpolation should wrap in [x,y,z] (a list of'
                     ' 3 items which are integer int or long).')
        weight_img_desc = ('File name of the weighting image, to weight each'
                           ' voxel of the reference image.')
        write_which_desc = ('Determines which images to reslice (a list of'
                            ' items which are a value of class "int").')
        write_interp_desc = ('0 <= a long integer <= 7. The method by which the'
                             ' images are sampled when being written in a'
                             ' different space.')
        write_wrap_desc = ('A list of from 3 to 3 items which are an integer'
                           ' (int or long).')
        write_mask_desc = 'Mask output image (a boolean)'
        out_prefix_desc = 'Realigned output prefix (a string).'

        # Outputs description
        realigned_files_desc = ('If write_which is [2, 0], [1, 0] or [2, 1] and'
                                ' jobtype is write or estwrite, these will be '
                                'the resliced files (a pathlike object or '
                                'string representing a file or a list of '
                                'items which are a pathlike object or string '
                                'representing a file or a list of list of '
                                'items which are a pathlike object or string '
                                'representing a file).')
        modified_in_files_desc = ('If the jobtype parameter is estimate or '
                                  'estwrite, these will be copies of the '
                                  'in_files with a rewritten header (a '
                                  'pathlike object or string representing a '
                                  'file or a list of items which are a '
                                  'pathlike object or string representing a '
                                  'file or a list of list of items which are a '
                                  'pathlike object or string representing a '
                                  'file).')
        mean_image_desc = ('If write_which is [2, 1] or [0, 1] and jobtype is '
                           'write or estwrite, this will be the mean image '
                           'file from the realignment (a pathlike object or '
                           'string representing a file).')
        realignment_parameters_desc = ('If the jobtype parameter is estimate '
                                       'or estwrite, this will be the '
                                       'estimated translation and rotation '
                                       'parameters (a pathlike object or '
                                       'string representing a file, or a list '
                                       'of pathlike objects or strings '
                                       'representing a file).')
        
        # Inputs traits
        self.add_trait('in_files',
                       InputMultiPath(Either(ImageFileSPM(),
                                             List(ImageFileSPM()),
                                             Undefined),
                                      copyfile=True,
                                      output=False,
                                      optional=False,
                                      desc=in_files_desc))
        
        self.add_trait("jobtype",
                       Enum('estwrite',
                            'estimate',
                            'write',
                            output=False,
                            optional=True,
                            desc=jobtype_desc))

        self.add_trait("quality",
                       Range(value=0.9,
                             low=0.0,
                             high=1.0,
                             output=False,
                             optional=True,
                             desc=quality_desc))

        self.add_trait("separation",
                       Range(value=4.0,
                             low=0.0,
                             high=None,
                             output=False,
                             optional=True,
                             desc=separation_desc))

        self.add_trait("fwhm",
                       Range(value=5.0,
                             low=0.0,
                             high=None,
                             output=False,
                             optional=True,
                             desc=fwhm_desc))

        self.add_trait("register_to_mean",
                       Bool(True,
                            output=False,
                            optional=True,
                            desc=register_to_mean_desc))

        self.add_trait("interp",
                       Range(value=2,
                             low=0,
                             high=7,
                             output=False,
                             optional=True,
                             desc=interp_desc))

        self.add_trait("wrap",
                       List(value=[0, 0, 0],
                            trait=Range(low=0, high=1),
                            minlen=3,
                            maxlen=3,
                            output=False,
                            optional=True,
                            desc=wrap_desc))

        self.add_trait("weight_img",
                       File(value=Undefined,
                            output=False,
                            optional=True,
                            desc=weight_img_desc))

        self.add_trait('write_which',
                       Enum([2, 1],
                            [1, 0],
                            [2, 0],
                            [0, 1],
                            output=False,
                            optional=True,
                            desc=write_which_desc))
 
        self.add_trait("write_interp",
                       Range(value=4,
                             low=0,
                             high=7,
                             output=False,
                             optional=True,
                             desc=write_interp_desc))

        self.add_trait("write_wrap",
                       List(value=[0, 0, 0],
                            trait=Range(low=0, high=1),
                            minlen=3,
                            maxlen=3,
                            output=False,
                            optional=True,
                            desc=write_wrap_desc))

        self.add_trait("write_mask",
                       Bool(default_value=True,
                            output=False,
                            optional=True,
                            desc=write_mask_desc))

        self.add_trait("out_prefix",
                       String(value='r',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        # Output traits
        self.add_trait("realigned_files",
                       OutputMultiPath(Either(List(File()),
                                              File()),
                                       output=True,
                                       optional=True,
                                       desc=realigned_files_desc))

        self.add_trait("modified_in_files",
                       OutputMultiPath(Either(List(File()),
                                              File()),
                                       output=True,
                                       optional=True,
                                       desc=modified_in_files_desc))

        self.add_trait("mean_image",
                       File(output=True,
                            optional=True,
                            desc=mean_image_desc))
        
        self.add_trait("realignment_parameters",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=realignment_parameters_desc))


        self.init_default_traits()
        self.init_process('nipype.interfaces.spm.Realign')

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
        super(Realign, self).list_outputs()
        
        # Outputs definition and tags inheritance (optional)
        if self.in_files and self.in_files != [Undefined]:
            self.process.in_files = self.in_files

            if self.out_prefix:
                self.process.out_prefix = self.out_prefix
                
            # The management of self.process.output_directory could be delegated
            # to the populse_mia.user_interface.pipeline_manager.process_mia
            # module. We can't do it at the moment because the
            # sync_process_output_traits() of the capsul/process/nipype_process
            # module raises an exception in nipype if the mandatory parameters
            # are not yet defined!
            if self.output_directory:
                self.process.output_directory = self.output_directory

            else:
                print('No output_directory was found...!\n')

            if self.jobtype:
                self.process.jobtype = self.jobtype

            if self.write_which:
                self.process.write_which = self.write_which

            for k in ('realigned_files', 'modified_in_files',
                      'mean_image', 'realignment_parameters'):
                self.outputs[k] = getattr(self.process, '_' + k)

        if self.outputs:

            for key, val in self.outputs.items():

                if key == 'realigned_files':

                    if ((self.jobtype == 'estimate') or
                     (self.write_which == [0, 1] and
                     (self.jobtype == 'write' or  self.jobtype == 'estwrite'))):
                        self.outputs['realigned_files'] = Undefined

                    else:
                        
                        if not isinstance(val, list):
                            val = [val]
                            
                        for in_val, out_val in zip(self.in_files, val):
                            _, fileOval = os.path.split(out_val)
                            _, fileIval = os.path.split(in_val)

                            if self.out_prefix:
                                fileOvalNoPref = fileOval[len(self.out_prefix):]
                        
                            else:
                                fileOvalNoPref = fileOval[len('r'):]

                            if fileOvalNoPref == fileIval:
                                self.inheritance_dict[out_val] = in_val

                if key == 'modified_in_files':

                    if self.jobtype == 'write':
                        self.outputs['modified_in_files'] = Undefined

                    else:

                        if not isinstance(val, list):
                            val = [val]
                            
                        for in_val, out_val in zip(self.in_files, val):
                            _, fileOval = os.path.split(out_val)
                            _, fileIval = os.path.split(in_val)
                            fileOvalNoPref = fileOval

                            if fileOvalNoPref == fileIval:
                                self.inheritance_dict[out_val] = in_val

                # Only the first data in in_files gives mean_image
                if key == 'mean_image':

                    if (
                      (self.jobtype == 'estimate') or
                    ((self.write_which == [2, 0] or self.write_which == [1, 0])
                                             and
                     (self.jobtype == 'write' or  self.jobtype == 'estwrite'))
                       ):
                        self.outputs['mean_image'] = Undefined

                    else:

                        if not isinstance(val, list):
                            values = [val]

                        for in_val, out_val in zip(self.in_files, val):
                            _, fileOval = os.path.split(out_val)
                            _, fileIval = os.path.split(in_val)
                            fileOvalNoPref = fileOval[len('mean'):]

                            if fileOvalNoPref == fileIval:
                                self.inheritance_dict[out_val] = in_val

                if key == 'realignment_parameters':

                    if self.jobtype == 'write':
                        self.outputs['realignment_parameters'] = Undefined

                    else:

                        if not isinstance(val, list):
                            values = [val]

                        for in_val, out_val in zip(self.in_files, val):
                            _, fileOval = os.path.split(out_val)
                            _, fileIval = os.path.split(in_val)
                            fileOvalNoPref = fileOval[len('rp_'):]

                            if ('.'.join(fileOvalNoPref.split('.')[:-1]) ==
                                            '.'.join(fileIval.split('.')[:-1])):
                                self.inheritance_dict[out_val] = in_val

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Realign, self).run_process_mia()
        self.process.in_files = self.in_files

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        self.process.jobtype = self.jobtype
        self.process.quality = self.quality
        self.process.separation = self.separation
        self.process.fwhm = self.fwhm
        self.process.register_to_mean = self.register_to_mean
        self.process.interp = self.interp
        self.process.wrap = self.wrap
        self.process.weight_img = self.weight_img
        self.process.write_which = self.write_which
        self.process.write_interp = self.write_interp
        self.process.write_wrap = self.write_wrap
        self.process.write_mask = self.write_mask
        return self.process.run(configuration_dict={})


class SliceTiming(ProcessMIA):
    """
    *Temporal correction to get back every slice at the same acquisition time*

    Please, see the complete documentation for the `SliceTiming brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/spm/SliceTiming.html>`

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(SliceTiming, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['spm', 'nipype']

        # Inputs description
        desc_in_file = ("The images to be brought back to the reference slice "
                        "time (a list of pathlike objects or string "
                        "representing a file or of list of items which are a "
                        "pathlike object or strings representing a file or a "
                        "_Undefined or None)")
        desc_acquisition = ('Type of image acquisition (one of "sequential '
                            'ascending", "sequential descending", "interleaved '
                            '(middle-top)", "interleaved (bottom-up)", '
                            '"interleaved (top-down)"). From the acquisition '
                            'parameter the slice_order parameter can be '
                            'calculated automatically')
        desc_num_slices = ("Number of slices in the in_file (an integer). "
                           "Can be automatically retrieved from the database")
        desc_TR = ("Repetition time in the in_file (a float, in seconds). "
                   "Can be automatically retrieved from the database")
        desc_TA = ("Time of volume acquisition (a float, in seconds). Can be "
                   "automatically calculated as TR-(TR/num_slices)")
        desc_slice_order = ('Order of how were recorded the slices during the '
                            'volumes aquisition (a list of items which are '
                            'an integer or a float or a _Undefined or None. '
                            'Can be automatically calculated according to the '
                            'acquisition" parameter)')
        desc_ref_slice = ('The chosen reference slice, will serve as the '
                          'reference for the other slices during the '
                          'correction (an integer or a float or an _Undefined '
                          'or None). A default value is calculated according '
                          'to the "acquisition" parameter')
        desc_out_prefix = ("SliceTiming routine output prefix (a string)")

        # Outputs description
        desc_timed_files = ("The image after the SliceTiming correction (a "
                            "pathlike object or string representing a file)")

        # Input traits
        self.add_trait("in_files",
                       InputMultiPath(Either(ImageFileSPM(),
                                             List(ImageFileSPM()),
                                             Undefined),
                                      output=False,
                                      optional=False,
                                      desc=desc_in_file))

        self.add_trait("acquisition",
                       Enum("sequential ascending",
                            "sequential descending",
                            "interleaved (middle-top)",
                            "interleaved (bottom-up)",
                            "interleaved (top-down)",
                            value="sequential ascending",
                            output=False,
                            optional=True,
                            desc=desc_acquisition))

        self.add_trait("num_slices",
                       Either(Int(),
                              Undefined,
                              output=False,
                              optional=True,
                              desc=desc_num_slices))
        self.num_slices = Undefined

        self.add_trait("TR",
                       Either(Float(),
                              Undefined,
                              output=False,
                              optional=True,
                              desc=desc_TR))
        self.TR = Undefined

        self.add_trait("TA",
                       Either(Float(),
                              Undefined,
                              output=False,
                              optional=True,
                              desc=desc_TA))
        self.TA = Undefined

        self.add_trait("slice_order",
                       Either(List(Either(Int(),
                                          Float())),
                              Undefined,
                              output=False,
                              optional=True,
                              desc=desc_slice_order))
        self.slice_order = Undefined

        self.add_trait("ref_slice",
                       Either(Int(),
                              Float(),
                              Undefined,
                              output=False,
                              optional=True,
                              desc=desc_ref_slice))
        self.ref_slice = Undefined

        self.add_trait("out_prefix",
                       String('a',
                              output=False,
                              optional=True,
                              desc=desc_out_prefix))

        # Output traits
        self.add_trait("timed_files",
                       OutputMultiPath(File(),
                                       output=True,
                                       desc=desc_timed_files))

        self.init_default_traits()
        self.init_process('nipype.interfaces.spm.SliceTiming')

    def _get_database_value(self):
        """sets default values for certain parameters.

        If the following parameters are undefined, this method attempts to give
        them a value from the database or by making a calculation from other
        known parameters:
        - num_slices (from database)
        - TR (from database)
        - TA (from TR-(TR/num_slices))
        - slice_order (from the "acquisition" parameter value)
        - ref_slice (from the "acquisition" parameter value)

        :returns: True if no problem is detected, otherwise False

        """
        # FIXME: If several data are entered in in_files parameter, we use the
        # first element to determine, from the database, few calculation
        # parameters. These calculation parameters must be identical for all
        # data in in_files. Currently there is no verification to see if it's
        # true (e.g. do all data have the same slice number?, etc.).
        result = True
        complete_path = self.in_files[0]
        file_position = (complete_path.find(self.project.getName())
                            + len(self.project.getName()) + 1)
        database_filename = complete_path[file_position:]
        temp = None
        
        if ('Dataset dimensions (Count, X,Y,Z,T...)' in
                     self.project.session.get_fields_names(COLLECTION_CURRENT)):
            temp = self.project.session.get_value(COLLECTION_CURRENT,
                                                  database_filename,
                                                  ('Dataset dimensions '
                                                   '(Count, X,Y,Z,T...)'))[3]

        if self.num_slices == Undefined:

            if temp:
                self.num_slices = temp

            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle("User_processes_ECdev - "
                                   "SliceTiming Error!")
                msg.setText("Warning: the value of the num_slices parameter "
                            "was not found in the database. "
                            "Please check the data.")
                msg.setStandardButtons(QMessageBox.Close)
                msg.buttonClicked.connect(msg.close)
                msg.exec()
                result = False

        else:

            if temp:

                if self.num_slices != temp:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Warning)
                    msg.setWindowTitle("User_processes_ECdev - "
                                       "SliceTiming Error!")
                    msg.setText("Warning: The value for the num_slices "
                                "parameter does not match the value in the "
                                "database. Please check the data.")
                    msg.setStandardButtons(QMessageBox.Close)
                    msg.buttonClicked.connect(msg.close)
                    msg.exec()
                    result = False

        temp = None

        if ('RepetitionTime' in
                     self.project.session.get_fields_names(COLLECTION_CURRENT)):
            temp = self.project.session.get_value(COLLECTION_CURRENT,
                                                  database_filename,
                                                  'RepetitionTime')[0]/1000

        if self.TR == Undefined:

            if temp:
                self.TR = temp

            else:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle("User_processes_ECdev - "
                                   "SliceTiming Error!")
                msg.setText("Warning: the value of the TR parameter "
                            "was not found in the database. "
                            "Please check the data.")
                msg.setStandardButtons(QMessageBox.Close)
                msg.buttonClicked.connect(msg.close)
                msg.exec()
                result =  False

        else:

            if temp:

                if self.TR != temp:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Warning)
                    msg.setWindowTitle("User_processes_ECdev - "
                                       "SliceTiming Error!")
                    msg.setText("Warning: The value for the TR "
                                "parameter does not match the value in the "
                                "database. Please check the data.")
                    msg.setStandardButtons(QMessageBox.Close)
                    msg.buttonClicked.connect(msg.close)
                    msg.exec()
                    result = False

        # FIXME: in SPM, TA can be egal to 0 (If slice_order and
        # ref_slice are entered in milliseconds, TA will not be used and
        # can be set to 0). Here we consider that slice_order and
        # ref_slice paramters are only indices !

        if self.TA != 0:

            if self.TR != Undefined and self.num_slices != Undefined:
       
                if self.TA != Undefined:
             
                    if self.TA != self.TR - (self.TR / self.num_slices):
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Warning)
                        msg.setWindowTitle("User_processes_ECdev - "
                                           "SliceTiming Error!")
                        msg.setText("Warning: the value of the TA parameter "
                                    "doesn't exactly match the other "
                                    "parameters. Please check the data.")
                        msg.setStandardButtons(QMessageBox.Close)
                        msg.buttonClicked.connect(msg.close)
                        msg.exec()
                        result = False

                else:
                    self.TA = self.TR - (self.TR / self.num_slices)

            elif self.TA == Undefined:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle("User_processes_ECdev - "
                                   "SliceTiming Error!")
                msg.setText("The value of TR or num_slices parameters was not "
                            "found in the database, then the TA parameter "
                            "could not be calculated automatically. Please "
                            "check the data.")
                msg.setStandardButtons(QMessageBox.Close)
                msg.buttonClicked.connect(msg.close)
                msg.exec()
                result = False

        if self.slice_order != Undefined:

            if len(self.slice_order) != self.num_slices:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle("User_processes_ECdev - SliceTiming Error!")
                msg.setText("Warning: the slice_order parameter doesn't "
                            "match the num_slices parameter. Please check the "
                            "data.")
                msg.setStandardButtons(QMessageBox.Close)
                msg.buttonClicked.connect(msg.close)
                msg.exec()
                result = False

        elif self.num_slices != Undefined:
            
            if self.acquisition == "sequential ascending":
                self.slice_order = list(range(1, self.num_slices+1))

            if self.acquisition == "sequential descending":
                self.slice_order = list(range(self.num_slices, 0, -1))

            if self.acquisition == "interleaved (middle-top)":
                self.slice_order = []

                for (a, b) in (itertools.
                           zip_longest)(list(range(math.ceil(self.num_slices/2),
                                                   0,
                                                   -1)),
                                        list(range(math.ceil(self.num_slices),
                                                   math.ceil(self.num_slices/2),
                                                   -1))):

                    if (self.num_slices % 2) == 0:
                        self.slice_order.append(b)
                        self.slice_order.append(a)

                    else:
                        self.slice_order.append(a)
                        self.slice_order.append(b)

                if None in self.slice_order: self.slice_order.remove(None)

            if self.acquisition == "interleaved (bottom-up)":
                self.slice_order = (list(range(1, self.num_slices+1, 2))
                                    +list(range(2, self.num_slices+1, 2)))

            if self.acquisition == "interleaved (top-down)":
                self.slice_order = (list(range(self.num_slices, 0, -2))
                                    +list(range(self.num_slices-1, 0, -2)))

            if self.TA == 0 and self.TR != Undefined:
                self.slice_order = [(i-1) * (self.TR/self.num_slices) * 1000
                                    for i in self.slice_order]

        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle("User_processes_ECdev - SliceTiming Error!")
            msg.setText("Warning: As the num_slices parameter was not found in "
                        "the database, it was not possible to determine the "
                        "slice_order parameter automatically. Please check the "
                        "data.")
            msg.setStandardButtons(QMessageBox.Close)
            msg.buttonClicked.connect(msg.close)
            msg.exec()
            result = False

        if self.ref_slice == Undefined:

            if (self.acquisition == "sequential ascending" or
                                   self.acquisition == "sequential descending"):
                
                if self.num_slices != Undefined:
                    self.ref_slice = math.ceil(self.num_slices/2)

                else:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Warning)
                    msg.setWindowTitle("User_processes_ECdev - SliceTiming "
                                       "Error!")
                    msg.setText("Warning: It was not possible to determine the "
                                "ref_slice parameter automatically. Please "
                                "check the data.")
                    msg.setStandardButtons(QMessageBox.Close)
                    msg.buttonClicked.connect(msg.close)
                    msg.exec()
                    result = False

            else:
                self.ref_slice = 1

            if (self.TA == 0 and
                    self.TR != Undefined and
                    self.num_slices != Undefined):
                self.ref_slice = ((self.ref_slice-1) *
                                                 self.TR/self.num_slices * 1000)

        return result

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
        super(SliceTiming, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_files and self.in_files != [Undefined]:
            self.process.in_files = self.in_files
            # Definition of the necessary inputs for the spm fonction
            # acquired from the databse
            res = self._get_database_value()

            if not res:
                return self.make_initResult()

            if (self.num_slices != Undefined and
                        self.TR != Undefined and
                        self.TA != Undefined and
                        self.slice_order != Undefined and
                        self.ref_slice != Undefined):

                if self.out_prefix:
                    self.process.out_prefix = self.out_prefix

                # The management of self.process.output_directory could be delegated
                # to the populse_mia.user_interface.pipeline_manager.process_mia
                # module. We can't do it at the moment because the
                # sync_process_output_traits() of the capsul/process/nipype_process
                # module raises an exception in nipype if the mandatory parameters
                # are not yet defined!
                if self.output_directory:
                    self.process.output_directory = self.output_directory

                else:
                    print('No output_directory was found...!\n')

                self.outputs['timed_files'] = self.process._timecorrected_files

        if self.outputs:

            for key, val in self.outputs.items():

                if key == "timed_files":
                    
                    if not isinstance(val, list):
                        val = [val]

                    for in_val, out_val in zip(self.in_files, val):
                        _, fileOval = os.path.split(out_val)
                        _, fileIval = os.path.split(in_val)

                        if self.out_prefix:
                            fileOvalNoPref = fileOval[len(self.out_prefix):]
                        
                        else:
                           fileOvalNoPref = fileOval[len('s'):]

                        if fileOvalNoPref == fileIval:
                            self.inheritance_dict[out_val] = in_val

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(SliceTiming, self).run_process_mia()

        # in_files parameter are normally already in absolute path format.
        # So, the 3 next lines can be see as "in case of"...
        for idx, element in enumerate(self.in_files):
            full_path = os.path.abspath(element)
            self.in_files[idx] = full_path

        self.process.in_files = self.in_files
        self.process.num_slices = self.num_slices
        self.process.time_repetition = self.TR
        self.process.time_acquisition = self.TA
        self.process.slice_order = self.slice_order
        self.process.ref_slice = self.ref_slice

        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})


class Smooth(ProcessMIA):
    """
    *3D Gaussian smoothing of image volumes*

    Please, see the complete documention for the `Smooth brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/spm/Smooth.html>`_

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
        self.requirement = ['spm', 'nipype']

        # Inputs description
        in_files_desc = ('List of files to smooth. A list of items which are '
                         'an existing, uncompressed file '
                         '(valid extensions: [.img, .nii, .hdr]).')
        fwhm_desc = ('Full-width at half maximum (FWHM) of the Gaussian '
                     'smoothing kernel in mm (a float or a list of 3 items '
                     'which are a float).')
        data_type_desc = ('Data type of the output images '
                          '(an integer [int or long]).')
        implicit_masking_desc = ('A mask implied by a particular voxel value '
                                 '(a boolean).')
        out_prefix_desc = ('Specify the string to be prepended to the '
                           'filenames of the smoothed image file(s) '
                           '(a string).')

        # Outputs description
        smoothed_files_desc = ('The smoothed files (a pathlike object or '
                               'string representing a file, or a list of '
                               'pathlike objects or strings representing a '
                               'file).')

        # Input traits
        self.add_trait("in_files",
                       InputMultiPath(Either(ImageFileSPM(),
                                             Undefined),
                                      copyfile=False,
                                      output=False,
                                      optional=False,
                                      desc=in_files_desc))

        self.add_trait("fwhm",
                       Either(Float(),
                              List(Float(),
                                   minlen=3, maxlen=3),
                              default=[6.0,6.0,6.0],
                              output=False,
                              optional=True,
                              desc=fwhm_desc))

        self.add_trait("data_type",
                       Int(output=False,
                           optional=True,
                           desc=data_type_desc))
        
        self.add_trait("implicit_masking",
                       Bool(output=False,
                            optional=True,
                            desc=implicit_masking_desc))

        self.add_trait("out_prefix",
                       String('s',
                              output=False,
                              optional=True,
                              desc=out_prefix_desc))

        # Output traits 
        self.add_trait("smoothed_files",
                       OutputMultiPath(File(),
                                       output=True,
                                       desc=smoothed_files_desc))

        self.init_default_traits()
        self.init_process('nipype.interfaces.spm.Smooth')

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
        if self.in_files and self.in_files != [Undefined]:
            self.process.in_files = self.in_files

            if self.out_prefix:
                self.process.out_prefix = self.out_prefix

            # The management of self.process.output_directory could be delegated
            # to the populse_mia.user_interface.pipeline_manager.process_mia
            # module. We can't do it at the moment because the
            # sync_process_output_traits() of the capsul/process/nipype_process
            # module raises an exception in nipype if the mandatory parameters
            # are not yet defined!
            if self.output_directory:
                self.process.output_directory = self.output_directory

            else:
                print('No output_directory was found...!\n')

            self.outputs['smoothed_files'] = self.process._smoothed_files

        if self.outputs:
        
            for key, val in self.outputs.items():

                if key == "smoothed_files":

                    if not isinstance(val, list):
                        val = [val]

                    for in_val, out_val in zip(self.in_files, val):
                        _, fileOval = os.path.split(out_val)
                        _, fileIval = os.path.split(in_val)
                        
                        if self.out_prefix:
                            fileOvalNoPref = fileOval[len(self.out_prefix):]
                        
                        else:
                           fileOvalNoPref = fileOval[len('s'):]

                        if fileOvalNoPref == fileIval:
                            self.inheritance_dict[out_val] = in_val

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()
    
    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Smooth, self).run_process_mia()

        # in_files parameter are normally already in absolute path format.
        # So, the 3 next lines can be see as "in case of"...
        for idx, element in enumerate(self.in_files):
            full_path = os.path.abspath(element)
            self.in_files[idx] = full_path

        self.process.in_files = self.in_files
        self.process.fwhm = self.fwhm
        self.process.data_type = self.data_type
        self.process.implicit_masking = self.implicit_masking
            
        if self.out_prefix:
            self.process.out_prefix = self.out_prefix

        return self.process.run(configuration_dict={})
