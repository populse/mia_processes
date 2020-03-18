# -*- coding: utf-8 -*- #

"""The spm preprocess library of the mia_processes package.

The purpose of this module is to customise the main spm preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - Coregister
        - NewSegment
        - Normalize12
        - Realign
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
from .nipype_extension import NewSegmentMia

# populse_mia import
from populse_mia.software_properties import Config

# nipype import
from nipype.interfaces import spm
from nipype.interfaces.base import OutputMultiPath, InputMultiPath, File, traits
from nipype.interfaces.spm.base import ImageFileSPM

# Other import
import os
from traits.api import Undefined, Float
from pathlib import Path

class Coregister(Process_Mia):
    """
*Align together scans of different modalities*

-------------

.. [#label] Syntax: mia_processes/nipype Coregister <=> SPM12 Coregister.

    Useful links:
    `SPM12 Coregister <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=39>`_,
    `nipype Coregister <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.preprocess.html#coregister>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Coregister, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['matlab', 'spm']

        # Inputs description
        target_desc = '''*<=> ref* [#label]_

The reference file (remains stationary) while the source image is moved to match it. An existing, uncompressed file (valid extensions: [.img, .nii, .hdr]).

::

  ex. /home/ArthurBlair/data/downloaded_data/meanFunc.nii
'''
        source_desc = '''*<=> source* [#label]_

The image that is jiggled about to best match the target image. A list of items which are an existing, uncompressed file (valid extensions: [.img, .nii, .hdr]).

::

  ex. ['/home/ArthurBlair/data/raw_data/Anat.nii']
'''
        apply_to_files_desc = '''*<=> other* [#label]_

These are any images that need to remain in alignment with the source image (a list of items which are an existing file name) (valid extensions: [.img, .nii, .hdr]).

::

  ex. ['/home/ArthurBlair/data/raw_data/Func.nii']
'''
        jobtype_desc = '''[#label]_

One of 'estwrite' or 'estimate' or 'write'. If 'estimate' is selected, the registration parameters are stored in the headers of the 'source' and the 'apply_to_files' images. If 'write' is selected, the resliced images are named the same as the originals except that they are prefixed by out_prefix. if 'estwrite' is selected, the described procedures for 'estimate' and 'write' are performed and the output parameter (coregistered_source and/or coregistered_files) contains the resliced images and the one whose header has been rewritten. If it is necessary to choose one or the other for a subsequent calculation in a pipeline, the Auto_Filter_List brick (mia_processes library) can be used.

::

    ex. estimate
'''
        cost_function_desc = '''*<=> eoptions.cost_fun* [#label]_

One of 'mi' or 'nmi' or 'ecc' or 'ncc'. Registration involves finding parameters that either maximise or minimise some objective function. For inter-modal registration,  use 'Mutual Information', 'Normalised Mutual Information' or 'Entropy Correlation Coefficient'. For within modality, you could also use Normalised Cross Correlation.

  - 'mi': Mutual Information
  - 'nmi': Normalised Mutual Information
  - 'ecc': Entropy Correlation Coefficient
  - 'ncc': Normalised Cross Correlation

::

  ex. nmi
'''
        separation_desc = '''*<=> eoptions.sep* [#label]_

A list of items which are a float. The average distance between sampled points (in mm). Can be a vector to allow a coarse registration followed by increasingly fine ones.

::

  ex. [4, 2]
'''
        tolerance_desc = '''*<=> eoptions.tol* [#label]_

A list of 12 items which are a float. The acceptable tolerance for each of 12 params. Iterations stop when differences between
successive estimates are less than the required tolerance.

::

  ex. [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001]
'''
        fwhm_desc = '''*<=> eoptions.fwhm* [#label]_

A list of 2 items which are a float. Kernel of gaussian smooth to apply to the 256*256 joint histogram.

::

  ex. [7, 7]
'''
        write_interp_desc = '''*<=> roptions.interp* [#label]_

The method by which the images are sampled when being written in a different space. Nearest neighbour is fastest, but not recommended for image realignment. Trilinear Interpolation is probably OK for PET, or realigned and re-sliced fMRI, but not so suitable for fMRI with subject movemen because higher degree interpolation generally gives better results. Although higher degree methods provide better interpolation, but they are slower because they use more neighbouring voxels. (0 <= a long integer <= 7). Voxel sizes must all be identical and isotropic.

  - 0: Nearest neighbour
  - 1: Trilinear
  - 2: 2nd Degree B-Spline
  - 3: 3rd Degree B-Spline

  …

  - 7: 7th Degree B-Spline

::

  ex. 4
'''
        write_wrap_desc = '''*<=> roptions.wrap* [#label]_

Check if interpolation should wrap in [x,y,z] (a list of 3 items which are integer int or long). For example, in MRI scans, the images wrap around in the phase encode direction, so the subject’s nose may poke into the back of the subject’s head. These are typically:

  - No wrapping [0, 0, 0]: for PET or images that have already been spatially transformed (Also the recommended option if you are not really sure)
  - Wrap in Y [0, 1, 0], for (un-resliced) MRI where phase encoding is in the Y direction (voxel space)

::

    ex. [0 0 0]
'''
        write_mask_desc = '''*<=> roptions.mask* [#label]_

Mask output image (a boolean). Because of subject motion, different images are likely to have different patterns of zeros from where it was not possible to sample data. With masking enabled, the program searches through the whole time series looking for voxels which need to be sampled from outside the original images. Where this occurs, that voxel is set to zero for the whole set of images.

::

  ex. False
'''
        out_prefix_desc = '''*<=> roptions.prefix* [#label]_

Specify the string to be prepended to the filenames of the coregisterd image file(s).

::

  ex. r, capsul/nipype default value
'''
        
        # Outputs description
        coregistered_source_desc = '''Coregistered source files, corresponding to ‘source’ images (a list of items which are an existing file name).

::

  ex. /home/ArthurBlair/data/raw_data/meanFunc.nii
'''
        coregistered_files_desc = '''A list of items which are an existing file name. Coregistered other files, corresponding to 'apply_to_files' images.

::

  ex. /home/ArthurBlair/data/raw_data/Func.nii
'''

        # Inputs traits 
        self.add_trait("target",
                       ImageFileSPM(copyfile=False,
                                    output=False,
                                    optional=False,
                                    desc=target_desc))
        
        self.add_trait("source",
                       InputMultiPath(traits.Either(ImageFileSPM(),
                                                    Undefined),
                                      value=[Undefined],
                                      copyfile=True,
                                      output=False,
                                      optional=False,
                                      desc=source_desc))
        
        self.add_trait("apply_to_files",
                       InputMultiPath(traits.Either(File(),
                                                    Undefined),
                                      value=[Undefined],
                                      copyfile=True,
                                      output=False,
                                      optional=True,
                                      desc=apply_to_files_desc))

        self.add_trait("jobtype",
                       traits.Enum('estimate',
                                   'estwrite',
                                   'write',
                                   output=False,
                                   optional=True,
                                   desc=jobtype_desc))

        self.add_trait("cost_function",
                       traits.Enum('nmi',
                                   'mi',
                                   'ecc',
                                   'ncc',
                                   output=False,
                                   optional=True,
                                   desc=cost_function_desc))

        self.add_trait("separation",
                       traits.List(value=[4.0, 2.0],
                                   trait=traits.Range(low=0.0, high=None),
                                   minlen=1,
                                   maxlen=32,
                                   output=False,
                                   optional=True,
                                   desc=separation_desc))

        self.add_trait("tolerance",
                       traits.List(value=[.02, .02, .02, 0.001, 0.001, 0.001,
                                          .01, .01, .01, 0.001, 0.001, 0.001],
                                   trait=traits.Range(low=0.0, high=None),
                                   minlen=12,
                                   maxlen=12,
                                   output=False,
                                   optional=True,
                                   desc=tolerance_desc))

        self.add_trait("fwhm",
                       traits.List(value=[7.0, 7.0],
                                   trait=traits.Range(low=0.0, high=None),
                                   minlen=2,
                                   maxlen=2,
                                   output=False,
                                   optional=True,
                                   desc=fwhm_desc))

        self.add_trait("write_interp",
                       traits.Range(value=4,
                                    low=0,
                                    high=7,
                                    output=False,
                                    optional=True,
                                    desc=write_interp_desc))

        self.add_trait("write_wrap",
                       traits.List(value=[0, 0, 0],
                                   trait=traits.Range(low=0, high=1),
                                   minlen=3,
                                   maxlen=3,
                                   output=False,
                                   optional=True,
                                   desc=write_wrap_desc))

        self.add_trait("write_mask",
                       traits.Bool(default_value=False,
                                   output=False,
                                   optional=True,
                                   desc=write_mask_desc))

        self.add_trait("out_prefix",
                       traits.String('r',
                                     output=False,
                                     optional=True,
                                     desc=out_prefix_desc))

        # Output traits
        self.add_trait("coregistered_source",
                       OutputMultiPath(File(),
                                       output=True,
                                       desc=coregistered_source_desc))

        self.add_trait("coregistered_files",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=coregistered_files_desc))

        self.process = spm.Coregister()
        self.change_dir = True

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return 
        self.make_initResult() object.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Coregister, self).list_outputs()
        
        # Outputs definition and tags inheritance (optional)
        if self.outputs:
            self.outputs = {}

        if self.inheritance_dict:
            self.inheritance_dict = {}

        if (self.target and self.source and
              self.source != [Undefined] and self.jobtype):
            self.process.inputs.target = self.target
            self.process.inputs.source = self.source
            self.process.inputs.jobtype = self.jobtype

            if self.apply_to_files and self.apply_to_files != [Undefined]:
                self.process.inputs.apply_to_files = self.apply_to_files

            else:
                self.apply_to_files = Undefined

            if self.out_prefix:
                self.process.inputs.out_prefix = self.out_prefix

            self.outputs = self.process._list_outputs()

        if self.outputs:
            outputs_coregsource = None
            outputs_coregfiles = None
        
            for key, values in self.outputs.items():

                if key == "coregistered_source":
                    outputs_coregsource = values.copy()
                    
                    for fullname in values:
                        path, filename_out = os.path.split(fullname)

                        if not self.jobtype == "estimate":

                            if self.out_prefix:
                                filename_in = filename_out[
                                                  len(self.out_prefix):]
                        
                            else:
                                filename_in = filename_out[len('r'):]

                        else:
                            filename_in = filename_out
    
                        if os.path.join(path, filename_in) in self.source:
                            self.inheritance_dict[
                                           fullname] = os.path.join(path,
                                                                    filename_in)

                        if self.jobtype == "estwrite":

                            if isinstance(outputs_coregsource, list):
                                outputs_coregsource.append(
                                                os.path.join(path, filename_in))
                                self.inheritance_dict[os.path.join(path,
                                                                   filename_in)
                                                     ] = os.path.join(path,
                                                                    filename_in)
                                
                if (key == "coregistered_files" and
                    not values in ["<undefined>", traits.Undefined]):
                    outputs_coregfiles = values.copy()
                    
                    for fullname in values:
                        path, filename_out = os.path.split(fullname)

                        if not self.jobtype == "estimate":

                            if self.out_prefix:
                                filename_in = filename_out[
                                                  len(self.out_prefix):]

                            else:
                                filename_in = filename_out[len('r'):]

                        else:
                            filename_in = filename_out

                        if os.path.join(path,
                                        filename_in) in self.apply_to_files:
                            self.inheritance_dict[
                                           fullname] = os.path.join(path,
                                                                    filename_in)

                        if self.jobtype == "estwrite":

                            if isinstance(outputs_coregfiles, list):
                                outputs_coregfiles.append(
                                               os.path.join(path, filename_in))
                                self.inheritance_dict[os.path.join(path,
                                                                   filename_in)
                                                     ] = os.path.join(path,
                                                                    filename_in)

            if outputs_coregsource:
                self.outputs["coregistered_source"] = outputs_coregsource

            if outputs_coregfiles:
                self.outputs["coregistered_files"] = outputs_coregfiles

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Coregister, self).run_process_mia()
        self.process.inputs.target = self.target
        self.process.inputs.source = self.source
        self.process.inputs.apply_to_files = self.apply_to_files
        self.process.inputs.jobtype = self.jobtype
        self.process.inputs.cost_function = self.cost_function
        self.process.inputs.fwhm = self.fwhm
        self.process.inputs.separation = self.separation
        self.process.inputs.tolerance = self.tolerance
        self.process.inputs.out_prefix = self.out_prefix
        self.process.run()


class NewSegment(Process_Mia):
    """
*Segmentation: Segments,  bias  corrects  and  spatially normalises - all in the same model*

-------------

.. [#label] Syntax: mia_processes/nipype NewSegment <=> SPM12 Segment.

NOTE:
    - This interface currently supports single channel input only.
    - The warp.mrf, warp.cleanup and warp.fwhm, from  SPM12, are not used in this brick.

Useful links:
`SPM12 Segment <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=45>`_,
`nipype NewSegment <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.preprocess.html#newsegment>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(NewSegment, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['matlab', 'spm']

        # Inputs description
        
        channel_files_desc = ('Path of the scans for processing, valid '
                              'extensions: [.img, .nii, .hdr]. A list with one '
                              'string element corresponding to an existing '
                              'path file.')
        channel_info_desc = '''*<=> (channel.biasreg, channel.biasfwhm, (channel.write))* [#label]_

A tuple (consisting of a float, a float and a tuple consisting of a boolean, a boolean) with the following fields:

- bias reguralisation (a float between 0 and 10)
    The goal is to model, by different tissue classes, the intensity variations that arise due to different tissues, while model, with a bias field, those that occur because of the bias artifact due to the physics of MRI imaging. If the data have very little intensity non-uniformity artifact, then bias control should be increased. This effectively tells the algorithm that there is very little bias in the data, so it doesn't try to model it.

      | \- 0 No regularisation
      | \- 0.00001 extremely light regularisation
      | \-  ...
      | \- 1 very heavy regularisation
      | \- 10 extremely heavy regularisation

- bias FWHM (a float between 20 and infinity)
    Full Width at Half Maximum of Gaussian smoothness of bias. Smoother bias fields need fewer parameters to describe them. This means that the algorithm is faster for smoother intensity non-uniformities (e.g. 150 mm cutoff gives faster results than 20 mm cutoff).

- which maps to save (a tuple of two boolean values; (Field, Corrected))
    To save the estimated bias field or/and the bias corrected version of the processed image.

      | \- (False, False) save Nothing
      | \- (False, True) save bias corrected image only
      | \- (True, False) save estimated bias field only
      | \- (True, True) save estimated bias field and bias corrected image

::

  ex. (0.0001, 60, (False, True))
'''
        tissues_desc = '''*<=> [((tissue(i).tpm), tissue(i).ngaus, (tissue(i).native), (tissue(i).warped)), ((tissue(i+1).tpm), tissue(i+1).ngaus, (tissue(i+1).native), (tissue(i+1).warped)), ...]* [#label]_

A list of tuples (one per tissue, i from 1 to 6) with parameter values for each tissue types. Typically, the order of tissues is grey matter (i=1), white matter (i=2), CSF (i=3), bone (i=4), soft tissue (i=5) and air/background (i=6), if using tpm/TPM.nii from spm12.

Each tuple consists of the following fields:

    (tissue probability map (4D), 1-based index to frame),  number of gaussians, (which maps to save; Native, DARTEL), (which maps to save; Unmodulated, Modulated)

        * tissue probability map <=> tissue(i).tpm with i in (1, 2, 3, 4, 5, 6])
            The tissue probability image [.img, .nii, .hdr].

        * 1-based index to frame
            Index for the 4th dimension of the tissue probability map and then tissue type selection.

              | \- 1 to 6

        * number of gaussians <=> tissue(i).ngaus
            Typical numbers of Gaussians could be 2 for GM, WM, CSF, 3 for bone, 4 for other soft tissues and
            2 for air/background.

              | \- 1, 2, 3, 4, 5, 6 , 7, 8, inf -Non parametric-

        * which maps to save; Native, DARTEL <=> tissue(i).native
            To produce a tissue class image that is in alignment with the original (ci) or that can be used with
            the Dartel toobox (rci).

              | \- (False, False) Save Nothing
              | \- (True, False) save native only
              | \- (False, True ) save DARTEL only
              | \- etc.

        * which maps to save [Unmodulated, Modulated] <=> tissue(i).warped
            To produces spatially normalised versions of the tissue class, with (mcwi) and without (wci)
            modulation.

              | \ - (False, False) Save Nothing
              | \ - (True, False) save unmodulated only
              | \ - (False, True ) save modulated only
              | \ - etc.

::

  ex.  [(('/home/ArthurBlair/spm12/tpm/TPM.nii', 1), 2, (True, False), (False, False)),
        (('/home/ArthurBlair/spm12/tpm/TPM.nii', 2), 2, (True, False), (False, False)),
        (('/home/ArthurBlair/spm12/tpm/TPM.nii', 3), 2, (True, False), (False, False)),
        (('/home/ArthurBlair/spm12/tpm/TPM.nii', 4), 3, (True, False), (False, False)),
        (('/home/ArthurBlair/spm12/tpm/TPM.nii', 5), 4, (True, False), (False, False)),
        (('/home/ArthurBlair/spm12/tpm/TPM.nii', 6), 2, (True, False), (False, False))]
'''
        warping_regularization_desc = '''*<=> warp.reg* [#label]_

The measure of the roughness of the deformations for registration. Involve the sum of 5 elements (a float or list of floats; the latter is required by SPM12).

::

  ex. [0, 0.001, 0.5, 0.05, 0.2]
'''
        affine_regularization_desc = '''*<=> warp.affreg* [#label]_

Standard space for affine registration ('mni' or 'eastern' or 'subj' or 'none').

::

    ex. mni
'''

        sampling_distance_desc = '''*<=> warp.samp* [#label]_

Approximate distance between sampled points when estimating the model parameters (a float).

::

  ex. 3
'''
        write_deformation_fields = '''*<=> warp.write* [#label]_

Deformation fields can be saved to disk, and used by the deformation utility (a list of 2 booleans for which deformation fields to write; Inverse, Forward).

    | \- [False, False] Save nothing
    | \- [True, False] save Inverse only
    | \- [False, True] save Forward only
    | \- etc.

::

  ex. [False, True]
'''

        # Outputs description
        bias_corrected_images_desc = '''The bias corrected images (a list of items which are a pathlike object or string representing an existing file).

::

  ex. /home/ArthurBlair/data/raw_data/mAnat.nii
'''
        bias_field_images_desc = '''The estimated bias field (a list of items which are a pathlike object or string representing an existing file).

::

  ex. <undefined>
'''
        native_class_images_desc = '''Native space probability maps (a list of items which are a list of items which are a pathlike object or string representing an existing file).

::

  ex. [['/home/ArthurBlair/data/raw_data/c1Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/c2Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/c3Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/c4Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/c5Anat.nii']]
'''
        dartel_input_images_desc = '''"Imported" class images into a form that can be used with the Dartel toolbox (a list of items which are a list of items which are a pathlike object or string representing an existing file).

::

  ex. [['/home/ArthurBlair/data/raw_data/rc1Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/rc2Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/rc3Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/rc4Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/rc5Anat.nii']]
'''
        modulated_class_images_desc = '''Modulated and normalised class images (a list of items which are a list of items which are a pathlike object or string representing an existing file).

::

  ex. [['/home/ArthurBlair/data/raw_data/mwc1Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/mwc2Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/mwc3Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/mwc4Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/mwc5Anat.nii']]
'''
        normalized_class_images_desc = '''Normalised class images, without modulation (a list of items which are a list of items which are a pathlike object or string representing an existing file).

::

  ex. [['/home/ArthurBlair/data/raw_data/wc1Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/wc2Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/wc3Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/wc4Anat.nii'],
       ['/home/ArthurBlair/data/raw_data/wc5Anat.nii']]
'''
        inverse_deformation_field_desc = '''Inverse deformation field. Could be used for spatially normalising surface files as GIFTI (a list of items which are a pathlike object or string representing an existing file).

::

  ex. /home/ArthurBlair/data/raw_data/iy_Anat.nii
'''
        forward_deformation_field_desc = '''Forward deformation field. Could be used for spatially normalising images to MNI space (a list of items which are a pathlike object or string representing an existing file).

::

  ex. /home/ArthurBlair/data/raw_data/y_Anat.nii
'''
        transformation_mat_desc = '''Normalisation transformation (a list of items which are a pathlike object or string representing an existing file).

::

  ex. /home/ArthurBlair/data/raw_data/Anat_seg8.mat
'''

        # Tissues parameter definition
        config = Config()
        resources_path = os.path.join(config.get_mia_path(), 'resources')
        tpm_path = os.path.join(resources_path, 'spm12', 'tpm', 'TPM.nii')
        
        if not Path(tpm_path).exists():
            print('\n The {} file seems to not exists ...'.format(tpm_path))

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
                       traits.Tuple(traits.Range(value=0.0001,
                                                 low=0.,
                                                 high=10.),
                                    traits.Range(value=60.,
                                                 low=20.,
                                                 high=None),
                                    traits.Tuple(traits.Bool(False),
                                                 traits.Bool(True)),
                                    output=False,
                                    optional=True,
                                    desc=channel_info_desc))

        self.add_trait("tissues",
                       traits.List(
                           traits.Tuple(
                               traits.Tuple(ImageFileSPM(exists=True),
                                            traits.Int()),
                               traits.Int(),
                               traits.Tuple(traits.Bool, traits.Bool),
                               traits.Tuple(traits.Bool, traits.Bool)),
                           value=tissues_list,
                           output=False,
                           optional=True,
                           desc=tissues_desc))

        self.add_trait("warping_regularization",
                       traits.Either(traits.List(traits.Float(),
                                                 minlen=5,
                                                 maxlen=5),
                                     traits.Float(),
                                     default=[0, 0.001, 0.5, 0.05, 0.2],
                                     output=False,
                                     optional=True,
                                     desc=warping_regularization_desc))

        self.add_trait("affine_regularization",
                       traits.Enum('mni',
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
                       traits.List(traits.Bool(),
                                   minlen=2,
                                   maxlen=2,
                                   value=[False, True],
                                   output=False,
                                   optional=True,
                                   desc=write_deformation_fields))

        # Output traits
        self.add_trait("bias_corrected_images",
                       OutputMultiPath(output=True,
                                       optional=True,
                                       desc=bias_corrected_images_desc))

        self.add_trait("bias_field_images",
                       OutputMultiPath(output=True,
                                       optional=True,
                                       desc=bias_field_images_desc))

        self.add_trait("native_class_images",
                       traits.List(traits.List(File()),
                                   output=True,
                                   optional=True,
                                   desc=native_class_images_desc))

        self.add_trait("dartel_input_images",
                       traits.List(traits.List(File()),
                                   output=True,
                                   optional=True,
                                   desc=dartel_input_images_desc))

        self.add_trait("modulated_class_images",
                       traits.List(traits.List(File()),
                                   output=True,
                                   optional=True,
                                   desc=modulated_class_images_desc))

        self.add_trait("normalized_class_images",
                       traits.List(traits.List(File()),
                                   output=True,
                                   optional=True,
                                   desc=normalized_class_images_desc))

        self.add_trait("inverse_deformation_field",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=inverse_deformation_field_desc))

        self.add_trait("forward_deformation_field",
                       ImageFileSPM(output=True,
                                    optional=True,
                                    desc=forward_deformation_field_desc))

        self.add_trait("transformation_mat",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=transformation_mat_desc))

        # process instanciation
        self.process = NewSegmentMia() # workaround to try to decrease the
                                       # instantiation time of the NewSegment
                                       # class
        #self.process = spm.NewSegment()
        self.change_dir = True

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return 
        self.make_initResult() object.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(NewSegment, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.outputs:
            self.outputs = {}

        if self.inheritance_dict:
            self.inheritance_dict = {}

        if self.channel_files:
            self.process.inputs.channel_files = self.channel_files
            
            if self.channel_info:
                self.process.inputs.channel_info = self.channel_info

            if self.write_deformation_fields:
                (self.process.inputs.
                       write_deformation_fields) = self.write_deformation_fields

            if self.tissues:
                    self.process.inputs.tissues = self.tissues
                    self.outputs = self.process._list_outputs()

        """
         When there is only one image at the input, the tags inheritance is
         directly managed in PipelineManagerTab.add_plug_value_to_database().
         However the inheritance of 2 output parameters is defined below, as
         an example
        """
        if self.outputs:
        
            for key, values in self.outputs.items():

                if key == "native_class_images":

                    if values[0]:
                        path, filename = os.path.split(values[0][0])
                        filename_without_prefix = filename[2:]

                        if (os.path.join(path,
                                         filename_without_prefix)
                                                         in self.channel_files): 

                            for fullname in values:
                                self.inheritance_dict[fullname[0]
                                                     ] = os.path.join(path,
                                                        filename_without_prefix)

                if key == "forward_deformation_field":

                    for fullname in values:
                        path, filename = os.path.split(fullname)
                        filename_without_prefix = filename[2:]

                        if (os.path.join(path,
                                     filename_without_prefix)
                                                         in self.channel_files):

                            self.inheritance_dict[fullname
                                                 ] = os.path.join(path,
                                                        filename_without_prefix)

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()
    
    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(NewSegment, self).run_process_mia()
        self.process.inputs.channel_files = self.channel_files
        self.process.inputs.channel_info = self.channel_info
        self.process.inputs.tissues = self.tissues
        self.process.inputs.warping_regularization = self.warping_regularization
        self.process.inputs.affine_regularization = self.affine_regularization
        self.process.inputs.sampling_distance = self.sampling_distance
        (self.process.inputs.
                       write_deformation_fields) = self.write_deformation_fields
        self.process.run()


class Normalize12(Process_Mia):
    """    
    *Computes the warp that best aligns the template (atlas) to the individual’s image*

    Please, see the complete documention for the `Normalize12 brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/html/documentation/preprocess/spm/Normalize12.html>`_

    """
    
    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Normalize12, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['matlab', 'spm']

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
                                   'dimension.')
        normalized_image_desc = ('Normalised file that needed to be aligned (a '
                                'list of items which are an existing file '
                                 'name).')
        normalized_files_desc = ('Normalised other files (a list of items '
                                 'which are an existing file name).')


        # Tpm parameter definition
        config = Config()
        resources_path = os.path.join(config.get_mia_path(), 'resources')
        tpm_path = os.path.join(resources_path, 'spm12', 'tpm', 'TPM.nii')
        
        if not Path(tpm_path).exists():
            print('\n The {} file seems to not exists ...'.format(tpm_path))
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
                       InputMultiPath(traits.Either(ImageFileSPM(),
                                                    traits.List(ImageFileSPM()),
                                                    Undefined),
                                      value=[Undefined],
                                      output=False,
                                      optional=True,
                                      desc=apply_to_files_desc))

        self.add_trait("jobtype",
                       traits.Enum("write",
                                   "est",
                                   "estwrite",
                                   output=False,
                                   optional=True,
                                   desc=jobtype_desc))

        self.add_trait("bias_regularization",
                       traits.Enum(0.0001,
                                   0,
                                   0.00001,
                                   0.001,
                                   0.01,
                                   0.1,
                                   1,
                                   10,
                                   output=False,
                                   optional=True,
                                   desc=bias_regularization_desc))

        self.add_trait("bias_fwhm",
                       traits.Enum(60,
                                   30,
                                   40,
                                   50,
                                   70,
                                   80,
                                   90,
                                   100,
                                   110,
                                   120,
                                   130,
                                   140,
                                   150,
                                   "Inf",
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
                       traits.Enum("mni",
                                   "size",
                                   "none",
                                   output=False,
                                   optional=True,
                                   desc=affine_regularization_type_desc))

        self.add_trait("warping_regularization",
                       traits.List(value=[0, 0.001, 0.5, 0.05, 0.2],
                                   trait=traits.Float(),
                                   minlen=5,
                                   maxlen=5,
                                   output=False,
                                   optional=True,
                                   desc=warping_regularization_desc))

        self.add_trait("smoothness",
                       traits.Float(0.,
                                   output=False,
                                   optional=True,
                                   desc=smoothness_desc))

        self.add_trait("sampling_distance",
                       traits.Float(3,
                                    output=False,
                                    optional=True,
                                    desc=sampling_distance_desc))

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

        # Outputs traits
        self.add_trait("deformation_field",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=deformation_field_desc))
                       
        self.add_trait("normalized_image",
                        OutputMultiPath(File(),
                                        output=True,
                                        optional=True,
                                        desc=normalized_image_desc))

        self.add_trait("normalized_files",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=normalized_files_desc))

        # process instanciation
        self.process = spm.Normalize12()
        self.change_dir = True

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return 
        self.make_initResult() object.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Normalize12, self).list_outputs()

        # Outputs definition
        if self.outputs:
            self.outputs = {}

        _flag = False
        self.process.inputs.jobtype = self.jobtype
        
        if ((self.apply_to_files) and (self.apply_to_files != [Undefined]) and
                         (self.deformation_file) and (self.jobtype == 'write')):
            self.process.inputs.apply_to_files = self.apply_to_files
            self.process.inputs.deformation_file = self.deformation_file
            _flag = True

            if self.image_to_align:
                self.image_to_align = Undefined

            if self.tpm:
                self.tpm = Undefined
                
        elif (self.image_to_align) and (self.tpm) and (self.jobtype == 'est'):
            self.process.inputs.image_to_align = self.image_to_align
            self.process.inputs.tpm = self.tpm
            _flag = True

            if self.apply_to_files and (self.apply_to_files != [Undefined]):
                self.apply_to_files = Undefined

            if self.deformation_file:
                self.deformation_file = Undefined
                
        elif ((self.image_to_align) and (self.tpm) and
                                                (self.jobtype == 'estwrite')):
            self.process.inputs.image_to_align = self.image_to_align
            self.process.inputs.tpm = self.tpm
            _flag = True

            if self.deformation_file:
                self.deformation_file = Undefined

            if (self.apply_to_files) and (self.apply_to_files != [Undefined]):
                self.process.inputs.apply_to_files = self.apply_to_files

        if _flag:
            self.outputs = self.process._list_outputs()

        #Tags inheritance (optional)
        if self.inheritance_dict:
            self.inheritance_dict = {}

        if self.outputs:

            for key, values in self.outputs.items():

                if (key == "normalized_files") and (values != Undefined):

                    for fullname in values:
                        path, filename = os.path.split(fullname)
                        filename_without_prefix = filename[len('w'):]

                        if (os.path.join(path,
                                         filename_without_prefix)
                                                        in self.apply_to_files):
                            self.inheritance_dict[fullname] = os.path.join(path,
                                                        filename_without_prefix)

                if (key == "deformation_field") and (values != Undefined):

                    path, filename = os.path.split(values)
                    filename_without_prefix = filename[len('y_'):]

                    if (os.path.join(path,
                                     filename_without_prefix)
                                                        in self.image_to_align):
                        self.inheritance_dict[values] = os.path.join(path,
                                                        filename_without_prefix)

                if (key == "normalized_image") and (values != Undefined):

                    path, filename = os.path.split(values)
                    filename_without_prefix = filename[len('w'):]
                    
                    if (os.path.join(path,
                                     filename_without_prefix)
                                                        in self.image_to_align):
                        self.inheritance_dict[values] = os.path.join(path,
                                                        filename_without_prefix)

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Normalize12, self).run_process_mia()

        if self.jobtype == 'write':
            self.process.inputs.write_bounding_box = self.write_bounding_box
            self.process.inputs.write_voxel_sizes = self.write_voxel_sizes
            self.process.inputs.write_interp = self.write_interp

        if self.jobtype == 'est':
            self.process.inputs.bias_regularization = self.bias_regularization
            self.process.inputs.bias_fwhm = self.bias_fwhm
            self.process.inputs.tpm = self.tpm
            (self.process.inputs.
                   affine_regularization_type) = self.affine_regularization_type
            (self.process.inputs.
                           warping_regularization) = self.warping_regularization
            self.process.inputs.smoothness = self.smoothness
            self.process.inputs.sampling_distance = self.sampling_distance

        if self.jobtype == 'estwrite':
            self.process.inputs.bias_regularization = self.bias_regularization
            self.process.inputs.bias_fwhm = self.bias_fwhm
            self.process.inputs.tpm = self.tpm
            (self.process.inputs.
                   affine_regularization_type) = self.affine_regularization_type
            (self.process.inputs.
                           warping_regularization) = self.warping_regularization
            self.process.inputs.smoothness = self.smoothness
            self.process.inputs.sampling_distance = self.sampling_distance
            self.process.inputs.write_bounding_box = self.write_bounding_box
            self.process.inputs.write_voxel_sizes = self.write_voxel_sizes
            self.process.inputs.write_interp = self.write_interp

        self.process.run()


class Realign(Process_Mia):
    """
    *Realigns a time-series of images acquired from the same subject*

    Please, see the complete documention for the `Realign brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/html/documentation/preprocess/spm/Realign.html>`_

    """
    
    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Realign, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['matlab', 'spm']
        
        # Inputs description
        in_files_desc = ('A list  of items with string elements corresponding'
                         ' to existing path files.')
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
                                ' jobtype is write or estwrite, these will be'
                                ' the resliced files (a list of items which are'
                                ' a list of items which are an existing file'
                                ' name).')
        modified_in_files_desc = ('If the jobtype parameter is estimate or'
                                  ' estwrite, these will be copies of the'
                                  ' in_files with a rewritten header (a list of'
                                  ' items which are a list of items which are'
                                  ' an existing file name).')
        mean_image_desc = ('If write_which is [2, 1] or [0, 1] and jobtype is'
                           ' write or estwrite, this will be the mean image'
                           ' file from the realignment (an existing file'
                           ' name).')
        realignment_parameters_desc = ('If the jobtype parameter is estimate'
                                       ' or estwrite, this will be the'
                                       ' estimated translation and rotation'
                                       ' parameters (a list of items which are'
                                       ' an existing file name).')
        
        # Inputs traits
        self.add_trait('in_files',
                       InputMultiPath(traits.Either(ImageFileSPM(),
                                                    traits.List(ImageFileSPM()),
                                                    Undefined),
                                      value=[Undefined],
                                      copyfile=True,
                                      output=False,
                                      optional=False,
                                      desc=in_files_desc))
        
        self.add_trait("jobtype",
                       traits.Enum('estwrite',
                                   'estimate',
                                   'write',
                                   output=False,
                                   optional=True,
                                   desc=jobtype_desc))

        self.add_trait("quality",
                       traits.Range(value=0.9,
                                    low=0.0,
                                    high=1.0,
                                    output=False,
                                    optional=True,
                                    desc=quality_desc))

        self.add_trait("separation",
                       traits.Range(value=4.0,
                                    low=0.0,
                                    high=None,
                                    output=False,
                                    optional=True,
                                    desc=separation_desc))

        self.add_trait("fwhm",
                       traits.Range(value=5.0,
                                    low=0.0,
                                    high=None,
                                    output=False,
                                    optional=True,
                                    desc=fwhm_desc))

        self.add_trait("register_to_mean",
                       traits.Bool(True,
                                   output=False,
                                   optional=True,
                                   desc=register_to_mean_desc))

        self.add_trait("interp",
                       traits.Range(value=2,
                                    low=0,
                                    high=7,
                                    output=False,
                                    optional=True,
                                    desc=interp_desc))

        self.add_trait("wrap",
                       traits.List(value=[0, 0, 0],
                                   trait=traits.Range(low=0, high=1),
                                   minlen=3,
                                   maxlen=3,
                                   output=False,
                                   optional=True,
                                   desc=wrap_desc))

        self.add_trait("weight_img",
                       traits.File(value=Undefined,
                                   output=False,
                                   optional=True,
                                   desc=weight_img_desc))

        self.add_trait('write_which',
                       traits.Enum([2, 1],
                                   [1, 0],
                                   [2, 0],
                                   [0, 1],
                                   output=False,
                                   optional=True,
                                   desc=write_which_desc))
 
        self.add_trait("write_interp",
                       traits.Range(value=4,
                                    low=0,
                                    high=7,
                                    output=False,
                                    optional=True,
                                    desc=write_interp_desc))

        self.add_trait("write_wrap",
                       traits.List(value=[0, 0, 0],
                                   trait=traits.Range(low=0, high=1),
                                   minlen=3,
                                   maxlen=3,
                                   output=False,
                                   optional=True,
                                   desc=write_wrap_desc))

        self.add_trait("write_mask",
                       traits.Bool(default_value=True,
                                   output=False,
                                   optional=True,
                                   desc=write_mask_desc))

        self.add_trait("out_prefix",
                       traits.String(value='r',
                                     output=False,
                                     optional=True,
                                     desc=out_prefix_desc))

        # Output traits
        self.add_trait("realigned_files",
                       OutputMultiPath(traits.Either(traits.List(File()),
                                                     File()),
                                       output=True,
                                       optional=True,
                                       desc=realigned_files_desc))

        self.add_trait("modified_in_files",
                       OutputMultiPath(traits.Either(traits.List(File()),
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
                                       desc=realignment_parameters_desc))  # rp_

        self.process = spm.Realign()
        self.change_dir = True

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return 
        self.make_initResult() object.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Realign, self).list_outputs()
        
        # Outputs definition and tags inheritance (optional)
        if self.outputs:
            self.outputs = {}

        if self.inheritance_dict:
            self.inheritance_dict = {}
        
        if self.in_files and self.in_files != [Undefined]:
            self.process.inputs.in_files = self.in_files

            if self.out_prefix:
                self.process.inputs.out_prefix = self.out_prefix

            self.outputs = self.process._list_outputs()

        if self.outputs:

            for key, values in self.outputs.items():

                if key == 'realigned_files':

                    if ((self.jobtype == 'estimate') or
                     (self.write_which == [0, 1] and
                     (self.jobtype == 'write' or  self.jobtype == 'estwrite'))):
                        self.outputs['realigned_files'] = Undefined

                    else:

                        for fullname in values:
                            path, filename = os.path.split(fullname)

                            if self.out_prefix:
                                filename_without_prefix = filename[len(
                                                              self.out_prefix):]

                            else:
                                filename_without_prefix = filename[len('r'):]

                            if (os.path.join(path,
                                             filename_without_prefix)
                                                 in self.in_files):
                                self.inheritance_dict[fullname] = (
                                    os.path.join(path, filename_without_prefix))

                if key == 'modified_in_files':

                    if self.jobtype == 'write':
                        self.outputs['modified_in_files'] = Undefined

                if key == 'mean_image':

                    if (
                      (self.jobtype == 'estimate') or
                    ((self.write_which == [2, 0] or self.write_which == [1, 0])
                                             and
                     (self.jobtype == 'write' or  self.jobtype == 'estwrite'))
                       ):
                        self.outputs['mean_image'] = Undefined

                    else:

                        if isinstance(values, str):
                            values = [values]

                        for fullname in values:
                            path, filename = os.path.split(fullname)
                            filename_without_prefix = filename[len('mean'):]

                            if (os.path.join(path,
                                             filename_without_prefix)
                                           in self.in_files):
                                self.inheritance_dict[fullname] = (
                                    os.path.join(path,
                                                 filename_without_prefix)
                                                                  )
                if key == 'realignment_parameters':

                    if self.jobtype == 'write':
                        self.outputs['realignment_parameters'] = Undefined

                    else:

                        if isinstance(values, str):
                            values = [values]

                        for fullname in values:
                            path, filename = os.path.split(fullname)
                            filename_without_prefix = filename[len('rp_'):]

                            if filename_without_prefix[-4:] == '.txt':
                                filename_without_prefix = (
                                          filename_without_prefix[:-4] + '.nii')

                            if (os.path.join(path,
                                             filename_without_prefix)
                                           in self.in_files):
                                self.inheritance_dict[fullname] = (
                                    os.path.join(path,
                                                 filename_without_prefix)
                                                                  )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Realign, self).run_process_mia()
        self.process.inputs.in_files = self.in_files
        self.process.inputs.jobtype = self.jobtype
        self.process.inputs.quality = self.quality
        self.process.inputs.separation = self.separation
        self.process.inputs.fwhm = self.fwhm
        self.process.inputs.register_to_mean = self.register_to_mean
        self.process.inputs.interp = self.interp
        self.process.inputs.wrap = self.wrap
        self.process.inputs.weight_img = self.weight_img
        self.process.inputs.write_which = self.write_which
        self.process.inputs.write_interp = self.write_interp
        self.process.inputs.write_wrap = self.write_wrap
        self.process.inputs.write_mask = self.write_mask
        self.process.inputs.out_prefix = self.out_prefix
        self.process.run()

        
class Smooth(Process_Mia):
    """
    *3D Gaussian smoothing of image volumes*

    Please, see the complete documention for the `Smooth brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/html/documentation/preprocess/spm/Smooth.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
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
                       InputMultiPath(traits.Either(ImageFileSPM(),
                                                    Undefined),
                                      value=[Undefined],
                                      copyfile=False,
                                      output=False,
                                      optional=False,
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
        if self.outputs:
            self.outputs = {}

        if self.inheritance_dict:
            self.inheritance_dict = {}

        if self.in_files and self.in_files != [Undefined]:
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
