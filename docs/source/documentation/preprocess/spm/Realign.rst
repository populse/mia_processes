:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+





.. line break
.. |br| raw:: html

   <br />

.. thin space
.. |ws1| raw:: html

   &thinsp;

.. em space

.. |ws2| raw:: html

   &emsp;

.. en space

.. |ws3| raw:: html

   &ensp;

.. non-breakable space

.. |ws4| raw:: html

   &nbsp;

=============
Realign brick
=============

Realigns a time-series of images acquired from the same subject
---------------------------------------------------------------

>>> from mia_processes.preprocess.spm import Realign
>>> Realign.help()

**Inputs parameters:**

- *in_files <=> data* [#label]_
    The images to realign (a list of items which are an existing, uncompressed file (valid extensions: [.img, .nii, .hdr]).

    ::

      ex. ['/home/ArthurBlair/data/raw_data/Func.nii']

- *jobtype*
    One of 'estwrite', 'estimate' or 'write':

      | - estimate: generates realignment_parameters and modified_in_files
      | - write: |ws4| |ws4| with write_which == [2, 0] or [1, 0] generates realigned_files
      | |ws4| |ws4| |ws4| |ws4| |ws4| |ws4| with write_which == [2, 1] generates mean_image and realigned_files
      | |ws4| |ws4| |ws4| |ws4| |ws4| |ws4| with write_which == [0, 1] generates mean_image
      | - estwrite:|ws1| with write_which == [2, 0] or [1, 0] generates realignment_parameters, modified_in_files and realigned_files
      | |ws4| |ws4| |ws4| |ws4| |ws4| |ws4| with write_which == [2, 1] generates realignment_parameters, modified_in_file, mean_image and realigned_files
      | |ws4| |ws4| |ws4| |ws4| |ws4| |ws4| with write_which == [0, 1] generates realignment_parameters, modified_in_file and mean_image

    ::

      ex. estwrite

- *quality <=> eoptions.quality* [#label]_
    Quality versus speed trade-off (0.0 <= a floating point number <= 1.0). Highest quality (1) gives most precise results, whereas lower
    qualities gives faster realignment.

    ::

     ex. 0.9

- *separation <=> eoptions.sep* [#label]_
    Sampling separation in mm in the reference image (a floating point number >= 0.0). Smaller sampling distances gives more accurate
    results, but will be slower.

    ::

      ex. 4.0

- *fwhm <=> eoptions.fwhm* [#label]_
    The gaussian smoothing kernel width (mm, a floating point number >= 0.0) applied to the images before estimating the realignment
    parameters.

    ::

      ex. 5.0

- *register_to_mean <=> eoptions.rtm* [#label]_
    Indicate whether realignment is done to the mean image (True) or to the first image (False).

    ::

      ex. True

- *interp <=> eoptions.interp* [#label]_
    Degree of b-spline (1 <= a long integer <= 7) used for interpolation. Higher degree interpolation methods provide the better
    interpolation, but they are slower because they use more neighbouring voxels.

    ::

      ex. 2

- *wrap <=> eoptions.wrap* [#label]_
    Check if interpolation should wrap in [x,y,z] (a list of 3 items which are integer int or long). For example, in MRI scans, the images wrap
    around in the phase encode direction, so the subject's nose may poke into the back of the subject's head. These are typically:

      | \- No wrapping [0, 0, 0]: for PET or images that have already been spatially transformed. (Also the recommended option if
      | |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| you are not really sure).
      | \- Wrap in Y [0, 1, 0], for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).

    ::

      ex. [0, 0, 0]

- *weight_img <=> eoptions.weight* [#label]_
    Filename of optional weighting image, to weight each voxel of the reference image differently when estimating the realignment
    parameters. This would be used, for example, when there is a lot of extra-brain motion or when there are serious artifacts in a
    particular region of the images.

    ::

      ex. ''

- *write_which <=> roptions.which* [#label]_
    Determines which images to reslice (a list of items which are a value of class 'int'):
    
      | \- [2,0]: Reslices all the images (1..n), including the first image selected, which will remain in its original position.
      | \- [1,0]: Reslices images (2..n) only. Useful for if you wish to reslice (for example) a PET image to fit a structural MRI, without
      | |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| creating a second identical MRI volume.
      | \- [2,1]: All Images + Mean Image. In addition to reslicing the images, it also creates a mean of the resliced image.
      | \- [0,1]: Mean Image Only. Creates the mean resliced image only.

    ::

      ex. [2, 1]

- *write_interp <=> roptions.interp* [#label]_
    The method by which the images are sampled when being written in a different space. Nearest neighbour is fastest, but not
    recommended for image realignment. Trilinear Interpolation is probably OK for PET, but not so suitable for fMRI because higher degree
    interpolation generally gives better results. Although higher degree methods provide better interpolation, but they are slower because
    they use more neighbouring voxels. (0 <= a long integer <= 7). Voxel sizes must all be identical and isotropic.

      | \- 0: Nearest neighbour
      | \- 1: Trilinear
      | \- 2: 2nd Degree B-Spline
      | \- 3: 3rd Degree B-Spline
      | \...
      | \- 7: 7th Degree B-Spline

    ::

      ex. 4

- *write_wrap <=> roptions.wrap* [#label]_
    A list of from 3 items which are an integer (int or long). See the wrap parameter that is used if the jobtype parameter is equal to
    estimate.

      | \- [0, 0, 0]: No wrap
      | \- [1, 0, 0]: wrap X
      | \- [0, 1, 0]: wrap Y
      | \- [0, 0, 1]: Wrap Z
      | \- [1, 1, 1]: Wrap X, Y & Z
      | \...

    ::
     
      ex. [0, 0, 0]

- *write_mask <=> roptions.mask* [#label]_
    Mask output image (a boolean). Because of subject motion, different images are likely to have different patterns of zeros from where it
    was not possible to sample data. With masking enabled, the program searches through the whole time series looking for voxels which
    need to be sampled from outside the original images. Where this occurs, that voxel is set to zero for the whole set of images.

    ::

      ex. True
 
- *out_prefix <=> roptions.prefix* [#label]_
    Realigned output prefix (a string).

    ::

      ex. r

**Outputs parameters:**

- *realigned_files*
    If the write_which parameter is equal to [2, 0], [1, 0] or [2, 1] and jobtype parameter is equal to write or estwrite, these will be the
    resliced files (a list of items which are a list of items which are a pathlike object or string representing an existing file or a pathlike
    object or string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/derived_data/rFunc.nii

- *modified_in_files*
    If the jobtype parameter is equal to estimate or estwrite, these will be copies of the in_files with a rewritten header (a list of items
    which are a list of items which are a pathlike object or string representing an existing file or a pathlike object or string representing an
    existing file).

    ::

      ex. /home/ArthurBlair/data/derived_data/Func.nii

- *mean_image*
    If the write_which parameter is equal to [2, 1] or [0, 1] and jobtype parameter is equal to write or estwrite, this will be the Mean image
    file from the realignment (a pathlike object or string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/derived_data/meanFunc.nii

- *realignment_parameters*
    If the jobtype parameter is equal to estimate or estwrite, this will be the Estimated translation and rotation parameters (a list of items
    which are a pathlike object or string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/derived_data/rp_Func.txt

- *spm_script_file*
    A file necessary for the internal workings of nipype and automatically generated at the run step time (a string representing a file).

    ::

      ex./home/ArthurBlair/scripts/pyscript_realign_459a6c0d-d060-406c-80ce-f40de77692f9.m
 
-------------

.. [#label] Syntax: mia_processes/nipype Realign <=> SPM12 Realign.
	    
	    Usefull links:
	    `SPM12 Realign <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=25>`_, 
	    `nipype Realign <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.preprocess.html#realign>`_
..
	    `nipype <https://nipype.readthedocs.io/en/latest/interfaces/generated/interfaces.spm/preprocess.html#realign>`_
