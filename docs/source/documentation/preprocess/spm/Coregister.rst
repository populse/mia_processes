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

================
Coregister brick
================

Align together scans of different modalities
--------------------------------------------

>>> from mia_processes.preprocess.spm import Coregister
>>> Coregister.help()

**Inputs parameters:**

- *target <=> ref* [#label]_
    The reference file (remains stationary) while the source image is moved to match it. An existing, uncompressed file (valid extensions:
    [.img, .nii, .hdr]).

    ::

      ex. /home/ArthurBlair/data/downloaded_data/meanFunc.nii

- *source <=> source* [#label]_
    The image that is jiggled about to best match the target image. A list of items which are an existing, uncompressed file (valid
    extensions: [.img, .nii, .hdr]).

    ::

      ex. ['/home/ArthurBlair/data/raw_data/Anat.nii']

- *apply_to_files <=> other* [#label]_
    These are any images that need to remain in alignment with the source image (a list of items which are an existing file name).
    
    ::

      ex. ['/home/ArthurBlair/data/raw_data/Func.nii']



- *jobtype* [#label]_
    One of 'estwrite' or 'estimate' or 'write'. If 'estimate' is selected, the registration parameters are stored in the headers of the 'source'
    and the 'apply_to_files' images. If 'write' is selected, the resliced images are named the same as the originals except that they are
    prefixed by out_prefix. if 'estwrite' is selected, the described procedures for 'estimate' and 'write' are performed and the output parameter
    (coregistered_source and/or coregistered_files) contains the resliced images and the one whose header has been rewritten. If it is necessary
    to choose one or the other for a subsequent calculation in a pipeline, the Filter_Files_List brick (mia_processes library) can be used.

    ::

       ex. estimate

- *cost_function <=> eoptions.cost_fun* [#label]_
    One of 'mi' or 'nmi' or 'ecc' or 'ncc'. Registration involves finding parameters that either maximise or minimise some objective
    function. For inter-modal registration,  use 'Mutual Information', 'Normalised Mutual Information' or 'Entropy Correlation Coefficient'. For
    within modality, you could also use Normalised Cross Correlation. 

      | \- 'mi': Mutual Information
      | \- 'nmi': Normalised Mutual Information
      | \- 'ecc': Entropy Correlation Coefficient
      | \- 'ncc': Normalised Cross Correlation

    ::

      ex. nmi

- *separation <=> eoptions.sep* [#label]_
    A list of items which are a float. The average distance between sampled points (in mm). Can be a vector to allow a coarse registration
    followed by increasingly fine ones.

    ::

      ex. [4, 2]

- *tolerance <=> eoptions.tol* [#label]_
    A list of 12 items which are a float. The acceptable tolerance for each of 12 params. Iterations stop when differences between
    successive estimates are less than the required tolerance.

    ::

      ex. [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001]

- *fwhm <=> eoptions.fwhm* [#label]_
    A list of 2 items which are a float. Kernel of gaussian smooth to apply to the 256*256 joint histogram.

    ::

      ex. [7, 7] 

- *write_interp <=> roptions.interp* [#label]_
    The method by which the images are sampled when being written in a different space. Nearest neighbour is fastest, but not
    recommended for image realignment. Trilinear Interpolation is probably OK for PET, or realigned and re-sliced fMRI, but not so suitable
    for fMRI with subject movemen because higher degree interpolation generally gives better results. Although higher degree methods
    provide better interpolation, but they are slower because they use more neighbouring voxels. (0 <= a long integer <= 7). Voxel sizes
    must all be identical and isotropic.

      | \- 0: Nearest neighbour
      | \- 1: Trilinear
      | \- 2: 2nd Degree B-Spline
      | \- 3: 3rd Degree B-Spline
      | …
      | \- 7: 7th Degree B-Spline

    ::

      ex. 4

- *write_wrap <=> roptions.wrap* [#label]_
    Check if interpolation should wrap in [x,y,z] (a list of 3 items which are integer int or long). For example, in MRI scans, the images wrap
    around in the phase encode direction, so the subject’s nose may poke into the back of the subject’s head. These are typically:

      | \- No wrapping [0, 0, 0]: for PET or images that have already been spatially transformed (Also the recommended option if
      | |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| you are not really sure)
      | \- Wrap in Y [0, 1, 0], for (un-resliced) MRI where phase encoding is in the Y direction (voxel space)

    ::

        ex. [0 0 0]

- *write_mask <=> roptions.mask* [#label]_
    Mask output image (a boolean). Because of subject motion, different images are likely to have different patterns of zeros from where it
    was not possible to sample data. With masking enabled, the program searches through the whole time series looking for voxels which
    need to be sampled from outside the original images. Where this occurs, that voxel is set to zero for the whole set of images.

    ::

      ex. False

- *out_prefix <=> roptions.prefix* [#label]_
    Specify the string to be prepended to the filenames of the coregisterd image file(s).

    ::

      ex. r, capsul/nipype default value

**Outputs parameters:**

- *coregistered_source*
    A list of items which are an existing file name. Coregistered source files, corresponding to 'source' images.

    ::

      ex. /home/ArthurBlair/data/raw_data/Anat.nii

- *coregistered_files*
    A list of items which are an existing file name. Coregistered other files, corresponding to 'apply_to_files' images.

    ::

      ex. /home/ArthurBlair/data/raw_data/Func.nii

-------------

.. [#label] Syntax: mia_processes/nipype Coregister <=> SPM12 Coregister.
	    
	    Usefull links:
	    `SPM12 Coregister <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=39>`_,
	    `nipype Coregister <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.preprocess.html#coregister>`_
..
  `nipype <https://nipype.readthedocs.io/en/latest/interfaces/generated/interfaces.spm/preprocess.html#coregister>`_
