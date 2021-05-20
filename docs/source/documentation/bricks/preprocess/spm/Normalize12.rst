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

=================
Normalize12 brick
=================

Computes the warp that best aligns the template (atlas) to the individual's image
---------------------------------------------------------------------------------

>>> from mia_processes.bricks.preprocess.spm import Normalize12
>>> Normalize12.help()

**Inputs parameters:**

- *image_to_align <=> subj.vol* [#label]_
    The image that the template (atlas) data is warped into alignment with. Mutually exclusive with *deformation_file* parameter. A
    pathlike object or string representing an existing file (valid extensions in [.img, .nii, .hdr]).

    ::

      ex. /home/ArthurBlair/data/raw_data/Anat.nii

- *deformation_file <=> subj.def*  [#label]_
    File y_*.nii containing 3 deformation fields for the deformation in x, y and z dimension. Mutually exclusive with *image_to_align* and
    *tpm* parameters. A pathlike object or string representing an existing file (valid extensions in [.img, .nii, .hdr]).

    ::

      ex. /home/ArthurBlair/data/downloaded_data/y_Anat.nii


- *apply_to_files <=> subj.resample* [#label]_
    Files to apply transformation to. They can be any images that are in register with the image used to generate the deformation. A list
    of items which are an existing, uncompressed file (valid extensions in [.img, .nii, .hdr]).

    ::

      ex. ['/home/ArthurBlair/data/raw_data/Anat.nii']

- *jobtype*
    One of 'write' (write) or  'est' (estimate) or 'estwrite' (estimate and write).

      | \- write: Needs the deformation_file and apply_to_files input parameters. The normalized_files output parameter comes
      | \ |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| from the apply_to_files parameter.

      | \- est: Needs the tpm and image_to_align input parameters. The deformation_field output parameter comes from the
      | \ |ws1| |ws1| |ws1| |ws1| |ws4| image_to_align parameter.

      | \- estwrite: Needs, at least, the tpm and image_to_align input parameters. The deformation_field and normalized_image
      | \ |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| output parameters come from the image_to_align parameter. If the apply_to_files input parameter is defined,
      | \ |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| |ws1| the normalized_files output parameter is also created.
      
    ::

      ex. write

- *bias_regularization <=> eoptions.biasreg*  [#label]_
    The goal is to model, by different tissue classes, the intensity variations that arise due to different tissues, while model, with a bias field, those that occur because of the bias artifact due to the physics of MRI imaging. If the data have very little intensity non-uniformity artifact, then bias control should be increased. This effectively tells the algorithm that there is very little bias in the data, so it doesn't try to model it (a float between 0 and 10).

      | \- 0 No regularisation
      | \- 0.00001 extremely light regularisation
      | \-  ...
      | \- 1 very heavy regularisation
      | \- 10 extremely heavy regularisation

    ::

      ex. 0.0001

- *bias_fwhm <=> eoptions.biasfwhm* [#label]_
    Full Width at Half Maximum of Gaussian smoothness of bias (a value in [30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, ‘Inf’). Smoother bias fields need fewer parameters to describe them. This means that the algorithm is faster for smoother intensity non-uniformities (e.g. 150 mm cutoff gives faster results than 30 mm cutoff).
  
    ::

      ex. 60

- *tpm <=> eoptions.tpm* [#label]_
    The template in form of tissue probability atlas (a pathlike object or string representing an existing file). Mutually exclusive with the *deformation_file* parameter.

    ::

      ex. /home/ArthurBlair/spm12/tpm/TPM.nii

- *affine_regularization_type <=> eoptions.affreg* [#label]_
    Standard space for affine registration (‘mni’ or ‘size’ or ‘none’).

    ::

      ex. mni

- *warping_regularization <=> eoptions.reg* [#label]_
    The measure of the roughness of the deformations for registration. Involve the sum of 5 elements (list of floats).
            
    ::

      ex. [0.0, 0.001, 0.5, 0.05, 0.2]

- *smoothness <=> eoptions.fwhm* [#label]_
    Value to smooth the data before normalization (a float; in mm). 0 is a good value for MRI.

    ::

      ex. 0.0

- *sampling_distance <=> eoptions.samp* [#label]_
    Approximate distance between sampled points when estimating the model parameters (a float).

    ::

      ex. 3.0

- *write_bounding_box <=> woptions.bb* [#label]_
    A list of 2 items which are a list of items which are a float. This is the bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).

    ::

      ex. [[-78, -112, -50], [78, 76, 85]]

- *write_voxel_sizes <=> woptions.vox* [#label]_
    A list of 3 items which are a float. This is the voxel sizes of the written normalised images.

    ::

      ex. [1.0, 1.0, 1.0]

- *write_interp  <=> woptions.interp* [#label]_
    This is the method by which the images are sampled when being written in a different space (0 <= a long integer <= 7).

      | \- 0 Nearest neighbour
      | \- 1 Trilinear (OK for PET, realigned fMRI, or segmentations)
      | \- 2 2nd Degree B-spline
      | \-  ...
      | \- 7 7nd Degree B-spline.

    ::

      ex. 1

**Outputs parameters**


- *deformation_field*
    NIfTI file containing 3 deformation fields for the deformation in x, y and z dimension (a list of items which are a pathlike object or string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/derived_data/y_Anat.nii

- *normalized_files*
    Normalised other files (a list of items which are a pathlike object or string representing an existing file).

    ::
       
      ex. /home/ArthurBlair/data/derived_data/wAnat.nii

-------------

.. [#label] Syntax: mia_processes/nipype Normalize12 <=> SPM12 Normalise.

	    NOTE:
	        - Contrary to what is stated in the nipype documentation, this brick does not accept the out_prefix parameter which would allow to use a prefix other than "w" for the normalised image.

	    Usefull links:
	    `SPM12 Normalise <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=51>`_,
	    `nipype Normalize12 <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.preprocess.html#normalize12>`_
