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
NewSegment brick
================

Segments,  bias  corrects  and  spatially normalises - all in the same model
---------------------------------------------------------------------------- 

>>> from mia_processes.preprocess.spm import NewSegment
>>> NewSegment.help()

**Inputs parameters:**

- *channel_files <=> channel.vols* [#label]_
    Path of the scans for processing (valid extensions, .img, .nii, .hdr).

    ::

      ex. ['/home/ArthurBlair/data/raw_data/Anat.nii']

- *channel_info <=> (channel.biasreg, channel.biasfwhm, (channel.write))* [#label]_
    A tuple (consisting of a float, a float and a tuple consisting of a boolean, a boolean) with the following fields:

    - bias reguralisation (a float between 0 and 10)
        The goal is to model, by different tissue classes, the intensity variations that arise due to different tissues, while model, with a
	bias field, those that occur because of the bias artifact due to the physics of MRI imaging. If the data have very little intensity
	non-uniformity artifact, then bias control should be increased. This effectively tells the algorithm that there is very little bias in
	the data, so it doesn't try to model it.

          | \- 0 No regularisation
          | \- 0.00001 extremely light regularisation
          | \-  ...
          | \- 1 very heavy regularisation
          | \- 10 extremely heavy regularisation

    - bias FWHM (a float between 20 and infinity)
        Full Width at Half Maximum of Gaussian smoothness of bias. Smoother bias fields need fewer parameters to describe them.
	This means that the algorithm is faster for smoother intensity non-uniformities (e.g. 150 mm cutoff gives faster results than 20
	mm cutoff).

    - which maps to save (a tuple of two boolean values; (Field, Corrected))
        To save the estimated bias field or/and the bias corrected version of the processed image.

          | \- (False, False) save Nothing
          | \- (False, True) save bias corrected image only
          | \- (True, False) save estimated bias field only
          | \- (True, True) save estimated bias field and bias corrected image

    ::

      ex. (0.0001, 60, (False, True))

- *tissues <=> [((tissue(i).tpm), tissue(i).ngaus, (tissue(i).native), (tissue(i).warped)), ((tissue(i+1).tpm), tissue(i+1).ngaus,
  (tissue(i+1).native), (tissue(i+1).warped)), ...]* [#label]_
    A list of tuples (one per tissue, i from 1 to 6) with parameter values for each tissue types. Typically, the order of tissues is grey
    matter (i=1), white matter (i=2), CSF (i=3), bone (i=4), soft tissue (i=5) and air/background (i=6), if using tpm/TPM.nii from
    spm12.

    Each tuple consists of the following fields:

        (tissue probability map (4D), 1-based index to frame),  number of gaussians, (which maps to save; Native, DARTEL),
	(which maps to save; Unmodulated, Modulated)

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

- *warping_regularization <=> warp.reg* [#label]_
    The measure of the roughness of the deformations for registration. Involve the sum of 5 elements (a float or list of floats; the latter is
    required by SPM12).

    ::

      ex. [0, 0.001, 0.5, 0.05, 0.2]

- *affine_regularization <=> warp.affreg* [#label]_
    Standard space for affine registration ('mni' or 'eastern' or 'subj' or 'none').

    ::

       ex. mni
 
- *sampling_distance <=> warp.samp* [#label]_
    Approximate distance between sampled points when estimating the model parameters (a float).

    ::

      ex. 3

- *write_deformation_fields <=> warp.write* [#label]_
    Deformation fields can be saved to disk, and used by the deformation utility (a list of 2 booleans for which deformation
    fields to write; Inverse, Forward).

        | \- [False, False] Save nothing
        | \- [True, False] save Inverse only
        | \- [False, True] save Forward only
        | \- etc.

    ::

      ex. [False, True]

**Outputs parameters:**

- *bias_corrected_images*
    The bias corrected images (a list of items which are a pathlike object or string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/raw_data/mAnat.nii

- *bias_field_images*
    The estimated bias field (a list of items which are a pathlike object or string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/raw_data/BiasField_Anat.nii

- *native_class_images*
    Native space probability maps (a list of items which are a list of items which are a pathlike object or string representing an existing
    file).

    ::

      ex. [['/home/ArthurBlair/data/raw_data/c1Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/c2Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/c3Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/c4Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/c5Anat.nii']]

- *dartel_input_images*
    "Imported" class images into a form that can be used with the Dartel toolbox (a list of items which are a list of items which are a pathlike object or string representing an existing file).

    ::

      ex. [['/home/ArthurBlair/data/raw_data/rc1Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/rc2Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/rc3Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/rc4Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/rc5Anat.nii']]

- *modulated_class_images*
    Modulated and normalised class images (a list of items which are a list of items which are a pathlike object or string representing an existing file).

    ::

      ex. [['/home/ArthurBlair/data/raw_data/mwc1Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/mwc2Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/mwc3Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/mwc4Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/mwc5Anat.nii']]

- *normalized_class_images*
    Normalised class images, without modulation (a list of items which are a list of items which are a pathlike object or string representing an existing file).

    ::

      ex. [['/home/ArthurBlair/data/raw_data/wc1Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/wc2Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/wc3Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/wc4Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/wc5Anat.nii']]

- *inverse_deformation_field*
    Inverse deformation field. Could be used for spatially normalising surface files as GIFTI (a list of items which are a pathlike object or
    string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/raw_data/iy_Anat.nii

- *forward_deformation_field*
    Forward deformation field. Could be used for spatially normalising images to MNI space (a list of items which are a pathlike object or
    string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/raw_data/y_Anat.nii

- *transformation_mat*
    Normalisation transformation (a list of items which are a pathlike object or string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/raw_data/Anat_seg8.mat

-------------

.. [#label] Syntax: mia_processes/nipype NewSegment <=> SPM12 Segment.
	    
	    NOTE:
	        - This interface currently supports single channel input only.
	        - The warp.mrf, warp.cleanup and warp.fwhm, from  SPM12, are not used in this brick.
		  
	    Usefull links:
	    `SPM12 Segment <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=45>`_,
	    `nipype NewSegment <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.preprocess.html#newsegment>`_
..
  `nipype <https://nipype.readthedocs.io/en/latest/interfaces/generated/interfaces.spm/preprocess.html#newsegment>`_
