:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=================
DWIDenoise brick
=================

Denoise DWI data
------------------

DWI data denoising and noise map estimation using Marchenko-Pastur PCA.

If used, this steps should be performed as the first step of the image processing pipeline (before interpolation or smoothing).

(mrtrix dwidenoise command)

*Please note that, in mia_processes, MRtrix bricks required diffusion data in MRtrix .mif format.
The MRConvert brick could be used to convert diffusion data in NIfTI format into MRtrix .mif format.
In this case, bvec and bval file should be specified.*

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input DWI image (valid extensions: [.mif]).

    ::

      ex. '/home/username/data/raw_data/DWI.mif'

**Optional inputs with default value parameters:**

- *extend* (a tuple of the form: (an integer, an integer, an integer), default value is (5, 5, 5), optional)
    Set the window size of the denoising filter

    ::

      ex. (5, 5, 5)

- *noise* (a boolean, default value is True, optional)
    Get noise map.

    ::

      ex. True

**Optional inputs parameters:**

- *in_mask* (a string representing an existing file, optional)
    Input mask image, only process voxels within those voxels (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/DWI_brainmask.mif'

**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    The output denoised DWI image

    ::

      ex. '/home/username/data/derived_data/DWI_denoised.mif'

- *noise_map* (a pathlike object or string representing a file)
    The output noise map

    ::

      ex. '/home/username/data/derived_data/DWI_noise_map.mif'

-------------

Usefull links:

`mrtrix dwidenoise <https://mrtrix.readthedocs.io/en/latest/reference/commands/dwidenoise.html>`_

`mrtrix dwidenoise - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.preprocess.html#dwidenoise>`_
