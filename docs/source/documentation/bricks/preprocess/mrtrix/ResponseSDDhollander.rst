:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==========================
ResponseSDDhollander brick
==========================

Estimate response function(s) for spherical deconvolution using the Dhollander algorithm
----------------------------------------------------------------------------------------

Unsupervised estimation of WM, GM and CSF response functions that does not require a T1 image.

The aime is to decompose the diffusion signal into a set of smaller individual fiber orientations (ie deconvolve the fiber orientation distributions).

(mrtrix dwi2response dhollander command)

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a pathlike object or a string representing an existing file)
    Input DWI image (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/DWI.mif'


**Optional inputs with default value parameters:**

- *erode* (an integer, default value is 3, optional)
    Number of erosion passes to apply to initial (whole brain) mask.
    Set to 0 to not erode the brain mask.

    ::

      ex. 3

- *fa_thresh* (a float, default value is 0.2, optional)
    FA threshold for crude WM versus GM-CSF separation

    ::

      ex. [15, 7]

- *get_final_voxels* (a boolean, default value is False, optional)
    Get an image showing the final voxel selection(s).
    The output image shows which voxels from the image were used to construct the basus functions for each tisuue type.

    ::

      ex. False


**Optional inputs:**

- *in_mask* (a string representing an existing file, optional)
    Provide initial mask image (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/DWI_brainmask.nii'

- *max_sh* (a list of items which are an integer, optional)
    Maximum harmonic degree of response function (a single value for single-shell response and a list for multi-shell response)

    ::

      ex. [40]

- *wm_algo* (fa or tax or tournier, default value is False, optional)
    Use external algorithm for WM single-fibre voxel selection.
    By default this option is not used and the algorithm used is the built-in Dhollander 2019

    ::

      ex. fa


**Outputs parameters:**

- *csf_file* (a pathlike object or string representing a file)
    Output CSF response text file

    ::

      ex. '/home/username/data/derived_data/DWI_response_csf.txt'

- *wm_file* (a pathlike object or string representing a file)
    Output WM response text file

    ::

      ex. '/home/username/data/derived_data/DWI_response_wm.txt'

- *gm_file* (a pathlike object or string representing a file)
    Output GM response text file

    ::

      ex. '/home/username/data/derived_data/DWI_response_gm.txt'

- *voxels_image* (a pathlike object or string representing a file)
    Output CSF response text file

    ::

      ex. '/home/username/data/derived_data/DWI_response_voxels.mif'

-------------

Usefull links:

`mrtrix dwi2response dhollander <https://mrtrix.readthedocs.io/en/latest/reference/commands/dwi2response.html#dwi2response-dhollander>`_

`mrtrix dwi2response - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.preprocess.html#responsesd>`_
