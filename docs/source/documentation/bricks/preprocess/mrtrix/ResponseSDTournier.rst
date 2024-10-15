:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==========================
ResponseSTrournie brick
==========================

Estimate response function for spherical deconvolution using the Tournier algorithm
----------------------------------------------------------------------------------------

Use the Tournier iterative algorithm for single-fibre voxel selection and response function estimation

(mrtrix dwi2response tournier command)

*Please note that, in mia_processes, MRtrix bricks required diffusion data in MRtrix .mif format.
The MRConvert brick could be used to convert diffusion data in NIfTI format into MRtrix .mif format.
In this case, bvec and bval file should be specified.*

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a pathlike object or a string representing an existing file)
    Input DWI image (valid extensions: [.mif]).

    ::

      ex. '/home/username/data/derived_data/DWI.mif'


**Optional inputs with default value parameters:**

- *get_final_voxels* (a boolean, default value is False, optional)
    Get an image showing the final voxel selection(s).
    The output image shows which voxels from the image were used to construct the basus functions for each tisuue type.

    ::

      ex. False


**Optional inputs:**
- *number* (an integer, optional)
    Number of single-fibre voxels to use when calculating response function

    ::

      ex.6

- *iter_voxels* (an integer, optional)
    Number of single-fibre voxels to select when preparing for the next iteration

    ::

      ex.6

- *dilate* (an integer, optional)
    Number of mask dilation steps to apply when deriving voxel mask to test in the next iteration

    ::

      ex.6

- *max_iters* (an integer, optional)
    Maximum number of iterations (set to 0 to force convergence)

    ::

      ex.6

- *in_mask* (a string representing an existing file, optional)
    Provide initial mask image (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/DWI_brainmask.nii'

- *max_sh* ( an integer, optional)
    Maximum harmonic degree of response function (a single value for single-shell response)

    ::

      ex.6


**Outputs parameters:**

- *wm_file* (a pathlike object or string representing a file)
    Output WM response text file

    ::

      ex. '/home/username/data/derived_data/DWI_response_wm.txt'

- *voxels_image* (a pathlike object or string representing a file)
    Image showing the final voxel selection

    ::

      ex. '/home/username/data/derived_data/DWI_response_voxels.mif'

-------------

Useful links:

`mrtrix dwi2response tournier <https://mrtrix.readthedocs.io/en/dev/reference/commands/dwi2response.html#dwi2response-tournier>`_

`mrtrix dwi2response - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.preprocess.html#responsesd>`_
