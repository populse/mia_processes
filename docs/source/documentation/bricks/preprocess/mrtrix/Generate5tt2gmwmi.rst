:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=======================
Generate5tt2gmwmi brick
=======================

Generate a mask image appropriate for seeding streamlines on the grey matter-white matter interface
---------------------------------------------------------------------------------------------------

(mrtrix 5tt2gmwmi command)

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input 5TT segmented anatomical image (valid extensions: [.mif, .nii, .nii.gz]).
    From Generate5ttfsl brick for example.

    ::

      ex. '/home/username/data/raw_data/T1w_5tt.mif'

**Optional inputs parameters:**

- *in_mask* (a string representing an existing file)
    Filter an input mask image (valid extensions: [.mif, .nii, .nii.gz]) according to those voxels that lie upon the grey matter - white matter boundary.

    If no input mask is provided, the output will be a whole-brain mask image calculated using the anatomical image only.

    ::

      ex. '/home/username/data/derived_data/brainmask.mif'


**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    The output mask image.

    ::

      ex. '/home/username/data/derived_data/T1w_5tt_gmwmSeed.mif'

-------------

Usefull links:

`mrtrix 5tt2gmwmi <https://mrtrix.readthedocs.io/en/latest/reference/commands/5tt2gmwmi.html>`_
