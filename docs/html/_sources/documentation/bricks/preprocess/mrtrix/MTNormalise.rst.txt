:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=================
MTNormalise brick
=================

Multi-tissue informed log-domain intensity normalisation
--------------------------------------------------------

(mrtrix mtnormalise command)

*Please note that, in mia_processes, MRtrix bricks required diffusion data in MRtrix .mif format.
The MRConvert brick could be used to convert diffusion data in NIfTI format into MRtrix .mif format.
In this case, bvec and bval file should be specified.*

--------------------------------------

**Mandatory inputs parameters:**

- *in_files* (a list of items which are a pathlike object or a string representing an existing file)
    Input tissue component image (valid extensions: [.mif]).

    ::

      ex. ['/home/username/data/derived_data/DWI_wmfod.mif', '/home/username/data/derived_data/DWI_gmfod.mif', '/home/username/data/derived_data/DWI_csffod.mif']


- *mask* (a string representing an existing file)
    The mask defines the data used to compute the intensity normalisation (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/DWI_brainmask.nii'

**Optional inputs with default value parameters:**

- *order_number* (an integer, default value is 3, optional)
    The maximum order of the polynomial basis used to fit the normalisation field in the log-domain.

    ::

      ex. 3

- *niter_number* (an integer or a list of integer, default value is [15, 7], optional)
    Number of iteration. The first (and potentially only) entry applies to the main loop.
    If supplied, the second entry applies to the inner loop to update the balance factors

    ::

      ex. [15, 7]

- *reference_number* (a float, default value is 0.282095, optional)
    The (positive) reference value to which the summed tissue compartments will be normalised

    ::

      ex. 0.282095

- *balanced_number* (a boolean, default value is False, optional)
    Incorporate the per-tissue balancing factors into scaling of the output images

    ::

      ex. False

**Outputs parameters:**

- *out_file* (a list of items which are pathlike object or string representing a file)
    Normalised outputs images
    ::

      ex. ['/home/username/data/derived_data/DWI_wmfod_norm.mif', '/home/username/data/derived_data/DWI_gmfod_norm.mif', '/home/username/data/derived_data/DWI_csffod_norm.mif']

-------------

Useful links:

`mrtrix mtnormalise <https://mrtrix.readthedocs.io/en/latest/reference/commands/mtnormalise.html>`_
