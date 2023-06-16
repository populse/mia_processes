:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==============
Binarize brick
==============

Binarizes the input image applying the given lower threshold
------------------------------------------------------------

All voxels lower than the threshold are set to 0 and all the others voxels are set to 1.

Be carrefull, if the suffix and prefix input parameters are not defined or consist only of one or more white spaces, the input parameter will be overwritten.

Adapted from `niworkflows <https://github.com/nipreps/niworkflows>`_.

--------------------------------------


**Mandatory inputs parameters:**

- *in_files* (A list of items with string elements corresponding to existing path files.)
    Input images (valid extensions: [.nii, .nii.gz]).

    ::

      ex. ['/home/username/data/raw_data/T1w.nii']

**Optional inputs with default value parameters:**

- *prefix* (a string, optional, default value is '')
    Prefix of the output images.

    ::

      ex. ''


- *suffix* (a string, optional, default value is '_bin')
    Suffix of output images.

    ::

      ex. '_bin'

- *thresh_low* (a float, default value is 0.0)
    Lower threshold for binarization.

    ::

      ex. 0.0


**Outputs parameters:**

- *out_files* (a list of pathlike objects or strings representing a file)
    Images after application of the binarization (extensions: [.nii, .nii.gz]).

    ::

      ex. ['/home/username/data/derived_data/T1w_bin.nii']

-------------

Usefull links:

`Binarize niworkflows <https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/nibabel.py#L92>`_
