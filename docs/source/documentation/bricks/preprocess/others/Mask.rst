:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Mask brick
============

Apply a binary mask to an image.
Be carrefull, if the suffix and prefix input parameters are not defined or consist only of one or more white spaces, the input parameter will be overwritten.

Adapted from `niworkflows <https://github.com/nipreps/niworkflows>`_.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *mask_file* (a string representing an existing file)
    Mask image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/mask.nii'

**Optional inputs with default value parameters:**

- *prefix* (a string, optional, default value is '')
    Prefix of the output image. Default is ''.

    ::

      ex. ''


- *suffix* (a string, optional, default value is '_masked')
    Suffix of output image. Default is '_masked'.

    ::

      ex. '_masked'



**Outputs parameters:**

- *out_file* (a strings representing a file)
    Path of the masked image (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/T1w_masked.nii'

-------------

Usefull links:
`Mask niworflow <https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/norm.py#L474>`_
