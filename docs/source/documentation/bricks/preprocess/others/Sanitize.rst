:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Sanitize brick
============

Sanitize input bold image. 

Adapted from  https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/header.py#L394

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import Sanitize

>>> Sanitize.help()

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image which intensity will be clipped (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

**Optional inputs with default value parameters:**

- *max_32bits* (a boolean, optional)
    Cast data to float32 if higher precision is encountered.

    ::

      default_value. False

- *n_volumes_to_discard* (an int, optional)
    Number of non steady-state volumes.
    
    ::

      default_value. 0

- *prefix* (a string, optional)
    Prefix of the output image.
    
    ::

      default_value. ''

- *suffix* (a string, optional)
   Suffix of output image.
   
    ::

      default_value. '_valid'

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Path of the sanitize image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/func_valid.nii'

-------------

Usefull links:
`SanitizeImage niworflow <https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/header.py#L394>`_
