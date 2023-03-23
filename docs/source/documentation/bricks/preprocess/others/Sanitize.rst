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

**Inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image which intensity will be clipped (valid extensions: [.nii, .nii.gz])

    ::

      ex. '/home/username/data/raw_data/func.nii'

- *prefix* (a string, optional=True)
    Prefix of the output image. Default is ''.
    
    ::

      ex. ''

- *suffix* (a string, optional=True)
   Suffix of output image. Default is '_valid'.
    
    ::

      ex. '_valid'

- *n_volumes_to_discard* (an int, optional=True)
    Number of non steady-state volumes. Default is 0
    
    ::

      ex. 0

- *max_32bits* (a boolean, optional=True)
    Cast data to float32 if higher precision is encountered. Default is False

    ::

      ex. False



**Outputs parameters:**

- *out_file* (a strings representing a file)
    Path of the sanitize image 
    
    ::

      ex. '/home/username/data/raw_data/func_valid.nii'

-------------

Usefull links:
`SanitizeImage niworflow <https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/header.py#L394>`_
