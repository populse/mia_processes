:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

================
Resample_1 brick
================

Resamples an image to the resolution of a reference image.

/!\ Documentation in progress

--------------------------------------

**Mandatory inputs parameters:**

- *files_to_resample* (list of pathlike object or string reprsenting  a file or a list of items which are a pathlike object or string representing a file)
    The 3D images that will be resampled (valid extensions: [.nii, .nii.gz]).

    ::

      ex. ['/home/username/data/raw_data/T1w_002.nii']


- *reference_image* (a string representing an existing file)
    A 3D or 4D image used as reference to resample the files_to_resample images (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/template.nii'


**Optional inputs with default value parameters:**

- *interp* (an integer between 0 and 5, default value is 3)
    The order of the spline interpolation.
    Default value is trilinear (3).

    ::

      ex. 3

- *prefix* (a string, optional, default value is '')
    Prefix of the output image.

    ::

      ex. ''

- *suffix* (a string, optional, default value is '_003')
    Suffix of output image.

    ::

      ex. '_003'

- *suffix_to_delete* (a string, optional, default value is '_002')
    The suffix to delete from the files_to_resample, when creating the out_files.
    ::

      ex. '_002'


**Outputs parameters:**

- *out_files* (list of pathlike object or string reprsenting  a file or a list of items which are a pathlike object or string representing a file)
    TThe resulting image after resampling ( extensions: [.nii, .nii.gz]).

    ::

      ex. ['/home/username/data/raw_data/T1w_003.nii']

-------------
