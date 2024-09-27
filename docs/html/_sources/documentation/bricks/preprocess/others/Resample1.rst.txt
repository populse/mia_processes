:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

================
Resample1 brick
================

Resamples an image to the resolution of a reference image using nilearn.image.resample_to_img()
-----------------------------------------------------------------------------------------------


**Mandatory inputs parameters:**

- *files_to_resample*
    The 3D images that will be resampled (list of pathlike object or string representing a file or a list of items which are a pathlike object or string representing a file ; valid extensions: [.nii, .nii.gz]).

    ::

      ex. ['/home/username/data/raw_data/T1w_002.nii']


- *reference_image*
    A 3D or 4D image used as reference to resample the `files_to_resample` images (a string representing an existing file ; valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/template.nii'


**Optional inputs with default value parameters:**

- *interp*
    The resample method (a string in [‘continuous’, ‘linear’, ‘nearest’]).

    ::

      ex. 'linear'

- *prefix*
    Prefix of the output image (a string, default value is '').

    ::

      ex. ''

- *suffix*
    Suffix of output image (a string, default value is '_003').

    ::

      ex. '_003'

- *suffix_to_delete*
    The suffix to delete from the files_to_resample, when creating the out_files (a string, default value is '_002').
    ::

      ex. '_002'


**Outputs parameters:**

- *out_files*
    The resulting image after resampling ( list of pathlike object or string representing a file or a list of items which are a pathlike object or string representing a file ; extensions: [.nii, .nii.gz]).

    ::

      ex. ['/home/username/data/derived_data/T1w_003.nii']

-------------

Useful links:

 `nilearn.image.resample_to_img() <https://nilearn.github.io/dev/modules/generated/nilearn.image.resample_to_img.html>`_
