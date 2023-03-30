:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
TSNR brick
============

Computes the time-course SNR for a time series

Adapted from  `Nipype Cofounds <https://github.com/nipy/nipype/blob/f662acfce8def4717e0c3414618f3a5de5913b31/nipype/algorithms/confounds.py#L974>`_

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import TSNR

>>> TSNR.help()

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz])

    ::

      ex. '/home/username/data/raw_data/func.nii'

**Optional inputs with default value parameters:**

- *prefix_stddev* (a string, optional)
    Prefix of the stddev output. Default is ''.
    
    ::

      ex. ''

- *prefix_tsnr* (a string, optional)
    Prefix of the TSNR output image. Default is ''.
    
    ::

      ex. ''


- *suffix_stddev* (a string, optional)
  Suffix of the STD dev output. Default is '_stddev'.
    
    ::

      ex. '_stddev'
      
- *suffix_tnsr* (a string, optional)
  Suffix of the TSNR output image. Default is '_tsnr'.
    
    ::

      ex. '_tsnr'


**Outputs parameters:**

- *out_tsnr_file* (a strings representing a file)
    Path of the tsnr image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/func_tsnr.nii'

- *out_tsnr_file* (a strings representing a file)
    Path of the  STD dev image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/func_stddev.nii'

-------------

Usefull links:
`Cofunds TSNR nipype <https://github.com/nipy/nipype/blob/f662acfce8def4717e0c3414618f3a5de5913b31/nipype/algorithms/confounds.py#L899>`_
