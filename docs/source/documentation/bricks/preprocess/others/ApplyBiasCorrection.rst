:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
ApplyBiasCorrection brick
============

Apply bias field correction to an input file using the bias image

>>> from mia_processes.bricks.preprocess.others import ApplyBiasCorrection
>>> ApplyBiasCorrection.help()

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz])

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *bias_image* (a string representing an existing file)
    Bias field image (valid extensions: [.nii, .nii.gz])
    
    ::

      ex. '/home/username/data/raw_data/T1w_bias.nii'


**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    Bias corrected file

    ::

      ex. '/home/username/data/derived_data/T1w_inu.nii'

-------------
