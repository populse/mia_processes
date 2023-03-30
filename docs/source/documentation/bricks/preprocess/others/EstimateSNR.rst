:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
EstimateSNR brick
============

Estimate SNR using a segmentation file 

Adapted from https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/workflows/anatomical.py#L970
--------------------------------------

>>> from mia_processes.bricks.preprocess.others import EstimateSNR
>>> EstimateSNR.help()

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *seg_file* (a string representing an existing file)
    Segmentation file used ti calculate SNR (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

**Outputs parameters:**

- *out_snr* (a string representing a file)
    Output Images (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/T1w_estimated_snr.txt

-------------


Usefull links:
`mriqc <https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/workflows/anatomical.py#L970>`_