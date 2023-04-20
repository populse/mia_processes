:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
EstimateSNR brick
============

Estimate SNR using a segmentation file

Adapted from `mriqc <https://mriqc.readthedocs.io/en/22.0.6//>`_ anatomical workflow.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *seg_file* (a string representing an existing file)
    Segmentation file used to compute SNR (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

**Outputs parameters:**

- *out_snr* (a float)
    Computed SNR

    ::

      ex. 10.96

-------------


Usefull links:
`mriqc <https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/workflows/anatomical.py#L970>`_
