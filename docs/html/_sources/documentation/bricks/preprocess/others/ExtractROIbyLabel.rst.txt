:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=======================
ExtractROIbyLabel brick
=======================

Extract a specific ROI from a segmentation file using labels
------------------------------------------------------------


**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input segmented image (valid extensions: [.nii, .nii.gz]).
    For example a binary brain mask or an image of tissues or structures segmentation
    (each region should be associated to an integer).

    ::

      ex. '/home/username/data/raw_data/native_structures_sub-001.nii'

**Optional inputs with default value parameters:**

- *labels* (a list of integer, default is [1])
    List of the labels of the ROI to extract.

    ::

      ex. ['47', '48']


- *save_each_roi* (a boolean, optional, default value is False)
    Save a each regions as a NIfTI file.

    ::

      ex. True

- *save_concate_roi* (a boolean, optional, default value is False)
    Save a concatenation of all the ROI extrated as a NIfTI file.

    ::

      ex. True


**Outputs parameters:**

- *out_files* (list of pathlike object or string representing a file)
    Extracted ROI files

    ::

      ex. ['/home/username/data/raw_data/native_structures_sub-001_47.nii',
           '/home/username/data/raw_data/native_structures_sub-001_48.nii'
        ]

- *out_concate* (a pathlike object or string representing a file)
    Image with all ROI

    ::

      ex. '/home/username/data/raw_data/native_structures_sub-001_concate_47_48.nii'
