:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

======================
ExtractSignalROI brick
======================

Extract signals from ROI using a segmentation file with label
-------------------------------------------------------------


**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz]).
    Signals will be extracted from this image.

    ::

      ex. '/home/username/data/raw_data/sub-001_bold.nii'

- *in_seg* (a string representing an existing file)
    Input segmented image (valid extensions: [.nii, .nii.gz]).
    For example a binary brain mask or an image of tissues or structures segmentation
    (each region should be associated to an integer).

    ::

      ex. '/home/username/data/raw_data/native_structures_sub-001_bold_space.nii'

- *labels* (a list of integer)
    List of the labels of the ROI to extract.

    ::

      ex. ['47', '48']

**Outputs parameters:**

- *signals* (a pathlike object or string representing a file)
    Extracted signal for each ROI in a csv file.

    ::

      ex. '/home/username/data/raw_data/sub-001_bold_extracted_signals_47_48.csv'


-------------

Usefull links:

`nilearn NiftiLabelsMasker <https://nilearn.github.io/stable/modules/generated/nilearn.maskers.NiftiLabelsMasker.html>`_
