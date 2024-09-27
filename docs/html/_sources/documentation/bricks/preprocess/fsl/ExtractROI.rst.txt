:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

================
ExtractROI brick
================

Extract region of interest (ROI) from an image using fslroi (FSL)
-----------------------------------------------------------------

| It is possible to:
| - take a 3D ROI from a 3D data set (or if it is 4D, the same ROI is taken from each time point and a new 4D data set is created)
| - extract just some time points from a 4D data set
| - control time and space limits to the ROI.

Note that the arguments are minimum index and size (not maximum index).
So to extract voxels 10 to 12 inclusive you would specify 10 and 3 (not 10 and 12).

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w_5tt.nii'

**Optional inputs with default value parameters:**

- *suffix* (a string, optional, default value is roi)
    Output suffix

    ::

      ex. 'roi'

- *output_type* ('NIFTI' or 'NIFTI_GZ', optional, default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   - NIFTI: \*.nii
    |   - NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

**Optional inputs**

- *t_min* (an integer, optional)
    t minimum index.

    ::

      ex. 0

- *t_size* (an integer, optional)
    t size.

    ::

      ex. 1

- *x_min* (an integer, optional)
    x minimum index.

    ::

      ex. 0

- *x_size* (an integer, optional)
    x size.

    ::

      ex. 1

- *y_min* (an integer, optional)
    y minimum index.

    ::

      ex. 0

- *y_size* (an integer, optional)
    y size.

    ::

      ex. 1

- *z_min* (an integer, optional)
    z minimum index.

    ::

      ex. 0

- *z_size* (an integer, optional)
    z size.

    ::

      ex. 1

**Outputs parameters:**

- *roi_file* (a pathlike object or string representing a file)
    Output file

    ::

      ex. '/home/username/data/derived_data/5tt_roi.nii'

-------------

Useful links:

`FSL fslroi <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Fslutils>`_

`FSL fslroi - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.fsl.preprocess.html#fast>`_
