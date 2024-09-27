:orphan:

.. toctree::

+--------------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+-------------------------------------------+----------------------------------------------------+

==================
BoldIQMsPlot brick
==================

Plot a figure showing the slice-wise signal intensity at the extremes for the identification of spikes,
the outliers metric, the DVARS, the FD and the carpetplot.
The carpet plot rows correspond to voxelwise time series, and are separated into regions:
cortical gray matter, deep gray matter, white matter and cerebrospinal fluid, cerebellum and the brain-edge or “crown".
The crown corresponds to the voxels located on a closed band around the brain.

Used `Niworkflow FMRIplot <https://github.com/nipreps/niworkflows/blob/83bb10ee21d89546b1c5b05f2c7e4ae29ad1bece/niworkflows/viz/plots.py#L38>`_

--------------------------------------

**Mandatory inputs parameters:**

- *in_func* (a string representing an existing file)
    Bold image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func_valid.nii'

**Optional inputs with default value parameters:**

- *drop_trs* (an int, optional, default value is 0)
    Number of dummy scans drops.

    ::

      ex. 2

- *fd_thresh* (a float, optional, default value is 0.2)
    Motion threshold for FD computation.

    ::

      ex. 0.2

**Optional inputs:**

- *carpet_seg*
    Carpet segmentation.

    ::

      ex. '/home/username/data/derived_data/cseg_t_tpl-MNI152NLin2009cAsym_res-01_desc-carpet_dseg.nii.gz'


- *in_dvars_file* (a string representing an existing file, optional)
    DVARS file.

    ::

      ex. '/home/username/data/derived_data/dvars_reg_func_valid.out'


- *in_fd_file* (a string representing an existing file, optional)
    FD file.

    ::

      ex. '/home/username/data/derived_data/fd_reg_func_valid_oned.out'


- *in_outliers_file* (a string representing an existing file, optional)
    Outliers file.

    ::

      ex. '/home/username/data/derived_data/outliers_reg_func_valid.out'


- *in_spikes_file* (a string representing an existing file, optional)
    Spikes file.

    ::

      ex. '/home/username/data/derived_data/spikes_reg_func_valid.out'

- *tr* (a float, optional)
    Repetition time. Default value is Undefined, the TR in the database will be used.

    ::

        ex. 3000.0


**Outputs parameters:**

- *out_file* (a strings representing a file)
    A figure with carpet and outliers/dvars/FD/spikes plot.

    ::

      ex. '/home/username/data/raw_data/func_valid_fmriplot.png'

-------------

Useful links:

`Niworkflow FMRIplot <https://github.com/nipreps/niworkflows/blob/83bb10ee21d89546b1c5b05f2c7e4ae29ad1bece/niworkflows/viz/plots.py#L38>`_
