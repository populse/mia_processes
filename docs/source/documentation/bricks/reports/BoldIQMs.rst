:orphan:

.. toctree::

+--------------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+-------------------------------------------+----------------------------------------------------+

==============
BoldQMs brick
==============

Computes the functional Image Quality Metrics (IQMs) as defined in mriqc v22.06
-------------------------------------------------------------------------------

A no-reference IQM is a measurement of some aspect of the actual image which cannot be compared to a reference value for
the metric since there is no ground-truth about what this number should be.
All the computed IQMs corresponding to an image are saved in a JSON file.

Adapted from `mriqc <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/functional.py#L243>`_.

--------------------------------------

**Mandatory inputs parameters:**

- *in_epi* (a string representing an existing file)
    Mean bold image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/mean_reg_func_valid.nii'

**Optional inputs with default value parameters:**

- *in_fd_thresh* (a float, optional, default value is 0.2)
    Motion threshold for FD computation

    ::

      ex. 0.2

**Optional inputs:**

- *in_dummy_TRs* (an integer , optional)
    Number of dummy scans.

    ::

      ex. 2

- *in_dvars_file* (a string representing an existing file, optional)
    DVARS file.

    ::

      ex. '/home/username/data/derived_data/dvars_reg_func_valid.out'


- *in_fd_file* (a string representing an existing file, optional)
    FD file.

    ::

      ex. '/home/username/data/derived_data/fd_reg_func_valid_oned.out'


- *in_fwhm_file* (a string representing an existing file, optional)
    A file with FWHM information.

    ::

      ex. '/home/username/data/derived_data/fwhm_mean_reg_func_valid.out'

- *in_gcor* (a float)
    Global correlation value.

    ::

      ex. 0.0612

- *in_hmc"* (a string representing an existing file, optional)
    Motion corrected input image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/reg_func_valid.nii'


- *in_mask* (a string representing an existing file, optional)
    Brain mask (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/automask_reg_func_valid.nii'



- *in_outliers_file* (a string representing an existing file, optional)
    Outliers file.

    ::

      ex. '/home/username/data/derived_data/outliers_reg_func_valid.out'


- *in_QI_file* (a string representing an existing file, optional)
    Quality index file.

    ::

      ex. '/home/username/data/derived_data/QI_reg_func_valid.out'


- *in_spikes_file* (a string representing an existing file, optional)
    Spikes file.

    ::

      ex. '/home/username/data/derived_data/spikes_reg_func_valid.out'

- *in_tsnr* (a string representing an existing file, optional)
    tSNR input volume. (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/reg_func_valid_tsnr.nii'

**Outputs parameters:**

- *out_file* (a strings representing a file)
    A json file with all the IQMs computed.

    ::

      ex. '/home/username/data/derived_data/mean_reg_func_valid_functional_qc.json'

-------------

Usefull links:

`mriqc func qc  <https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/qc/functional.py>`_
