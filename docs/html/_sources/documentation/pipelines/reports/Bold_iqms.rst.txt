:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

==================
Bold_iqms pipeline
==================

Compute no-reference IQMs (image quality metrics) from functional (BOLD).
Used in bold_mriqc pipeline.

--------------------------------------

**Pipeline insight**

| Bold_iqms pipeline combines the following pipelines and processes:
|   - `Outlier Count <../../bricks/reports/OutlierCount.html>`_
|   - `Carpet parcellation  <../../bricks/reports/CarpetParcellation.html>`_
|   - `Compute frame wise displacement  <../../bricks/reports/FramewiseDisplacement.html>`_
|   - `Spikes  <../../bricks/reports/Spikes.html>`_
|   - `Compute GCor  <../../bricks/reports/GCor.html>`_
|   - `Compute DVARS  <../../bricks/reports/ComputeDVARS.html>`_
|   - `Compute FWHMX  <../../bricks/reports/FWHMX.html>`_
|   - `Compute quality index  <../../bricks/reports/QualityIndex.html>`_
|   - `Compute IQMs <../../bricks/reports/BoldIQMS.html>`_

**Mandatory inputs parameters**

- *brainmask* (a string representing an existing file)
    Brain mask. An existing file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/automask_mean_reg_func_valid.nii'

- *epi_mean* (a string representing an existing file)
    Mean functional image (BOLD). An existing file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/mean_reg_func_valid.nii'

- *epi_parc* (a string representing an existing file)
    Template data in bold space (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/t_tpl-MNI152NLin2009cAsym_res-01_desc-carpet_dseg.nii'

- *hmc_motion* (a string representing an existing file)
    Movement parameters file (extensions: [.txt]).

    ::

      ex. '/home/username/data/derived_data/reg_d_func_oned.txt'

- *ras_epi* (a string representing an existing file)
    Functional image (BOLD). An existing file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func_valid.nii'

**Optional inputs parameters**

- *dummy_TRS* (an interger, optional)
    Volume to discard

    ::

      ex. 3

- *dvars_intensity_normalization* (a float, optional)
    See `dvars brick <../../bricks/reports/ComputeDvars.html>`_

    ::

      ex. 1000.0

- *dvars_remove_zero_variance* (a boolean, optional)
    See `dvars brick <../../bricks/reports/ComputeDvars.html>`_

    ::

      ex. True

- *epi_tsnr* (a string representing an existing file, optional)
    See `TSNR brick <../../bricks/others/TSNR.html>`_

    ::

      ex. '/home/username/data/derived_data/reg_func_valid_tsnr.nii'

- *fd_normalize* (a boolean, optional)
    See `Frame displacement brick <../../bricks/others/FramewiseDisplacement.html>`_

    ::

      ex. False

- *fd_parameter_source* (optional)
    See `Frame displacement brick <../../bricks/others/FramewiseDisplacement.html>`_

    ::

      ex. 'AFNI'

- *fd_radius* (a float, optional)
    See `Frame displacement brick <../../bricks/others/FramewiseDisplacement.html>`_

    ::

      ex. 50.0

- *fd_thresh* (a float, optional)
    Frame displacement threshold

    ::

      ex. 0.2

- *fwhm_combine* (a boolean, optional)
    See `FWHMx brick <../../bricks/others/FWHMx.html>`_

    ::

      ex. True

- *fwhm_detrend* (a boolean, optional)
    See `FWHMx brick <../../bricks/others/FWHMx.html>`_

    ::

      ex. True

- *hmc_epi* (a string representing an existing file, optional)
    Head motion corrected functional image

    ::

      ex. '/home/username/data/raw_data/reg_func_valid.nii'


- *quality_index_automask* (a boolean, optional)
    See `quality index brick <../../bricks/reports/QualityIndex.html>`_

    ::

      ex. True

- *outlier_fraction* (a boolean, optional)
    See `Outlier Count brick <../../bricks/preprocess/others/OutlierCount.html>`_

    ::

      ex. True

- *spikes_detrend* (a boolean, optional)
    See `Spikes brick <../../bricks/reports/Spikes.html>`_

    ::

      ex. False

- *spikes_no_zscore* (a boolean, optional)
    See `Spikes brick <../../bricks/reports/Spikes.html>`_

    ::

      ex. True

- *spikes_skip_frames* (an integer, optional)
    See `Spikes brick <../../bricks/reports/Spikes.html>`_

    ::

      ex. 0

- *spikes_skipe_thresh* (a float, optional)
    See `Spikes brick <../../bricks/reports/Spikes.html>`_

    ::

      ex. 6.0

**Outputs parameters:**

- *BoldQC_out_file*
    Functional Image-Quality Metrics summary json file.

    ::

      ex. '/home/username/data/derived_data/mean_reg_func_valid_bold_qc.json'

- *carpet_seg*
    Carpet segmentation.

    ::

      ex. '/home/username/data/derived_data/cseg_t_tpl-MNI152NLin2009cAsym_res-01_desc-carpet_dseg.nii.gz'


-------------

Usefull links:

`mriqc 22.06 <https://mriqc.readthedocs.io/en/22.0.6/>`_
