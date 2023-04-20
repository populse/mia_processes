:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

================================
Bold_mriqc pipeline
================================

Get no-reference IQMs (image quality metrics) from functional (BOLD)
data using mriqc functional workflow (`mriqc v22.06 <https://mriqc.readthedocs.io/en/22.0.6/>`_).

**Disclaimer**: A special effort has been made to provide a pipeline in Mia that gives as much as possible
the same result as when computing with the native MRIQC. The variation in results between multiple runs
of the same inputs is the result of random sampling and floating point precision errors. These variations
are usually small, but if an identical result is sought between multiple runs, the environment variable
ANTS_RANDOM_SEED should be set (e.g. ANTS_RANDOM_SEED = 1). Further discussion is available
`in a mia_processes ticket <https://github.com/populse/mia_processes/issues/16>`_.

--------------------------------------

**Pipeline insight**

| Bold_mriqc pipeline combines the following pipelines and processes:
|   - `Non steady state detection <../../bricks/preprocess/others/NonSteadyDetector.html>`_
|   - `Sanitize <../../bricks/preprocess/others/Sanitize.html>`_
|   - `TSNR  <../../bricks/preprocess/others/TSNR.html>`_
|   - `Mean <../../bricks/preprocess/afni/TStatMean.html>`_
|   - `Mask <../../bricks/preprocess/afni/Automask.html>`_
|   - `Volume registration <../../bricks/preprocess/afni/Volreg.html>`_
|      (default values : twopass = True, nterpolation = 'Fourier')
|   - `MIN align <../../pipelines/preprocess/Bold_mni_align.html>`_
|   - `IQMS computation <../../pipelines/reports/Bold_iqms.html>`_
|   - `Func report  <../../bricks/reports/ReportFuncMriqc.html>`_

**Mandatory inputs parameters**

- *func_file* (a string representing an existing file)
    A functional image (BOLD). An existing file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'


**Outputs parameters:**

- *func_report*
    Functional Image-Quality Metrics summary report.

    ::

      ex. '/home/username/data/derived_data/func_ref_functional_mriqcReport_2023_03_31_11_42_10_75.pdf'

- *carpet_seg*
    Carpet segmentation.

    ::

      ex. '/home/username/data/derived_data/cseg_t_tpl-MNI152NLin2009cAsym_res-01_desc-carpet_dseg.nii.gz'


-------------

Usefull links:

`mriqc 22.06 <https://mriqc.readthedocs.io/en/22.0.6/>`_
