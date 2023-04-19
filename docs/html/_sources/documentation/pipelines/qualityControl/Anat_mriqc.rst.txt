:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

================================
Anat_mriqc pipeline
================================

Get no-reference IQMs (image quality metrics) from structural (T1w and T2w)
data using mriqc anatomical workflow (`mriqc v22.06 <https://mriqc.readthedocs.io/en/22.0.6/>`_).

**Disclaimer**: A special effort has been made to provide a pipeline in Mia that gives as much as possible
the same result as when computing with the native MRIQC. The variation in results between multiple runs
of the same inputs is the result of random sampling and floating point precision errors. These variations
are usually small, but if an identical result is sought between multiple runs, the environment variable
ANTS_RANDOM_SEED should be set (e.g. ANTS_RANDOM_SEED = 1). Further discussion is available
`in a mia_processes ticket <https://github.com/populse/mia_processes/issues/16>`_.

--------------------------------------

**Pipeline insight**

| Anat_mriqc pipeline combines the following pipelines and processes:
|   - `ConformImage <../../bricks/preprocess/others/ConformImage.html>`_
|   - `Harmonize <../../bricks/preprocess/others/Harmonize.html>`_
|   - `Skull Stripping <../../pipelines/preprocess/Anat_skullstrip_synthstrip.html>`_ (using SynthStrip from Freesurfer)
|   - `Segmentation <../../bricks/preprocess/fsl/FastSegment.html>`_ (using Fast from FSL)
|   - `Spatial Normalisation <../../pipelines/preprocess/Anat_spatial_norm.html>`_
|   - `Head mask <../../bricks/preprocess/fsl/FastSegment.html>`_ (using Bet from FSL)
|   - `Air mask  <../../pipelines/preprocess/Anat_airmask.html>`_
|   - `Anat_mni_tpms <../../pipelines/preprocess/Anat_mni_tpms.html>`_
|   - `IQMs computation  <../../bricks/reports/AnatIQMs.html>`_
|   - `FWHMx computation  <../../bricks/reports/FWHMx.html>`_
|   - `Anat report  <../../bricks/reports/ReportFuncMriqc.html>`_

**Mandatory inputs parameters**

- *anat_file* (a string representing an existing file)
    An anatomical image (T1w or T2w). An existing file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'


**Outputs parameters:**

- *anat_report*
    Anatomical Image-Quality Metrics summary report. 

    ::

      ex. '/home/username/data/derived_data/T1w_ref_anatomical_mriqcReport_2023_03_31_11_42_10_75.pdf'


-------------

Usefull links:

`mriqc 22.06 <https://mriqc.readthedocs.io/en/22.0.6/>`_
