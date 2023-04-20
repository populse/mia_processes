:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

================================
Anat_headmask pipeline
================================

Compute head mask from skull stripped structural image using "Denoise" from Dipy.

Adapted from `mriqc <https://mriqc.readthedocs.io/en/22.0.6//>`_ anatomical workflow.

--------------------------------------

**Pipeline insight**

| Anat_headmask pipeline combines the following pipelines and processes:
|   - `Estimate SNR <../../bricks/preprocess/others/EstimateSNR.html>`_
|   - `Enhance <../../bricks/preprocess/others/Enhance.html>`_
|   - `Denoise <../../bricks/preprocess/dipy/Denoise.html>`_
|   - `Gradient Threshold <../../bricks/preprocess/others/GradientThreshold.html>`_

.. image:: ../../images/Anat_headmask_pipeline.png
  :width: 800
  :alt: Anat headmask pipeline

**Mandatory inputs parameters**

- *in_file* (a string representing an existing file)
    Skull stripped anatomical image (T1w or T2w) (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/ss_T1w.nii'

- *seg_file* (a string representing an existing file)
    Segmented image of the in_file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/ss_T1w_seg.nii'


**Outputs parameters:**

- *out_file*
    Head mask.

    ::

      ex. '/home/username/data/derived_data/denoise_ss_T1w_enh_grad.nii'

-------------

Usefull links:
`mriqc <https://mriqc.readthedocs.io/en/22.0.6//>`_
