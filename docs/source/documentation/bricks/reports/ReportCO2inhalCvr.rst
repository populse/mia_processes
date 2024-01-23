:orphan:

.. toctree::

+--------------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+-------------------------------------------+----------------------------------------------------+

=======================
ReportCO2inhalCvr brick
=======================

Generates the report for CVR evaluation using hypercapnia challenge pipeline
----------------------------------------------------------------------------


**Inputs parameters:**

- *anat* (a string representing an existing file)
    Input raw anatomical image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *norm_anat* (a string representing an existing file)
    Normalised anatomical image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/wT1w.nii'

- *norm_anat_fig_rows* (an integer, optional, default value is 5)
    The number of columns for the normalised anatomical slice planes plot.

    ::

      ex. '/home/username/data/derived_data/art_T1w_mask.nii'

- *norm_anat_fig_cols* (a pathlike object or string representing a file)
    Air mask ("hat" mask without artifacts) image.

    ::

      ex. '/home/username/data/derived_data/air_T1w_mask.nii'


- *norm_anat_inf_slice_start* (a string representing an existing file)
    Mask image.

    ::

      ex. '/home/username/data/derived_data/pre_n4c_T1W_clipped_desc-brain.nii'

- *norm_anat_slices_gap* (a pathlike object or string representing a file)
    Outskin mask (extensions: [.nii, .nii.gz])

    ::

      ex. '/home/username/data/derived_data/T1w_brain_outskin_mask.nii'

- *norm_func* (a string representing an existing file)
    Segmentation mask (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/ss_n4c_T1w_clippep_seg.nii'

- *norm_func_fig_rows* (an integer, optional, default value is 5)
    The number of lines for the anatomical slice planes plot.

    ::

      ex. 5

- *norm_func_fig_cols* (an integer, optional, default value is 5)
    The number of columns for the anatomical slice planes plot.

    ::

      ex. 5

- *norm_func_inf_slice_start* (a string, optional, default value is 'spikes')
    Prefix of the output image.

    ::

      ex. 'spikes_'

- *norm_func_slices_gap* (a string representing an existing file)
    An existing, uncompressed normalised anatomical image file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *realignment_parameters_desc* (an integer, optional, default value is 5)
    The number of lines for the normalised anatomical slice planes plot.

    ::

      ex. 5

- *patient_info* (an integer, optional, default value is 5)
    The number of columns for the normalised anatomical slice planes plot.

    ::

      ex. 5


**Outputs parameters:**

- *report* (a strings representing a file)
    The generated report (pdf).

    ::

      ex. '/home/username/data/derived_data/T1w_ref_anatomical_mriqcReport_2023_03_31_11_42_10_75.pdf'


-------------

Usefull links:

`mriqc v22.06 <https://mriqc.readthedocs.io/en/22.0.6/>`_
