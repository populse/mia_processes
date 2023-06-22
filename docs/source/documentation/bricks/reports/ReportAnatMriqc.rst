:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=====================
ReportAnatMriqc brick
=====================

Generates the report for anatomical data in mriqc pipeline
----------------------------------------------------------

**Mandatory inputs parameters:**

- *anat* (a string representing an existing file)
    Input anat image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *IQMs_file* (a string representing an existing file)
    A .json file containing the IQMs (valid extensions: [.json]).

    ::

      ex. '/home/username/data/derived_data/T1w_anat_qc.json'

- *art_mask* (a pathlike object or string representing a file)
    Artifact mask image (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/art_T1w_mask.nii'

- *air_mask* (a pathlike object or string representing a file)
    Air mask ("hat" mask without artifacts) image.

    ::

      ex. '/home/username/data/derived_data/air_T1w_mask.nii'


- *brain_mask* (a string representing an existing file)
    Mask image.

    ::

      ex. '/home/username/data/derived_data/pre_n4c_T1W_clipped_desc-brain.nii'

- *head_mask* (a pathlike object or string representing a file)
    Outskin mask (extensions: [.nii, .nii.gz])

    ::

      ex. '/home/username/data/derived_data/T1w_brain_outskin_mask.nii'

- *segmentation* (a string representing an existing file)
    Segmentation mask (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/ss_n4c_T1w_clippep_seg.nii'

**Optional inputs with default value parameters:**

- *anat_fig_cols* (an integer, optional, default value is 5)
    The number of lines for the anatomical slice planes plot.

    ::

      ex. 5

- *anat_fig_rows* (an integer, optional, default value is 5)
    The number of columns for the anatomical slice planes plot.

    ::

      ex. 5

- *out_prefix* (a string, optional, default value is 'spikes')
    Prefix of the output image.

    ::

      ex. 'spikes_'

- *norm_anat* (a string representing an existing file)
    An existing, uncompressed normalised anatomical image file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *norm_anat_fig_cols* (an integer, optional, default value is 5)
    The number of lines for the normalised anatomical slice planes plot.

    ::

      ex. 5

- *norm_anat_fig_rows* (an integer, optional, default value is 5)
    The number of columns for the normalised anatomical slice planes plot.

    ::

      ex. 5

**Optional inputs:**

- *anat_inf_slice_start* (an integer, optional)
    The first index displayed in anatomical slice planes plot.
    Default is Undefined (ie parameter not used).

    ::

      ex. 0

- *anat_slices_gap* (an integer, optional)
    Gap between slices in anatomical slice planes plot.
    Default is Undefined (ie parameter not used).

    ::

      ex. 0


- *norm_anat_inf_slice_start* (an integer, optional)
    The first index displayed in normalised anatomical slice planes plot.
    Default is Undefined (ie parameter not used).

    ::

      ex. 0

- *norm_anat_slices_gap* (an integer, optional)
    Gap between slices in normalised anatomical slice planes plot.
    Default is Undefined (ie parameter not used).

    ::

      ex. 0



**Outputs parameters:**

- *report* (a strings representing a file)
    The generated report (pdf).

    ::

      ex. '/home/username/data/derived_data/T1w_ref_anatomical_mriqcReport_2023_03_31_11_42_10_75.pdf'


-------------

Usefull links:

`mriqc v22.06 <https://mriqc.readthedocs.io/en/22.0.6/>`_
