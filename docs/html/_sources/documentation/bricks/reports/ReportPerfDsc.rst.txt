:orphan:

.. toctree::

+--------------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+-------------------------------------------+----------------------------------------------------+

===================
ReportPerfDsc brick
===================

Generates report for Perfusion study using DSC-MRI
--------------------------------------------------


**Inputs parameters:**

- *norm_anat* (a string representing an existing file)
    Normalised anatomical image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/wAnat.nii'

- *norm_anat_fig_rows* (an integer, optional, default value is 5)
    The number of lines for the normalised anatomical slice planes plot.

    ::

      ex. 5

- *norm_anat_fig_cols* (an integer, optional, default value is 5)
    The number of columns for the normalised anatomical slice planes plot.

    ::

      ex. 5


- *norm_anat_inf_slice_start* (an integer, optional)
    The first index displayed in the normalised anatomical slice planes plot.
    Default is `Undefined` (ie automatically determined).

    ::

      ex. 10

- *norm_anat_slices_gap* (an integer, optional)
    Gap between slices in normalised anatomical slice planes plot.
    Default is `Undefined` (ie automatically determined).

    ::

      ex. 13

- *norm_anat_cmap* (a string, optional)
    The name of the colors range that are used for the normalised anatomical
    slice planes plot (default is `Greys_r`).

    ::

      ex. Greys_r

- *norm_anat_vmin* (a float, optional)
    Specifies the minimum value of the data range to be mapped to the
    colormap. Any value lower than vmin will be clipped to vmin before
    mapping to colors.

    ::

      ex. <undefined>

- *norm_anat_vmax* (a float, optional)
    Specifies the maximum value of the data range to be mapped to the
    colormap. Any value higher than vmax will be clipped to vmax before
    mapping to colors.

    ::

      ex. <undefined>

- *norm_func* (a string representing an existing file)
    Normalised functional image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/wrFunc.nii'

- *norm_func_fig_rows* (an integer, optional, default value is 5)
    The number of lines for the normalised functional slice planes plot.

    ::

      ex. 5

- *norm_func_fig_cols* (an integer, optional, default value is 5)
    The number of columns for the normalised functional slice planes plot.

    ::

      ex. 5

- *norm_func_inf_slice_start* (an integer, optional)
    The first index displayed in the normalised functional slice planes plot.
    Default is `Undefined` (ie automatically determined).

    ::

      ex. <undefined>

- *norm_func_slices_gap* (an integer, optional)
    Gap between slices in normalised functional slice planes plot.
    Default is `Undefined` (ie automatically determined).


    ::

      ex. <undefined>

- *norm_func_cmap* (a string, optional)
    The name of the colors range that are used for the normalised functional
    slice planes plot (default is `Greys_r`).

    ::

      ex. Greys_r

- *norm_func_vmin* (a float, optional)
    Specifies the minimum value of the data range to be mapped to the
    colormap. Any value lower than vmin will be clipped to vmin before
    mapping to colors.

    ::

      ex. <undefined>

- *norm_func_vmax* (a float, optional)
    Specifies the maximum value of the data range to be mapped to the
    colormap. Any value higher than vmax will be clipped to vmax before
    mapping to colors.

    ::

      ex. <undefined>

- *realignment_parameters* (a string representing an existing file)
    Estimated translation and rotation parameters during fMRI recording
    (valid extensions: .txt).

    ::

      ex. /home/username/data/derived_data/rp_Func.txt


- *beta_image* (a string representing an existing file)
    The estimated effect size or regression coefficients associated with the
    1st predictor (e.g., conditions or regressors) included in the GLM model.
    A spatial map that represent the estimated effect size or intensity of
    activation/deactivation associated with the 1st predictor/regressor at
    each voxel in the brain (valid extensions:.nii)

    ::

      ex. /home/username/data/derived_data/patient_ref_data/beta_0001.nii

- *beta_cmap* (a string, optional)
    The name of the colors range that are used for the beta map slice
    planes plot (default: `rainbow`).

    ::

      ex. rainbow

- *beta_vmin* (a float, optional)
    Specifies the minimum value of the data range to be mapped to the
    colormap. Any value lower than vmin will be clipped to vmin before
    mapping to colors.

    ::

      ex. <undefined>

- *beta_vmax* (a float, optional)
    Specifies the maximum value of the data range to be mapped to the
    colormap. Any value higher than vmax will be clipped to vmax before
    mapping to colors.

    ::

      ex. <undefined>

- *spmT_image* (a string representing an existing file)
    A file containing t-statistic values per voxel, which indicate the
    strength of the effect of interest at each brain voxel, derived from the
    general linear model (GLM) analysis performed in SPM (valid
    extensions:.nii)

    ::

      ex. /home/username/data/derived_data/patient_ref_data/spmT_0001.nii

- *spmT_cmap* (a string, optional)
    The name of the colors range that are used for the spmT map slice
    planes plot (default: `rainbow`).

    ::

      ex. rainbow

- *spmT_vmin* (a float, optional)
    Specifies the minimum value of the data range to be mapped to the
    colormap. Any value lower than vmin will be clipped to vmin before
    mapping to colors.

    ::

      ex. <undefined>

- *spmT_vmax* (a float, optional)
    Specifies the maximum value of the data range to be mapped to the
    colormap. Any value higher than vmax will be clipped to vmax before
    mapping to colors.

    ::

      ex. <undefined>

- *patient_info*
    Optional dictionary with information about the patient (keys in
    ['PatientRef', 'Pathology', 'Age', 'Sex', MR', 'Gas', 'GasAdmin'])

    ::

      ex. {'PatientRef': 'sub-1', 'Pathology': 'ACMD', 'Age': '56', 'Sex': 'M', 'MR': '3T', 'Gas': 'Bactal', 'GasAdmin': 'Mask'}


**Outputs parameters:**

- *report* (a strings representing a file)
    The generated report (pdf).

    ::

      ex. /home/username/data/derived_data/sub-1_CO2_inhal_CVR_Report_2024_01_24_09_34_58_08.pdf
