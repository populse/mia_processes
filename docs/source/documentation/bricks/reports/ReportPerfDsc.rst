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


- *CBV_image* (a string representing an existing file)
    Cerebral Blood Volume map (a file with .nii format).
    ::

      ex. /home/username/data/derived_data/Func_CBV_deconv.nii

- *CBV_cmap* (a string, optional)
    The name of the colors range that are used for the CBV map slice
    planes plot (default: `rainbow`).

    ::

      ex. rainbow

- *CBV_vmin* (a float, optional)
    Specifies the minimum value of the data range to be mapped to the
    colormap. Any value lower than vmin will be clipped to vmin before
    mapping to colors.

    ::

      ex. <undefined>

- *CBV_vmax* (a float, optional)
    Specifies the maximum value of the data range to be mapped to the
    colormap. Any value higher than vmax will be clipped to vmax before
    mapping to colors.

    ::

      ex. <undefined>

- *CBF_image* (a string representing an existing file)
    Cerebral Blood Flow map (a file with .nii format).
    ::

      ex. /home/username/data/derived_data/Func_CBF_deconv.nii

- *CBF_cmap* (a string, optional)
    The name of the colors range that are used for the CBF map slice
    planes plot (default: `rainbow`).

    ::

      ex. rainbow

- *CBF_vmin* (a float, optional)
    Specifies the minimum value of the data range to be mapped to the
    colormap. Any value lower than vmin will be clipped to vmin before
    mapping to colors.

    ::

      ex. <undefined>

- *CBF_vmax* (a float, optional)
    Specifies the maximum value of the data range to be mapped to the
    colormap. Any value higher than vmax will be clipped to vmax before
    mapping to colors.

    ::

      ex. <undefined>

- *Tmax_image* (a string representing an existing file)
    The Time to Maximum map (a file with .nii format).
    ::

      ex. /home/username/data/derived_data/Func_Tmax_deconv.nii

- *Tmax_cmap* (a string, optional)
    The name of the colors range that are used for the Tmax map slice
    planes plot (default: `rainbow`).

    ::

      ex. rainbow

- *Tmax_vmin* (a float, optional)
    Specifies the minimum value of the data range to be mapped to the
    colormap. Any value lower than vmin will be clipped to vmin before
    mapping to colors.

    ::

      ex. <undefined>

- *Tmax_vmax* (a float, optional)
    Specifies the maximum value of the data range to be mapped to the
    colormap. Any value higher than vmax will be clipped to vmax before
    mapping to colors.

    ::

      ex. <undefined>

- *MTT_image* (a string representing an existing file)
    Mean Transit Time map (a file with .nii format).
    ::

      ex. /home/username/data/derived_data/Func_MTT_deconv.nii

- *MTT_cmap* (a string, optional)
    The name of the colors range that are used for the MTT map slice
    planes plot (default: `rainbow`).

    ::

      ex. rainbow

- *MTT_vmin* (a float, optional)
    Specifies the minimum value of the data range to be mapped to the
    colormap. Any value lower than vmin will be clipped to vmin before
    mapping to colors.

    ::

      ex. <undefined>

- *MTT_vmax* (a float, optional)
    Specifies the maximum value of the data range to be mapped to the
    colormap. Any value higher than vmax will be clipped to vmax before
    mapping to colors.

    ::

      ex. <undefined>

- *aif_file*
    The AIF (a file in .json format) typically obtained previously with the
    Make_AIF brick.

    ::

      ex. '/home/username/data/raw_data/swrfunc_aif.json'

- *patient_info*
    Optional dictionary with information about the patient (keys in
    ['PatientRef', 'Pathology', 'Age', 'Sex', MR'])

    ::

      ex. {'PatientRef': 'sub-1', 'Pathology': 'ACMD', 'Age': '56', 'Sex': 'M', 'MR': '3T'}

**Outputs parameters:**

- *report* (a strings representing a file)
    The generated report (pdf).

    ::

      ex. /home/username/data/derived_data/sub-1_CO2_inhal_CVR_Report_2024_01_24_09_34_58_08.pdf
