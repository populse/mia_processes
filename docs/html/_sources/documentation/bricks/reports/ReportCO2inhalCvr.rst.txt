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

      ex. '/home/username/data/raw_data/Anat.nii'

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

- *realignment_parameters_desc* (a string representing an existing file)
    Estimated translation and rotation parameters (valid extensions: .txt).

    ::

      ex. 5

- *patient_info* (a dictionary with keys in ['PatientRef', 'Pathology', 'Age', 'Sex', MR', 'Gas', 'GasAdmin'])
    Optional dictionary with information about the patient.

    ::

      ex. {'PatientRef': 'sub-1', 'Pathology': 'ACMD', 'Age': '56', 'Sex': 'M', 'MR': '3T', 'Gas': 'Bactal', 'GasAdmin': 'Mask'}


**Outputs parameters:**

- *report* (a strings representing a file)
    The generated report (pdf).

    ::

      ex. '/home/username/data/derived_data/sub-1_CO2_inhal_CVR_Report_2024_01_24_09_34_58_08.pdf'
