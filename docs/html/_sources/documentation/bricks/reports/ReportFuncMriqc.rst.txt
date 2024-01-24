:orphan:

.. toctree::

+--------------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+-------------------------------------------+----------------------------------------------------+

=====================
ReportFuncMriqc brick
=====================

Generates the report for functional data in mriqc pipeline
----------------------------------------------------------

**Mandatory inputs parameters:**

- *func* (a string representing an existing file)
    Input functional image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

- *IQMs_file* (a string representing an existing file)
    A .json file containing the IQMs (valid extensions: [.json]).

    ::

      ex. '/home/username/data/derived_data/func_func_qc.json'

**Optional inputs with default value parameters:**

- *func_fig_cols* (an integer, optional, default value is 5)
    The number of columns for the functional slice planes plot.

    ::

      ex. 5

- *func_fig_rows* (an integer, optional, default value is 5)
    The number of lines for the functional slice planes plot.

    ::

      ex. 5

- *out_prefix* (a string, optional, default value is 'spikes')
    Prefix of the output image.

    ::

      ex. 'spikes_'

- *norm_func_fig_cols* (an integer, optional, default value is 5)
    The number of columns for the normalised functional slice planes plot.

    ::

      ex. 5

- *norm_func_fig_rows* (an integer, optional, default value is 5)
    The number of lines for the normalised functional slice planes plot.

    ::

      ex. 5

**Optional inputs:**

- *brain_mask* (a string representing an existing file)
    Mask image.

    ::

      ex. '/home/username/data/derived_data/automask_mean_func_valid.nii'

- *func_inf_slice_start* (an integer, optional)
    The first index displayed in functional slice planes plot.
    Default is Undefined (ie parameter not used).

    ::

      ex. 0

- *func_slices_gap* (an integer, optional)
    Gap between slices in functional slice planes plot.
    Default is Undefined (ie parameter not used).

    ::

      ex. 0

- *IQMs_plot* (a string representing an existing file)
    A figure with carpet and outliers/dvars/FD/spikes plot.

    ::

      ex. '/home/username/data/raw_data/func_valid_fmriplot.png'

- *norm_func* (a string representing an existing file)
    An existing, normalised functional image file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/w_n4c_mean_reg_func_valid.nii'

- *mean_func* (a string representing an existing file)
    Mean functional image file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/mean_reg_func_valid.nii'


- *norm_func_inf_slice_start* (an integer, optional)
    The first index displayed in normalised functional slice planes plot.
    Default is Undefined (ie parameter not used).

    ::

      ex. 0

- *norm_func_slices_gap* (an integer, optional)
    Gap between slices in normalised functional slice planes plot.
    Default is Undefined (ie parameter not used).

    ::

      ex. 0

- *stddev_func* (a string representing an existing file)
    Functional standard deviation image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/reg_func_valid_stddev.nii'



**Outputs parameters:**

- *report* (a strings representing a file)
    The generated report (pdf).

    ::

      ex. ''/home/username/data/derived_data/func_ref_functional_mriqcReport_2023_03_31_11_42_10_75.pdf'


-------------

Usefull links:

`mriqc v22.06 <https://mriqc.readthedocs.io/en/22.0.6/>`_
