:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

======================
Mean_stdDev_calc brick
======================

Makes the mean and standard deviation of parametric maps.
---------------------------------------------------------

    - The `parametric_maps` are first resized, if necessary, to the size of the `rois_files`.
    - Next, the `parametric_maps` and the `rois_files` are convolved.
    - Finally, the mean and standard deviation are calculated for the corresponding ROIs.
    - The “PatientName_data/ROI_data/ROI_analysis” directory is created to receive the results.
      If this directory exists at runtime, it is overwritten.
    - Output file names are built like this:
        - ``roi``\_ ``calculation``\_ ``parameter``\_ ``contrast``.txt
	    - ``roi`` is deducted from each `rois_files` after deleting the extension. If `prefix_to_delete`
	      is defined and if it corresponds to the beginning of ``roi``, this beginning of string is deleted.
	    - ``calculation`` corresponds to "mean" (mean calculation) or "std" (standard deviation calculation).
	    - ``parameter`` is deducted from each `parametric_maps` file name. This is the string before the
	      first underscore. If there is no underscore, this is the file name after removing the extension.
	    - ``contrast`` is `contrast_type`.
    - To work correctly, the database entry for the first element of `parametric_maps` must have the
      `PatientName` tag filled in.

--------------------------------------

**Mandatory inputs parameters:**

- *parametric_maps*
    A list of uncompressed file.

    ::

      ex. ['/home/username/data/raw_data/spmT_0001.nii',
           '/home/username/data/raw_data/beta_0001.nii']

- *rois_files*
    A list of regions of interest (a list of uncompressed file), which will be applied to the parametric maps
    before calculating the mean and standard deviation of the parameters in the corresponding regions.

    ::

      ex. ['/home/username/data/raw_data/convACA_L.nii',
           '/home/username/data/raw_data/convACA_R.nii',
           '/home/username/data/raw_data/convACM_L.nii',
           '/home/username/data/raw_data/convACM_R.nii']

- *contrast_type*
    The contrast used (a string).

    ::

      ex. BOLD

- *prefix_to_delete*
    The string to be deleted from the deduced ROI name, when creating the `mean_out_files` and the `std_out_files`.

    ::

      ex. conv

**Outputs parameters:**

- *mean_out_files*
    A list of .txt files with the calculated mean for each ROI convolved with each parametric map.

    ::

      ex. ['/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACA_L_mean_spmT_BOLD.txt',
           '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACA_R_mean_spmT_BOLD.txt',
           '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACM_L_mean_spmT_BOLD.txt',
           '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACM_R_mean_spmT_BOLD.txt',
           '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACA_L_mean_beta_BOLD.txt',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACA_R_mean_beta_BOLD.txt',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACM_L_mean_beta_BOLD.txt',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACM_R_mean_beta_BOLD.txt']


- *std_out_files*
    A list of .txt files with the calculated standard deviation for each ROI convolved with each parametric map.

    ::

      ex. ['/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACA_L_std_spmT_BOLD.txt',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACA_R_std_spmT_BOLD.txt',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACM_L_std_spmT_BOLD.txt',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACM_R_std_spmT_BOLD.txt'
	   '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACA_L_std_beta_BOLD.txt',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACA_R_std_beta_BOLD.txt',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACM_L_std_beta_BOLD.txt',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/ROI_analysis/ACM_R_std_beta_BOLD.txt']
