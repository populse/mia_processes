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

    - The `parametric_maps` are first convolved with the ROIs corresponding
      to `doublet_list`. The mean and standard deviation are then calculated
      for the ROI in each parametric maps.
    - ROIs are defined from `doublet_list` parameter as
      doublet_list[0][0] + doublet_list[0][1] + '.nii',
      doublet_list[1][0] + doublet_list[1][1] + '.nii',
      etc.
    - The "/roi\_\ `PatientName`/ROI_analysis" directory is
      created to receive the results from the runtime (e.g.
      "/roi\_\ `PatientName`/ROI_analysis/" + doublet_list[0][0] + doublet_list[0][1] + "_mean_spmT_BOLD.txt").
    - To work correctly, the "/roi\_/`PatientName`/convROI_BOLD" directory
      must exist and contain the ROI files (e.g.
      "/roi\_\ `PatientName`/convROI_BOLD/conv" + doublet_list[0][0] + doublet_list[0][1] + ".nii").
    - To work correctly, the database entry for the first element of `parametric_maps` must
      have the `PatientName` tag filled in.

--------------------------------------

**Mandatory inputs parameters:**

- *doublet_list* (a list of lists)
    A list of lists containing doublets of strings.

    ::

      ex. [["ROI_OCC", "_L"], ["ROI_OCC", "_R"], ["ROI_PAR", "_l"], ["ROI_PAR", "_R"]]


- *parametric_maps* (a list of existing files)

    ::

      ex. ['/home/username/data/raw_data/spmT_0001.nii',
           '/home/username/data/raw_data/beta_0001.nii']


**Outputs parameters:**

- *mean_out_files* (a list of files)
    A list of .txt files with the calculated average for each ROI determined after convolution.

    ::

      ex. ['/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_OCC_L_mean_spmT_BOLD.txt',
           '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_OCC_R_mean_spmT_BOLD.txt',
	   '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_PAR_L_mean_spmT_BOLD.txt',
	   '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_PAR_R_mean_spmT_BOLD.txt',
	   '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_OCC_L_mean_beta_BOLD.txt',
	   '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_OCC_R_mean_beta_BOLD.txt',
	   '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_PAR_L_mean_beta_BOLD.txt',
	   '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_PAR_R_mean_beta_BOLD.txt']

- *std_out_files* (a list of files)
    A list of .txt files with the standard deviation for each ROI determined after convolution.

    ::

      ex. ['/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_OCC_L_std_spmT_BOLD.txt',
           '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_OCC_R_std_spmT_BOLD.txt',
	   '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_PAR_L_std_spmT_BOLD.txt',
	   '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_PAR_R_std_spmT_BOLD.txt',
	   '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_OCC_L_std_beta_BOLD.txt',
	   '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_OCC_R_std_beta_BOLD.txt',
	   '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_PAR_L_std_beta_BOLD.txt',
	   '/home/username/data/derived_data/roi_PatientName/ROI_analysis/ROI_PAR_R_std_beta_BOLD.txt']
