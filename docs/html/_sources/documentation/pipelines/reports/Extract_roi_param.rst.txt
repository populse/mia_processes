:orphan:

.. toctree::

+--------------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+-------------------------------------------+----------------------------------------------------+

==========================
Extract_roi_param pipeline
==========================

Produces gray matter masks for various ROIs and means, standard deviations, laterality indices for beta and spmT values in these ROIs
-------------------------------------------------------------------------------------------------------------------------------------

**Pipeline insight**

- Extract_roi_param pipeline combines the following bricks:
    - `Concat_to_list_of_list <../../bricks/tools/Concat_to_list_of_list.html>`_
    - `Import_Data <../../bricks/tools/Import_Data.html>`_
    - `Find_In_List <../../bricks/tools/Find_In_List.html>`_
    - `Files_To_List <../../bricks/tools/Files_To_List.html>`_
    - `ConvROI <../../bricks/preprocess/others/ConvROI.html>`_
    - `Resample1 <../../bricks/preprocess/others/Resample1.html>`_
    - `Resample2 <../../bricks/preprocess/others/Resample2.html>`_
    - `Mean_stdDev_calc <../../bricks/reports/Mean_stdDev_calc.html>`_
    - `Concat_to_list <../../bricks/tools/Concat_to_list.html>`_
    - `Result_collector <../../bricks/reports/Result_collector.html>`_

.. image:: ../../images/Extract_roi_param.png
  :width: 1000
  :alt: spatial mask pipeline

--------------------------------------

**Inputs parameters:**

- *spmT_images*
   A list of T-statistics images, previously obtained from the `EstimateContrast <../../bricks/stats/spm/EstimateContrast.html>`_ brick.

    ::

      ex. ['/home/username/data/raw_data/spmT_0001.nii']

- *beta_images*
    A list of estimated regression coefficients images (beta_000k.nii, where k indexes the kth regression coefficient), previously obtained from the `EstimateModel <../../bricks/stats/spm/EstimateModel.html>`_ brick.

    ::

      ex. ['/home/username/data/raw_data/beta_0001.nii',
           '/home/username/data/raw_data/beta_0002.nii',
           '/home/username/data/raw_data/beta_0003.nii',
           '/home/username/data/raw_data/beta_0004.nii',
           '/home/username/data/raw_data/beta_0005.nii',
           '/home/username/data/raw_data/beta_0006.nii',
           '/home/username/data/raw_data/beta_0007.nii',
           '/home/username/data/raw_data/beta_0008.nii']

- *mask_002*
    A grey matter mask with a resolution defined previously in the `Spatial_mask <../preprocess/Spatial_mask.html>`_.

    ::

      ex. /home/username/data/raw_data/mask_anat_002.nii

- *patient_info*
    A dictionary whose keys/values correspond to information about the patient. If the value for a key is not defined, the corresponding tag and its value will be searched for in the database.
    `patient_info` keys/values therefore take precedence over database tags, allowing the user to force patient data.

    ::

      ex. {'PatientName': <undefined>,
           'Pathology': <undefined>,
	   'Age': <undefined>,
	   'Sex': <undefined>,
	   'MR': <undefined>,
	   'Gas': <undefined>,
	   'GasAdmin': <undefined>}

**Outputs parameters:**

- *resample2_masks*
    A list of images, corresponding to `conv_roi_maks` after resampling to the resolution of `spmT_images`.

    ::

      ex. ['/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convACA_L_2.nii',
           '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convACA_R_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convACM_L_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convACM_R_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convACP_L_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convACP_R_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convPICA_L_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convPICA_R_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-CING_L_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-CING_R_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-FRON_L_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-FRON_R_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-INSULA_L_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-INSULA_R_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-OCC_L_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-OCC_R_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-PAR_L_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-PAR_R_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-STR_L_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-STR_R_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-TEMP_L_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-TEMP_R_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-THA_L_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convROI-THA_R_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convSCA_L_2.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convSCA_R_2.nii']

- *xls_files*
    A list of xls files containing the means and standard deviations of the parameters in the ROIs defined by the `Concat_to_list_of_list <../../bricks/tools/Concat_to_list_of_list.html>`_ brick.
    The laterality index between the two hemispheres is also calculated (files with ``IL`` in the name).

    ::

      ex. ['/home/username/data/derived_data/patient-name_data/results_aggregation/BOLD_IL_mean_beta.xls',
           '/home/username/data/derived_data/patient-name_data/results_aggregation/BOLD_mean_spmT.xls',
	   '/home/username/data/derived_data/patient-name_data/results_aggregation/BOLD_mean_beta.xls',
	   '/home/username/data/derived_data/patient-name_data/results_aggregation/BOLD_IL_std_spmT.xls',
	   '/home/username/data/derived_data/patient-name_data/results_aggregation/BOLD_IL_mean_spmT.xls',
	   '/home/username/data/derived_data/patient-name_data/results_aggregation/BOLD_std_spmT.xls',
	   '/home/username/data/derived_data/patient-name_data/results_aggregation/BOLD_std_beta.xls',
	   '/home/username/data/derived_data/patient-name_data/results_aggregation/BOLD_IL_std_beta.xls']

- *conv_roi_maks*
    A list of images, resulting from the convolution of `mask_002` with the ROIs defined by the `Concat_to_list_of_list <../../bricks/tools/Concat_to_list_of_list.html>`_ brick.

    ::

      ex. ['/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convACA_L.nii',
           '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convACA_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convACM_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_dat/convROI_BOLD/convACM_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convACP_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convACP_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convPICA_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convPICA_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-CING_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-CING_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-FRON_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-FRON_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-INSULA_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-INSULA_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-OCC_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-OCC_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-PAR_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-PAR_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-STR_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-STR_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-TEMP_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-TEMP_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-THA_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convROI-THA_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convSCA_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convSCA_R.nii']
