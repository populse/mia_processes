:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

================
Resample2 brick
================

Sets images to the resolution of a reference image using using nilearn.image.resample_to_img()
----------------------------------------------------------------------------------------------

- The "`PatientName`\_data/ROI_data/convROI_BOLD2" directory is created to receive the resampling
  results from the runtime. If this directory exists at runtime, it is deleted.
- To work correctly, the database entry for the `reference_image` parameter
  must have the `PatientName` tag filled in.

--------------------------------------

**Mandatory inputs parameters:**

- *files_to_resample*
    The images that will be resampled (a list of pathlike object or string representing a file or a list of items
    which are a pathlike object or string representing a file ; valid extensions: [.img, .nii, .hdr]).

    ::

      ex. ['/home/username/data/raw_data/convACA_L.nii',
           '/home/username/data/raw_data/convACA_R.nii',
	   '/home/username/data/raw_data/convACM_L.nii',
	   '/home/username/data/raw_data/convACM_R.nii']

- *reference_image*
    The reference image for resampling (a string representing an existing file ; valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/template.nii'

- *suffix*
    Suffix of output images (a string, default value is '_2').
    ::

      ex. _2


**Outputs parameters:**

- *out_images*
    The resampled  images (a list of pathlike object or string representing a file or a list of items which are a
    pathlike object or string representing a file ; valid  extensions: [.nii, .nii.gz]).

    ::

      ex. ['/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convACA_L_2.nii',
           '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convACA_R_2.nii',
           '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convACM_L_2.nii',
           '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD2/convACM_R_2.nii']

-------------

Usefull links:

 `nilearn.image.resample_to_img() <https://nilearn.github.io/dev/modules/generated/nilearn.image.resample_to_img.html>`_
