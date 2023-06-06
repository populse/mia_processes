:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

================
Resample2 brick
================

Setting regions of interest to the resolution of the *reference_image* using skimage.transform.resize().
--------------------------------------------------------------------------------------------------------

    - ROIs are defined from `doublet_list` parameter as
      doublet_list[0][0] + doublet_list[0][1] + '.nii',
      doublet_list[1][0] + doublet_list[1][1] + '.nii',
      etc.
    - The "/roi\_ `PatientName`/convROI_BOLD2" directory is
      created to receive the convolution results from the runtime (e.g.
      "/roi\_ `PatientName`/convROI_BOLD2/conv"doublet_list[0][0] + doublet_list[0][1]+ "2.nii"
    - To work correctly, the  ROIs must be in the  "/roi\_ `PatientName`/convROI_BOLD
      directory (e.g.
      "/roi\_ `PatientName`/convROI_BOLD/conv"doublet_list[0][0] + doublet_list[0][1]+ ".nii").
    - To work correctly, the database entry for the reference_image parameter
      must have the `PatientName` tag filled in.

--------------------------------------

**Mandatory inputs parameters:**

- *doublet_list* (a list of lists)
    A list of lists containing doublets of strings.

    ::

      ex. [["ROI_OCC", "_L"], ["ROI_OCC", "_R"], ["ROI_PAR", "_L"], ['ROI_PAR', '_R']]


- *reference_image* (a string representing an existing file)
    The reference image for resampling (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/template.nii'


**Outputs parameters:**

- *out_images* (list of pathlike object or string reprsenting a file or a list of items which are a pathlike object or string representing a file)
    The resampled  images ( extensions: [.nii, .nii.gz]).

    ::

      ex. ['/home/username/data/derived_data/roi_PatientName/convROI_BOLD2/convROI_OCC_L2.nii',
           '/home/username/data/derived_data/roi_PatientName/convROI_BOLD2/convROI_OCC_R2.nii',
	   '/home/username/data/derived_data/roi_PatientName/convROI_BOLD2/convROI_PAR_L2.nii',
	   '/home/username/data/derived_data/roi_PatientName/convROI_BOLD2/convROI_PAR_R2.nii']

-------------

Usefull links:
 `skimage.transform.resize() <https://scikit-image.org/docs/stable/api/skimage.transform.html#resize>`_
