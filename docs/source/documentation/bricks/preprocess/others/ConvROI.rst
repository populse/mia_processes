:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
ConvROI brick
============

Convolve regions of interest with a mask.

    - Resampling the mask to the size of the ROIs, using the first ROI.
    - Then convolve each ROI with resized mask.
    - ROIs are defined from doublet_list parameter as
      doublet_list[0][0] + doublet_list[0][1] + '.nii',
      doublet_list[1][0] + doublet_list[1][1] + '.nii',
      etc.

To work correctly, the database entry for the in_image parameter must have the "PatientName" tag filled in.
This tag is used to create the folder `output_directory/roi_PatientName` (that will contains the ROI files) and the folder `output_directory/roi_PatientName/convROI_BOLD` (that will contains the convolution results)


--------------------------------------

>>> from mia_processes.bricks.preprocess.others import ConvROI
>>> ConvROI.help()

**Inputs parameters:**

- *doublet_list* (a list of lists)
    A list of lists containing doublets of strings. 
    It is used to define ROI names as:
      doublet_list[0][0] + doublet_list[0][1] + '.nii',
      doublet_list[1][0] + doublet_list[1][1] + '.nii',
      ect

    ::

      ex.[["ROI_OCC", "_L"], ["ROI_OCC", "_R"], ["ROI_PAR", "_l"]]


- *in_image* (a string representing an existing file)
    Input mask 
    
    ::

      ex. '/home/username/data/raw_data/mask.nii'



**Outputs parameters:**

- *out_images* (a list of pathlike objects or strings representing a file)
    Conforme Image (extensions: [.nii, .nii.gz])
    
    ::

      ex. ['/home/username/data/raw_data/func.nii']

-------------

