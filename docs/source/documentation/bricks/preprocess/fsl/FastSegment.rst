:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
FastSegment brick
============

Brain tissue segmentation using FAST (FSL). 
This brick can be used to segment one input image. 
The tissue class segmented are Grey Matter, White Matter and CSF.
By default the image type is set to T1 but it is also possible to segment T2 and PD images by changing the img_type parameter. 

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    An image to be segmented. An existing file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. ['/home/username/data/raw_data/T1w.nii']

**Optional inputs with default value parameters:**

- *img_type* (an int, 1, 2 or 3, optional, default value is 1)
    Int specifying type of image: (1 = T1, 2 = T2, 3 = PD). Default image type is T1. 
    This option aids the segmentation in identifying which classes are which tissue type.

    ::

      ex. 1

- *output_type* ('NIFTI' or 'NIFTI_GZ', optional, default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *segments* (a boolean, optional, default value is True)
    Outputs a separate binary image for each tissue type.

    ::

      ex. True

**Outputs parameters:**

- *mixeltype* (a pathlike object or string representing a file)
    The mixeltype file represents the classification of each voxel's tissue mixture.  That is, voxels containing only one tissue type have a different mixeltype from that containing mixtures of two tissues, which is different again from those containing mixtures of all three tissues. 
   
    ::

      ex. '/home/username/data/derived_data/T1w_mixeltype.nii'

- *partial_volume_files* (a list of items which are file names.)
    A (non-binary) partial volume image for each class, where each voxel contains a value in the range 0-1 that represents the proportion of that class's tissue present in that voxel.

    ::

      ex.['/home/username/data/derived_data/T1w_pve_0.nii', '/home/username/data/derived_data/T1w_pve_1.nii', '/home/username/data/derived_data/T1w_pve_2.nii']

- *partial_volume_map* (a pathlike object or string representing a file)
    The pveseg map  contains the best hard segmentation that follows from the partial volume segmentation results. 
    That is, each voxel contains one number, representing the tissue type with the largest partial volume fraction.
    It can be useful for a quick visual assessment of the segmentation.
   
    ::

      ex. '/home/username/data/derived_data/T1w_pveseg.nii'

- *tissue_class_files* (a list of items which are file names.)
    Binary segmented volume files, one image per class. Values are only either 0 or 1.

    ::

      ex. ['/home/username/data/derived_data/T1w_seg_0.nii', '/home/username/data/derived_data/T1w_seg_1.nii', '/home/username/data/derived_data/T1w_seg_2.nii']

- *tissue_class_map* (a pathlike object or string representing a file)
    A binary segmented volume file where each voxel is classified into only one class.
    It is a single image that contains all the necessary information, with the first class taking intensity value 1 in the image, etc. 

    ::

      ex. '/home/username/data/derived_data/T1w_seg.nii'

-------------

Usefull links:
`FSL FAST <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FAST>`_
`FSL FAST - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.fsl.preprocess.html#fast>`_