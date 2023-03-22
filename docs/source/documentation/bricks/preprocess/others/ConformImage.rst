:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
ConformImage brick
============

Conform T1w image: 
- remove axes length 1 at end of image shapes if needed (nibabel.squeeze_image)
- reordered to be closest to canonical if needed (nibabel.as_closest_canonical)

Adapted from https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/common/conform_image.py#L75

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import ConformImage
>>> ConformImage.help()

**Inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz])

    ::

      ex. ['/home/username/data/raw_data/T1w.nii']


- *prefix* (a string)
    Prefix of the output image. Default is ''.
    
    ::

      ex. ''


- *suffix* (a string)
   Suffix of output image. Default is ''.
    
    ::

      ex. ''


**Outputs parameters:**

- *out_file* (a list of pathlike objects or strings representing a file)
    Conforme Image (extensions: [.nii, .nii.gz])
    
    ::

      ex. ['/home/username/data/raw_data/T1w_bin.nii']

-------------

Usefull links:
`Conform image mriqc <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/common/conform_image.py#L75>`_
`nibabel functions <https://nipy.org/nibabel/reference/nibabel.funcs.html>`_
