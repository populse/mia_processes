:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
GradientThreshold brick
============

Computes a threshold from the histogram of the magnitude gradient image using a segmentation of the input image.

Adapted from https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/anatomical.py#L1039

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import GradientThreshold
>>> GradientThreshold.help()

**Inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *seg_file* (a string representing an existing file)
    Segmentation file (valid extensions: [.nii, .nii.gz]). 
    For example, a tissue class map from the brick FastSegment (FSL). 

    ::

      ex. 'home/username/data/derived_data/T1w_seg.nii'

- *prefix* (a string, optional)
    Prefix of the output image. Default is ''.
    
    ::

      ex. ''


- *suffix* (a string, optional)
   Suffix of output image. Default is '_grad'.
    
    ::

      ex. '_grad'


**Outputs parameters:**

- *out_file* (a strings representing a file)
    Path of the thresholded scan (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/T1w_grad.nii'

-------------

Usefull links:
`Gradient threshold mriqc <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/anatomical.py#L1039>`_
