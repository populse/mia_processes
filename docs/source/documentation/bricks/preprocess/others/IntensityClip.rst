:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
IntensityClip brick
============

Clip the intensity range as prescribed by the percentiles
Remove outliers at both ends of the intensity distribution and fit into a given dtype.
To more robustly determine the clipping thresholds, data are removed of spikes with a median filter.
Once the thresholds are calculated, the denoised data are thrown away and the thresholds are applied on the original image.

Adapted from https://github.com/nipreps/niworkflows/blob/c2b8e0f7396c626df296a48217655298a82f3069/niworkflows/interfaces/nibabel.py#L460 

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import IntensityClip

>>> IntensityClip.help()

**Inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image which intensity will be clipped (valid extensions: [.nii, .nii.gz])

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *p_min* (a float, optional=True)
    Percentile for the lower bound. Default is 10.0. 

    ::

      ex. 10.0

- *p_max* (a float,  optional=True)
    Percentile for the upper bound. Default is 99.9

    ::

      ex. 99.9

- *nonnegative* (a boolean, optional=True)
    Whether input intensities must be positive. Default is True.
    
    ::

      ex. True

- *dtype* (one of 'int16', 'float32', 'uint8', optional=True)
    Output datatype. Default is 'int16'.
    
    ::

      ex. 'int16'

- *invert* (a boolean, optional=True)
   Finalize by inverting contrast. Default is False.
    
    ::

      ex. False



**Outputs parameters:**

- *out_file* (a strings representing a file)
    Path of the clipped image 
    
    ::

      ex. '/home/username/data/raw_data/T1w_clipped.nii'

-------------

Usefull links:
`IntensityClip niworflow <https://github.com/nipreps/niworkflows/blob/c2b8e0f7396c626df296a48217655298a82f3069/niworkflows/interfaces/nibabel.py#L460>`_
