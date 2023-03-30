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

Adapted from `niworkflows <https://github.com/nipreps/niworkflows>`_.

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import IntensityClip

>>> IntensityClip.help()

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image which intensity will be clipped (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

**Optional inputs with default value parameters:**

- *dtype* (one of 'int16', 'float32', 'uint8', optional)
    Output datatype.
    
    ::

      default value. 'int16'

- *invert* (a boolean, optional)
    Finalize by inverting contrast.
    
    ::

      default value. False
      
- *nonnegative* (a boolean, optional=True)
    Whether input intensities must be positive. 
    
    ::

      default value. True

- *p_max* (a float,  optional)
    Percentile for the upper bound.

    ::

      default value. 99.9

- *p_min* (a float, optional)
    Percentile for the lower bound.

    ::

      default value. 10.0

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Path of the clipped image (xtensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/T1w_clipped.nii'

-------------

Usefull links:
`IntensityClip niworflow <https://github.com/nipreps/niworkflows/blob/c2b8e0f7396c626df296a48217655298a82f3069/niworkflows/interfaces/nibabel.py#L460>`_
