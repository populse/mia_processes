:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===================
IntensityClip brick
===================

Clip the intensity range as prescribed by the percentiles
Remove outliers at both ends of the intensity distribution and fit into a given dtype.
To more robustly determine the clipping thresholds, data are removed of spikes with a median filter.
Once the thresholds are calculated, the denoised data are thrown away and the thresholds are applied on the original image.

Adapted from `niworkflows <https://github.com/nipreps/niworkflows>`_.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image which intensity will be clipped (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

**Optional inputs with default value parameters:**

- *dtype* (one of int16, float32, uint8, optional, default value is int16)
    Output datatype.

    ::

      ex. int16

- *invert* (a boolean, optional, default value is False)
    Finalize by inverting contrast.

    ::

      ex. False

- *nonnegative* (a boolean, optional, default value is True)
    Whether input intensities must be positive.

    ::

      ex. True

- *p_max* (a float,  optional, default value is 99.9)
    Percentile for the upper bound.

    ::

      ex. 99.9

- *p_min* (a float, optional, default value is 10.0)
    Percentile for the lower bound.

    ::

      ex. 10.0

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Path of the clipped image (xtensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/T1w_clipped.nii'

-------------

Usefull links:
`IntensityClip niworflow <https://github.com/nipreps/niworkflows/blob/c2b8e0f7396c626df296a48217655298a82f3069/niworkflows/interfaces/nibabel.py#L460>`_
