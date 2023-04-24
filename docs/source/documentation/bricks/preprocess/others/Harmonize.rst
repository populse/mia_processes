:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===============
Harmonize brick
===============

Harmonize input image using a white matter mask as decribed in the method to compute the arfifact mask described in [Mortamet2009].
Be carrefull, if the suffix and prefix input parameters are not defined or consist only of one or more white spaces, the input parameter will be overwritten.

Adapted from `mriqc <https://mriqc.readthedocs.io/en/22.0.6//>`_ anatomical workflow.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *wm_mask* (a string representing an existing file)
    Whithe matter mask file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. 'home/username/data/derived_data/T1w_seg.nii'

**Optional inputs with default value parameters:**

- *erodemask* (a boolean, optional, default value is True)
    Erode mask.

    ::

      ex. True

- *prefix* (a string, optional, default value is '')
    Prefix of the output image.

    ::

      ex. ''


- *suffix* (a string, optional, default value is '_harmonized')
    Suffix of output image.

    ::

      ex. '_harmonized'


**Outputs parameters:**

- *out_file* (a strings representing a file)
    Path of the harmonized image (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/T1w_harmonized.nii'

-------------

Usefull links:
`Harmonize mriqc <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/anatomical.py#L405>`_
