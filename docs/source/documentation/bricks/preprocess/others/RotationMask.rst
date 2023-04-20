:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
RotationMask brick
============

Compute the rotation mask image as decribed in the method to compute the arfifact mask described in [Mortamet2009].

Be carrefull, if the suffix and prefix input parameters are not defined or consist only of one or more white spaces, the input parameter will be overwritten.

Adapted from `mriqc <https://mriqc.readthedocs.io/en/22.0.6//>`_ anatomical workflow.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

**Optional inputs with default value parameters:**

- *prefix* (a string, optional, default value is '')
    Prefix of the output image.

    ::

      ex. ''


- *suffix* (a string, optional, default value is '_rotmasked')
    Suffix of output image.

    ::

      ex. '_rotmasked'


**Outputs parameters:**

- *out_file* (a strings representing a file)
    Path of the harmonized image (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/T1w_rotmasked.nii'

-------------

Usefull links:
`RotationMask mriqc <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/anatomical.py#L448>`_
