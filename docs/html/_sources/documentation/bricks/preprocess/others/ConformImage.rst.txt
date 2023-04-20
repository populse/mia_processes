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
    Prefix of the output image. Default is ''.

    ::

      ex. ''


- *suffix* (a string, optional, default value is '')
    Suffix of output image. Default is ''.

    ::

      ex. ''


**Outputs parameters:**

- *out_file* (a strings representing a file)
  Conforme Image (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/T1w.nii'

-------------

Usefull links:
`Conform image mriqc <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/common/conform_image.py#L75>`_
`nibabel functions <https://nipy.org/nibabel/reference/nibabel.funcs.html>`_
