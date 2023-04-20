:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

========================
CarpetParcellation brick
========================

Dilate brainmask, substract from itself then generate the union of obtained crown mask and epi parcellation.

Adapted from `mriqc <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/functional.py#L1022>`_ ,
`niworkflow Binary dilatation <https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/morphology.py#L46>`_,
`niworkflow Binary Subtraction <ttps://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/morphology.py#L79>`_.

--------------------------------------

**Mandatory inputs parameters:**

- *brain_mask* (a string representing an existing file)
    Brain mask (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/automask_mean_reg_func_valid.nii'

- *segmentation* (a string representing an existing file)
    Template carpet register in subject space (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/t_tpl-MNI152NLin2009cAsym_res-01_desc-carpet_dseg.nii.gz'


**Optional inputs with default value parameters:**

- *out_prefix* (a string, optional, default value is 'cseg')
    Specify the string to be prepended to the segmentation filename of the output file.

    ::

      ex. 'cseg_'

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Carpet parcelattion file (extensions: [.nii, .nii.gz])

    ::

      ex. '/home/username/data/derived_data/cseg_t_tpl-MNI152NLin2009cAsym_res-01_desc-carpet_dseg.nii.gz'

-------------

Usefull links:

`niworkflow Binary dilatation <https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/morphology.py#L46>`_,
`niworkflow Binary Subtraction <ttps://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/morphology.py#L79>`_
