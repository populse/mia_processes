:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

====================
DWIBiasCorrect brick
====================

Perform B1 field inhomogeneity correction for a DWI image
----------------------------------------------------------

(mrtrix dwibiascorrect command)

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input DWI image (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/DWI.mif'

**Optional inputs with default value parameters:**

- *use_fsl* (a boolean, default value is False, optional)
    Use FSL FAST to estimate the inhomogeneity field.
    FSL should be configured in populse.

    ::

      ex. False

- *use_ans* (a boolean, default value is True, optional)
    Use ANTS N4 to estimate the inhomogeneity field.
    ANTS should be configured in populse.

    ::

      ex. True

**Optional inputs parameters:**

- *in_mask* (a string representing an existing file, optional)
    Input mask image for bias field estimation (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/DWI_brainmask.nii'

**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    The output bias corrected DWI image

    ::

      ex. '/home/username/data/derived_data/DWI_unbias.mif'

- *bias_field_map* (a pathlike object or string representing a file)
    The output bias field map

    ::

      ex. '/home/username/data/derived_data/DWI_bias_field_map.mif'

-------------

Usefull links:

`mrtrix dwibiascorrect <https://mrtrix.readthedocs.io/en/latest/reference/commands/dwibiascorrect.html>`_

`mrtrix dwibiascorrect - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.preprocess.html#dwibiascorrect>`_
