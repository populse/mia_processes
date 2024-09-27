:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

====================
Generate5ttfsl brick
====================

Generate a five-tissues-type (5TT) segmented tissue image with FSL
-------------------------------------------------------------------

The generated 5TT image is suitable for use in Anatomically-Constrained Tractography (ACT).

(mrtrix 5ttgen fsl command)

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input T1w image (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.mif'

**Optional inputs with default value parameters:**

- *no_crop* (a boolean, default value is False, optional)
    Do NOT crop the resulting 5TT image to reduce its size (keep the same dimensions as the input image)

    ::

      ex. False


- *sgm_amyg_hipp* (a boolean, default value is False, optional)
    Represent the amygdalae and hippocampi as sub-cortical grey matter in the 5TT image

    ::

      ex. False

- *premasked* (a boolean, default value is False, optional)
    Indicate that brain masking has already been applied to the input image

    ::

      ex. False


**Optional inputs parameters:**

- *in_mask* (a string representing an existing file, optional)
    Manually provide a brain mask, rather than deriving one in the script (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/brainmask.mif'

- *t2_image* (a string representing an existing file, optional)
    Provide a T2-weighted image in addition to the default T1-weighted image;
    this will be used as a second input to FSL FAST(valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/T2.nii'

**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    The output 5TT image.

    ::

      ex. '/home/username/data/derived_data/T1w_5tt.mif'

-------------

Useful links:

`mrtrix 5ttgen <https://mrtrix.readthedocs.io/en/latest/reference/commands/5ttgen.html>`_

`mrtrix 5ttgen - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.utils.html#generate5tt>`_
