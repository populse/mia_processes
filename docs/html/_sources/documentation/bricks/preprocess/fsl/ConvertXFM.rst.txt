:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

================
ConvertXFM brick
================

Modify transformation matrix using convert_xfm (FSL)
-----------------------------------------------------


**Mandatory inputs parameters:**

- *in_transfo* (a string representing an existing file)
    Input transformation matrix (valid extensions: [.mat]).

    ::

      ex. '/home/username/data/raw_data/epi2struct.mat'

**Optional inputs with default value parameters:**

- *output_type* ('NIFTI' or 'NIFTI_GZ', optional, default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   - NIFTI: \*.nii
    |   - NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *invert_xfm* (a boolean, optional, default value is False)
    Invert input transformation (exclusive with fix_scale_skew, concat_xfm)

    ::

      ex. False

- *concat_xfm* (a boolean, optional, default value is False)
    Write joint transformation of two input matrices
    (exclusive with fix_scale_skew, invert_xfm).
    in_trasfo_2 parameter should be used for the second matrix.


    ::

      ex. False

- *fix_scale_skew* (a boolean, optional, default value is False)
    Use secondary matrix to fix scale and skew.
    (exclusive with concat_xfm, invert_xfm).
    in_trasfo_2 parameter should be used for the second matrix.


    ::

      ex. False

**Optional inputs :**

- *in_transfo_2* (a string representing an existing file, optional)
    Second input transformation matrix (valid extensions: [.mat]).

    ::

      ex. '/home/username/data/raw_data/struc2mni.mat'


**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    Final transformation matrix

    ::

      ex. '/home/username/data/derived_data/epi2struct_concate_struc2mni.mat

-------------

Usefull links:

`FSL convert_xfm <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT/UserGuide#convert_xfm>`_

`FSL ConvertXFM - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.fsl.utils.html#convertxfm>`_
