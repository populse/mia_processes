:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
SkullStrip brick
============

Extract the brain from surrounding tissue from MRI T1-weighted images (using AFNI 3dSkullStrip command).

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    'A 3D-T1 file to skull-strip (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'


**Optional inputs with default value parameters:**

- *output_type* (NIFTI or NIFTI_GZ, optional, default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *out_prefix* (a string, optional, default value is 'ss')
    Prefix of the output image.

    ::

        ex. 'ss_'

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out image (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/ss_T1w.nii'

-------------

Usefull links:

`AFNI 3dSkullStrip <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dSkullStrip.html>`_
`AFNI SkullStrip - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.preprocess.html#skullstrip>`_
