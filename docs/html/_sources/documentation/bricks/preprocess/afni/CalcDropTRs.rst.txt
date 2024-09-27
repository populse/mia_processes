:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==================
CalcDropTRs brick
==================

Drop volumes of bold datasets, using AFNI 3dcalc
------------------------------------------------

If "stop_idx" parameter is set to None or set to default (-1), "stop_idx" will be automatically set to the length of input file.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input 3D file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

**Optional inputs with default value parameters:**

- *output_type* (NIFTI or NIFTI_GZ, optional, , default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *out_prefix* (a string, optional, default value is 'cropped')
    Prefix of the output image.

    ::

      ex. 'cropped_'

- *start_idx* (an integer, optional, default value is 0)
    Start index for in_file_a.

    ::

      ex. 0

- *stop_idx* (an integer, optional, default value is -1)
    Stop index for in_file_a.
    If "stop_idx" parameter is set to None or set to default (-1), "stop_idx" will be automatically set to the length of input file.
    Cannot be lower than or equal to "start_idx" parameters.

    ::

      ex. 10

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out image (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/cropped_func.nii'

-------------

Useful links:

`AFNI 3dcalc <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html>`_

`AFNI Calc - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.utils.html#calc>`_
