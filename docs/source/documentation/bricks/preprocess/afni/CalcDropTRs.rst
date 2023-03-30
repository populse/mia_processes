:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
CalcDropTRs brick
============

Drop volumes of bold datasets (using AFNI 3dcalc command). 

If "stop_idx" parameter is set to None or set to default (-1), "stop_idx" will be automatically set to the length of input file. 

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import CalcDropTRs

>>> CalcDropTRs.help()

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input 3D file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

**Optional inputs with default value parameters:**

- *output_type* (NIFTI or NIFTI_GZ, optional)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz

    ::

      default value. NIFTI

- *out_prefix* (a string, optional)
    Prefix of the output image. Default is 'cropped_'.
    
    ::

      default value. 'cropped_'

- *start_idx* (an integer, optional)
    Start index for in_file_a.

    ::

      default value. 0

- *stop_idx* (an integer, optional)
    Stop index for in_file_a.
    If "stop_idx" parameter is set to None or set to default (-1), "stop_idx" will be automatically set to the length of input file. 
    Cannot be lower than or equal to "start_idx" parameters. 

    ::

      default value. -1

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/cropped_func.nii'

-------------

Usefull links:

`AFNI 3dcalc <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html>`_
`AFNI Calc - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.utils.html#calc>`_
