:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Calc brick
============

Voxel-by-voxel arithmetic on 3D datasets (using AFNI 3dcalc command). 
This bricks assumes that the voxel-by-voxel computations are being performed on datasets that occupy the same space and have the same orientations.    

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import Calc

>>> Calc.help()

**Mandatory inputs parameters:**

- *in_file_a* (a string representing an existing file)
    First input 3D file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/img_1.nii'

**Optional inputs with default value parameters:**

- *expr* (a string, optional)
    Arithmetic expression to apply between a, b and c. 

    Arithmetic expressions are allowed, using + - * / ** ^ and parentheses.
    It is also possible to use built in functions such as sin, cos, sqrt, mean, stdev...
    Some functions are also designed to help implement logical functions such as step(x)  = { 1 if x > 0; 0 if x <= 0 }. 

    For more example, refers to the AFNI documentation.

    ::

      default value. 'a*step(b)'

- *output_type* (NIFTI or NIFTI_GZ, optional)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz

    ::

      default value. NIFTI

- *out_prefix* (a string, optional)
    Prefix of the output image. Default is 'c_'.
    
    ::

        default value. 'c_'

**Optional inputs parameters:** 

- *in_file_b* (a string representing an existing file, optional)
    Second input 3D file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/img_2.nii'

- *in_file_c* (a string representing an existing file, optional)
    Third input 3D file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/img_3.nii'

- *single_idx* (an integer, optional)
    Volume index for in_file_a. The computation will be performed only on this volume.
    Default is Undefined (ie parameter not used)
    ::

      ex. '5'

- *start_idx* (an integer, optional)
    Start index for in_file_a. Requires stop_idx parameter.
    The computation will be performed only only on the volumes between start_idx and stop_idx.
    Default is Undefined (ie parameter not used).
    ::

      ex. '5'

- *stop_idx* (an integer, optional)
    Stop index for in_file_a. Requires start_idx parameter.
    The computation will be performed only only on the volumes between start_idx and stop_idx.
    Default is Undefined (ie parameter not used).
    ::

      ex. '35'

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/c_img_1.nii'

-------------

Usefull links:

`AFNI 3dcalc <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html>`_
`AFNI Calc - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.utils.html#calc>`_
