:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Volreg brick
============

Register an input volume to a base volume using AFNI 3dvolreg command

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import Volreg

>>> Volreg.help()

**Inputs parameters:**

- *in_file* (a string representing an existing file)
    Input file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

- *interpolation* (‘Fourier’ or ‘linear’ or ‘cubic’ or ‘quintic’ or ‘heptic’, optional)
    Different interpolation methods:
      -Fourier = Use a Fourier method (the default: most accurate; slowest).
      -linear  = Use linear (1st order polynomial) interpolation (least accurate).
      -cubic   = Use the cubic (3rd order) Lagrange polynomial interpolation.
      -quintic = Use the quintic (5th order) Lagrange polynomial interpolation.
      -heptic  = Use the heptic (7th order) Lagrange polynomial interpolation.
    Default is 'heptic'.

    ::

      ex. heptic

- *twopass* (a boolean, optional)
    Do two passes of the registration algorithm:
    (1) with smoothed base and data bricks, with linear interpolation, to get a crude alignment, then'
    (2) with the input base and data bricks, to get a fine alignment. 
    This method is useful when aligning high-resolution datasets that may need to be moved more than a few voxels to be aligned.
    Default is False. 

    ::

      ex. False

- *zpad* (a integer, optional)
    Zeropad around the edges by ‘n’ voxels during rotations. 
    Default is n = 4.
    ::

      ex. 4

- *copyorigin* (a boolean, optional)
    Copy base file origin coords to output.
    Default is False. 

    ::

      ex. False

- *timeshift* (a boolean, optional)
    Time shift to mean slice time offset.
    Default is False. 

    ::

      ex. False

- *in_weight_volume* (a tuple (a string representing an existing file, Integer), optional)
    Weights for each voxel specified by a file with an optional volume number (defaults to 0).
    Default is Undefined (ie parameter not used).

    ::

      ex. False

- *save_oned_matrix* (a boolean, optional)
    Save the transformation matrix  oned matrix. 

    ::

      ex. False

- *save_md1d_file* (a boolean, optional)
    Save max displacement outputfile (md1d) file. 

    ::

      ex. False

- *output_type* (NIFTI or NIFTI_GZ, optional)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *out_prefix* (a string, optional)
    Prefix of the output image. Default is 'reg_'.
    
    ::

      ex. 'reg_'


**Outputs parameters:**

- *out_file* (a strings representing a file)
    Register file (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/reg_func.nii'

- *oned_file* (a strings representing a file)
   The movement parameters file (extensions: [.txt]).
    
    ::

      ex. '/home/username/data/derived_data/reg_func_oned.txt'

- *oned_matrix* (a strings representing a file, optional)
   The transformation matrix (extensions: [.aff12.1D]).
    
    ::

      ex. '/home/username/data/derived_data/reg_func_oned_matrix.aff12.1D'

- *md1d_file* (a strings representing a file, optional)
   The transformation matrix (extensions: [.aff12.1D]).
    
    ::

      ex. '/home/username/data/derived_data/reg_func_md.1D'

-------------

Usefull links:

`AFNI 3dvolreg <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dvolreg.html>`_
`AFNI Volreg - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.preprocess.html#volreg>`_
