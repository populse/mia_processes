:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Volreg brick
============

Registers each 3D volume from the input dataset to the base volume using AFNI 3dvolreg command

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

**Optional inputs with default value parameters:**

- *copyorigin* (a boolean, optional, default value is False)
    Copy base file origin coords to output.

    ::

      ex. False

- *interpolation* (Fourier or linear or cubic or quintic or heptic, optional, default value is heptic)
    Different interpolation methods:
      - Fourier = Use a Fourier method (the default: most accurate; slowest).
      - linear  = Use linear (1st order polynomial) interpolation (least accurate).
      - cubic   = Use the cubic (3rd order) Lagrange polynomial interpolation.
      - quintic = Use the quintic (5th order) Lagrange polynomial interpolation.
      - heptic  = Use the heptic (7th order) Lagrange polynomial interpolation.

    ::

      ex. heptic

- *output_type* (NIFTI or NIFTI_GZ, optional, default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   - NIFTI: \*.nii
    |   - NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *out_prefix* (a string, optional, default value is 'reg')
    Prefix of the output image.

    ::

      ex. 'reg_'

- *save_oned_matrix* (a boolean, optional, default value is False)
    Save the transformation matrix  oned matrix.

    ::

      ex. False

- *save_md1d_file* (a boolean, optional, default value is False)
    Save max displacement outputfile (md1d) file.

    ::

      ex. False

- *timeshift* (a boolean, optional, default value is False)
    Time shift to mean slice time offset.

    ::

      ex. False

- *twopass* (a boolean, optional, default value is False)
    Do two passes of the registration algorithm:
    (1) with smoothed base and data bricks, with linear interpolation, to get a crude alignment, then'
    (2) with the input base and data bricks, to get a fine alignment.
    This method is useful when aligning high-resolution datasets that may need to be moved more than a few voxels to be aligned.

    ::

      ex. False

- *zpad* (a integer, optional, default value is 4)
    Zeropad around the edges by ‘n’ voxels during rotations.

    ::

      ex. 4

**Optional inputs with default value parameters:**

- *in_weight_volume* (a tuple (a string representing an existing file, Integer), optional)
    Weights for each voxel specified by a file with an optional volume number (defaults to 0).
    Default is Undefined (ie parameter not used).

    ::

      ex. ('/home/username/data/raw_data/mask.nii', 0)


**Outputs parameters:**

- *md1d_file* (a strings representing a file, optional)
    The transformation matrix (extensions: [.aff12.1D]).

    ::

      ex. '/home/username/data/derived_data/reg_func_md.1D'

- *oned_file* (a strings representing a file)
    Movement parameters file (extensions: [.txt]).

    ::

      ex. '/home/username/data/derived_data/reg_func_oned.txt'

- *oned_matrix* (a strings representing a file, optional)
    Transformation matrix (extensions: [.aff12.1D]).

    ::

      ex. '/home/username/data/derived_data/reg_func_oned_matrix.aff12.1D'

- *out_file* (a strings representing a file)
    Register file (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/reg_func.nii'

-------------

Usefull links:

`AFNI 3dvolreg <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dvolreg.html>`_
`AFNI Volreg - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.preprocess.html#volreg>`_
