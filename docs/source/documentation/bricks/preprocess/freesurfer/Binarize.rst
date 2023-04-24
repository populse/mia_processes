:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==============
Binarize brick
==============

Binarize an input volume (or volume-encoded surface file).

Binarization can be done based on threshold (using min or rmin and/or max or rmax options) or on matched values (using match option).
All voxels in the threshold/match are set to the binval value (1 by default) and all voxels out of range are set to the binvalnot value (0 by default)

By default all the option are set to Undefined.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input volume to be binarized.

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

**Optional inputs with default value  parameters:**

- *abs* (a boolean, optional, default value is False)
    Take absolute value of the input volume (ie, make input unsigned).

    ::

      ex. False

- *bin_col_nul* (a boolean, optional, default value is False)
    Set binarized voxel value to its column number.

    ::

      ex. False

- *get_count_file* (a boolean, optional, default value is False)
    Save number of hits in ascii file (hits, ntotvox, pct).
    Four numbers are saved: the number of voxels that match(nhits), the volume of the voxels that match, the total number of voxels in the volume (nvoxtot), and the percent matching (100*nhits/nvoxtot).

    ::

      ex. False

- *invert* (a boolean, optional, default value is False)
    Invert binval and binvalnot (set binval=0 and binvalnot=1).

    ::

      ex. False

- *max* (a float, optional, default value is 100.0)
    Maximum voxel threshold
    Cannot be used with match option.

    ::

      ex. 100.0

- *min* (a float, optional, default value is 0.0)
    Minimum voxel threshold.
    Cannot be used with match option.

    ::

      ex. 0.0

- *out_suffix* (a string, optional, default value is '_thresh')
    Output suffix.

    ::

      ex. _thresh


- *output_type* ('NIFTI' or 'NIFTI_GZ' or 'MGZ', optional, default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ, MGZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz
    |   MGZ: \*.mgz

    ::

      ex. NIFTI

- *zero_edges* (a boolean, optional, default value is False)
    Set the first and last planes in all dimensions to binvalnot(0 by default). This makes sure that all the voxels on the edge of the imaging volume are zero.

    ::

      ex. False

- *zero_slice_edges* (a boolean, optional, default value is False)
    Same as zero_edges but only for slices

    ::

      ex. False

**Optional inputs parameters:**

- *binval* (an integer, optional)
    Value to use for all voxels that are in the threshold/match.
    By default this value is 1.

    ::

      ex. 1

- *binvalnot* (an integer, optional)
    Value to use for all voxels that are out of range (not in the threshold/match.)
    By default this value is 0.

    ::

      ex. 0

- *dilate* (an integer, optional)
    Dilate binarization in 3D

    ::

      ex. 1

- *erode* (an integer, optional)
    Erode binarization in 3D (after any dilation)
    ::

      ex. 1

- *erode2d* (an integer, optional)
    Erode binarization in 2D (after any 3D erosion)
    ::

      ex. 1

- *frame_no* (an integer, optional)
    Use 0-based frame of input (default is 0)

    ::

      ex. 0

- *match* (a list of items which are an integer, optional)
    This option allows to use match instead of threshold to binarize.  Any number of match values can be specified.
    Cannot be used with min/rmin and max/rmax option.

    For example, it can be used to extract all the voxel with a label equal to 4 in the input file.

    ::

      ex. [4]

- *rmax* (a float, optional)
    Compute max based on rmax*globalmean.
    Cannot be used with match option.

    ::

      ex. 100.0

- *rmin* (a float, optional)
    Compute min based on rmin*globalmean.
    Cannot be used with match option.

    ::

      ex. 0.5

**Outputs parameters:**

- *count_file* (a pathlike object or string representing a file)
    Ascii file containing number of hits.

    ::

      ex. '/home/username/data/raw_data/T1w_binarize_count_file.txt'

- *out_file* (a pathlike object or string representing a file)
    Brain-extracted image

    ::

      ex. '/home/username/data/raw_data/T1w_desc-brain.nii'

-------------

Usefull links:
`Freesurfer Binarize <https://surfer.nmr.mgh.harvard.edu/fswiki/mri_binarize>`_
`Freesurfer Binarize - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.freesurfer.model.html#binarize>`_
