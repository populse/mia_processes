:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Binarize brick
============
 
 Binarize an input volume (or volume-encoded surface file). 

 Binarization can be done based on threshold (using min or rmin and/or max or rmax options) or on matched values (using match option). 
 All voxels in the threshold/match are set to the binval value (1 by default) and all voxels out of range are set to the binvalnot value (0 by default)

 By default all the option are set to Undefined.

>>> from mia_processes.bricks.preprocess.freesurfer import Binarize
>>> Binarize.help()

**Inputs parameters:**

- *in_file* (a string representing an existing file)
    Input volume to be binarized. 

    ::

      ex. /home/username/data/raw_data/T1w.nii


- *min* (float)
    Minimum voxel threshold.
    Cannot be used with match option.

    ::

      ex. 0.5

- *max* (float)
    Maximum voxel threshold
    Cannot be used with match option.

    ::

      ex. 100.0

- *rmin* (float)
    Compute min based on rmin*globalmean.
    Cannot be used with match option.

    ::

      ex. 0.5

- *rmax* (float)
    Compute max based on rmax*globalmean.
    Cannot be used with match option.

    ::

      ex. 100.0

- *match* (a list of items which are an integer)
    This option allows to use match instead of threshold to binarize.  Any number of match values can be specified. 
    Cannot be used with min/rmin and max/rmax option. 

    For example, it can be used to extract all the voxel with a label equal to 4 in the input file. 

    ::

      ex. [4]

- *get_count_file* (a boolean)
    Save number of hits in ascii file (hits, ntotvox, pct). 
    Four numbers are saved: the number of voxels that match(nhits), the volume of the voxels that match, the total number of voxels in the volume (nvoxtot), and the percent matching (100*nhits/nvoxtot). 
    Default is False

    ::

      ex. False

- *binval* (an integer)
    Value to use for all voxels that are in the threshold/match. 
    By default this value is 1.

    ::

      ex. 1

- *binvalnot* (an integer)
    Value to use for all voxels that are out of range (not in the threshold/match.)
    By default this value is 0.

    ::

      ex. 0

- *invert* (a boolean)
    Invert binval and binvalnot (set binval=0 and binvalnot=1). 
    Default is False.

    ::

      ex. False

- *bin_col_nul* (a boolean)
    Set binarized voxel value to its column number. 
    Default is False.

    ::

      ex. False

- *frame_no* (an integer)
    Use 0-based frame of input (default is 0)

    ::

      ex. 0

- *abs*  (a boolean)
    Take absolute value of the input volume (ie, make input unsigned). 
    Default is False
    ::

      ex. False

- *zero_edges*  (a boolean)
    Set the first and last planes in all dimensions to binvalnot(0 by default). This makes sure that all the voxels on the edge of the imaging volume are zero. 
    Default is False
    ::

      ex. False

- *zero_slice_edges*  (a boolean)
    Same as zero_edges but only for slices
    Default is False
    ::

      ex. False

- *dilate*  (an integer)
    Dilate binarization in 3D 

    ::

      ex. 1

- *erode*  (an integer)
    Erode binarization in 3D (after any dilation)
    ::

      ex. 1

- *erode2d*  (an integer)
    Erode binarization in 2D (after any 3D erosion)
    ::

      ex. 1

- *out_suffix* (a string)
    Output suffix. Default is '_thresh'

    ::

      ex. _thresh


- *output_type* ('NIFTI' or 'NIFTI_GZ' or 'MGZ')
    | Format of the output image (one of NIFTI, NIFTI_GZ, MGZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz
    |   MGZ: \*.mgz

    ::

      ex. NIFTI


**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    Brain-extracted image

    ::

      ex. /home/username/data/raw_data/T1w_desc-brain.nii

- *count_file* (a pathlike object or string representing a file)
    Ascii file containing number of hits.

    ::

      ex. /home/username/data/raw_data/T1w_binarize_count_file.txt

-------------

Usefull links:
`Freesurfer Binarize <https://surfer.nmr.mgh.harvard.edu/fswiki/mri_binarize>`_
`Freesurfer Binarize - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.freesurfer.model.html#binarize>`_