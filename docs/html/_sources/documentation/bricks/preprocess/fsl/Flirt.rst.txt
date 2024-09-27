:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===========
Flirt brick
===========

Linear (affine) intra- and inter-modal brain image registration using FLIRT (FSL)
---------------------------------------------------------------------------------

FLIRT (FMRIB's Linear Image Registration Tool) is a tool for linear (affine) intra and inter-modal brain image registration.

This brick can be used to calculate the affine transformation that registers the input file to the reference file and to obtain
an output volume where the transform is applied to the input volume to align it with the reference volume (use get_registered_file option).

It is also possible to apply a saved transformation to a volume (apply_xfm or apply_isoxfm options).
The reference volume must still be specified as this sets the voxel and image dimensions of the resulting volume.

Note that the reference image determines the Field of View (FOV) and voxel size of the output registered image.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input file (valid extensions: [.nii, .nii.gz]) to register to the reference file.

    ::

      ex. '/home/username/data/raw_data/b0mean.nii'

- *in_reference* (a string representing an existing file)
    Reference volume (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

**Optional inputs:**

- *get_registered_file* (a boolean, optional, default value is True)
    Obtain an output volume where the transform is applied to the input volume to align it with the reference volume.

    ::

      ex. True

- *apply_xfm* (a boolean, optional, default value is False)
    Apply transformation supplied by in_matrix_file
    Mutually exclusive with apply_isoxfm.

    ::

      ex. False

- *apply_isoxfm* (a float, optional)
    Apply transformation supplied by in_matrix_file but forces isotropic resampling (using the float provided)
    Mutually exclusive with apply_xfm.

    ::

      ex. 1.0

- *in_matrix_file* (a string representing an existing file, optional)
    Input 4x4 affine matrix transformation. Matrix used if apply_isoxfm or apply_xfm are used.

    ::

      ex. '/home/username/data/derived_data/matrix_struct_to_diff.mat'

- *angle_rep* (euler or quaternion, optional, default value is euler)
    Representation of rotation angles.

    ::

      ex. 'euler'

- *bbrslope* (a float, optional)
    Value of BBR slope, if BBR used for cost function.

    ::

      ex. 0.1

- *bbrtype* (signed, global_abs, local_abls, optional, default value is signed)
    Type of BBR cost function, if BBR used for cost function.

    ::

      ex. 'signed'

- *bgvalue* (a float, optional)
    Use specified background value for points outside FOV

    ::

      ex. 5.2

- *bins* (an integer, optional, default value is 256)
    Number of histogram bins

    ::

      ex. 256

- *coarse_search* (an integer, optional, default value is 60)
    Coarse search delta angle

    ::

      ex. 60

- *cost* (mutualinfo, corratio, normcor, normi, leastsq, labeldiff or bbr, optional, default value is corratio)
    | Cost function. There are two main types:
    | - intra-model (same modality): least squares and normalised correlation
    | - inter-model (different modalities): correlation ratio and mutual information-based

    For image from the same modality, both intra and inter-model type can be used.
    The recommended options (to try first) are correlation ratio for inter-modal and normalised correlation for intra-modal.

    ::

      ex. 'corratio'

- *cost_func* (mutualinfo, corratio, normcor, normi, leastsq, labeldiff or bbr, optional, default value is corratio)
    Cost function (maps to argument searchcost)

    ::

      ex. 'corratio'

- *data_type* (char, short, int, float or double, optional)
    Force output data type

    ::

      ex. 'char'

- *dof* (6,7, 9 or 12, optional, default is 12)
    Number of transform degrees of freedom.
    The number of degrees of freedom depends on whether the images to be registered are intra/inter-subject and small/large FOV.
    | If the images are from the same subject:
    | - large FOV --> 6 dof is appropriate (7 if the scanner voxel size may have changed)
    | - small FOV/single slice --> 3 dof is appropriate  (4 if the scanner voxel size may have changed)
    | If the images are from the different subjects:
    | - large FOV --> 12 dof is appropriate
    | - small FOV --> 6 dof is appropriate

    For 2D-to-2D registration, used rigid2D option.

    ::

      ex. 6

- *echospacing* (a float, optional)
    Value of EPI echo spacing - units of seconds.

    ::

      ex. 2.0

- *fine_search* (an integer, optional, default value is 18)
    Fine search delta angle

    ::

      ex. 18

- *fieldmap* (a string representing an existing file, optional)
    Fieldmap image in rads/s - must be already registered to the reference image

    ::

      ex. '/home/username/data/raw_data/field_map_registered.nii'

- *fieldmapmask* (a string representing an existing file, optional)
    Mask for fieldmap image

    ::

      ex. '/home/username/data/raw_data/field_map_registered_mask.nii'

- *force_scaling* (a boolean, optional, default value is False)
    Force rescaling even for low-res images.

    ::

      ex. False

- *interp* (trilinear, nearestneighbour, sinc,  or spline, optional, default value is trilinear),
    Final interpolation method used in reslicing.
    The interpolation is only used for the final transformation (and in applyxfm), not in the registration calculations.

    ::

      ex. 'trilinear'

- *in_weight* (a string representing an existing file, optional)
    File for input weighting volume.
    Cost function weighting volumes can be specified (either for input volume or for reference volume, see ref_weight option).
    This allows the cost function to have a different weighting at each voxel, which is useful for excluding areas (weight=0)
    of no interest, or increasing the weighting around important structures such as the ventricles.
    Note that this is different from masking the original images, as masking introduces artificial boundaries whereas weighting does not.


    ::

      ex. '/home/username/data/raw_data/in_weight.nii'

- *min_sampling* (a float, optional)
    Set minimum voxel dimension for sampling

    ::

      ex. 0.5

- *no_clamp* (a boolean, optional, default value is False)
    Do not use intensity clampinp.

    ::

      ex. False

- *no_resample* (a boolean, optional, default value is False)
    Do not change input sampling.

    ::

      ex. False

- *no_resample_blur* (a boolean, optional, default value is False)
    Do not use blurring on downsampling.

    ::

      ex. False

- *no_search* (a boolean, optional, default value is False)
    Set all angular searches to ranges 0 to 0.

    ::

      ex. False


- *output_type* ('NIFTI' or 'NIFTI_GZ', optional, default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   - NIFTI: \*.nii
    |   - NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *padding_size* (an integer, optional)
    For applyxfm: interpolates outside image by size

    ::

      ex. 2

- *pedir* (an integer, optional)
    Phase encode direction of EPI  with 1/2/3=x/y/z and -1/-2/-3=-x/-y/-z

    ::

      ex. 2

- *rigid2D* (a boolean, optional, default value is False)
    2D-to-2D registration mode. Use 2D rigid body mode and ignore dof

    ::

      ex. False

- *ref_weight* (a string representing an existing file, optional)
    File for reference weighting volume.
    Cost function weighting volumes can be specified (either for input volume or for reference volume, see also in_weight option).
    This allows the cost function to have a different weighting at each voxel, which is useful for excluding areas (weight=0)
    of no interest, or increasing the weighting around important structures such as the ventricles.
    Note that this is different from masking the original images, as masking introduces artificial boundaries whereas weighting does not.


    ::

      ex. '/home/username/data/raw_data/in_weight.nii'

- *save_log* (a boolean, optional, default value is False)
    Save fsl log

    ::

      ex. False

- *schedule* (a string representing an existing file, optional)
    Replaces default schedule. The schedule file specifies what transformations/DOF are allowed and how the optimisation is performed.


    ::

      ex. '{FSLDIR}/etc/flirtsch/sch2D_6dof'

- *searchr_x* (a list of 2 integer, optional, default value is [-90, 90])
    Search angles along x-axis, in degrees [min angle, max angle]

    ::

      ex. [-90, 90]

- *searchr_y* (a list of 2 integer, optional, default value is [-90, 90])
    Search angles along y-axis, in degrees [min angle, max angle]

    ::

      ex. [-90, 90]

- *searchr_z* (a list of 2 integer, optional, default value is [-90, 90])
    Search angles along z-axis, in degrees [min angle, max angle]

    ::

      ex. [-90, 90]

- *sinc_width* (an integer, optional)
    Full-width (window width) in voxels if sinc is used for final interpolation

    ::

      ex. 7

- *sinc_window* (rectangular, hanning or blackman, optional)
    Type of windowing function if sinc is used for final interpolation

    ::

      ex. 7

- *use_qform* (a boolean, optional, default value is False)
    Initialize using sform or qform

    ::

      ex. False

- *wm_seg* (a string representing an existing file, optional)
    White matter segmentation volume needed by BBR cost function


    ::

      ex. '/home/username/data/derived_data/wm_seg.nii'

- *wmcoords* (a string representing an existing file, optional)
    White matter boundary coordinates for BBR cost function


    ::

      ex. '/home/username/data/derived_data/wm_coors.nii'

- *wmnorms* (a string representing an existing file, optional)
    "White matter boundary normals for BBR cost function


    ::

      ex. '/home/username/data/derived_data/wm_norms.nii'

**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file, optional)
    Registered file

    ::

      ex. '/home/username/data/derived_data/b0mean_registered_with_T1w.nii'

- *out_log* (a pathlike object or string representing a file, optional)
    Log file

    ::

      ex. '/home/username/data/derived_data/b0mean_flirt_log.txt'

- *out_matrix* (a pathlike object or string representing a file, optional)
    Calculated affine transform

    ::

      ex. '/home/username/data/derived_data/b0mean_T1w_flirt.mat'

-------------

Useful links:

`FSL FLIRT <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT>`_

`FSL FLIRT FAQ <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT/FAQ>`_

`FSL FLIRT - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.fsl.preprocess.html#flirt>`_
