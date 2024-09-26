:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===========================
OneSampleTTestDesign brick
===========================

Create SPM design for one sample t-test (SPM12 fMRI factorial Design)
---------------------------------------------------------------------------------------------------------------

This brick can be used for second-level analysis. Level one design should be first performed.

--------------------------------------

**Mandatory inputs parameters:**

- *in_files* <=> scans* [#label]_ (a list of at least 2 items which are a string representing an existing file):
    Input contrasts files.

    ::

        ex. ['/home/username/data/derived_data/sub-001_data/con_001.nii', '/home/username/data/derived_data/sub-002_data/con_001.nii']

**Optional inputs parameters with a default value:**

- *out_dir_name* (a string, optional, default value is "spm_stat_2ndLevel"):
    Name of the directory where the SPM.mat file containing the specified design matrix will be written.
    This directory will be created in the derived_data folder of the project.

    ::

        ex. 'spm_stat_2ndLevel'

- *threshold_masking* <=> Threshold masking* [#label]_ (one of None, Absolute or Relative, optional, default value is None):
    Images are thresholded at a given value and only voxels at which all images exceed the threshold are included:

        | - None: No threshold masking
        | - Absolute: threshold masking done using an absolute threshold value, threshold_mask_value parameter should be filled
        | - Relative: threshold masking done using a proportion of the global value, threshold_mask_value parameter should be filled

    ::

        ex. 'Relative'

- *use_implicit_mask* <=> Implicit Mask* [#label]_ (a boolean, optional, default value is True):
    Use implicit mask NaNs or zeros to threshold.

    ::

        ex. True

- *global_calc* <=> Global Calculation * [#label]_ (one of Omit, Mean, User, optional, default value is Omit):
    This option is for PET or VBM data (not second level fMRI). There are three methods for estimating global effects:

        | - Omit: assuming no other options requiring the global value chosen
        | - User: enter your own vector of global values using the global_calc_values parameter
        | - Mean: SPM standard mean voxel value (within per image fullmean/8 mask)

    ::

        ex. 'Omit'

- *no_grand_mean_scaling* <=> Overall grand mean scalling* [#label]_ (a boolean, optional, default value is True):
    This option is for PET or VBM data (not second level fMRI).
    Do not perform grand mean scaling.

    ::

        ex. True

- *normalisation* <=> Normalisation* [#label]_ (one of 1, 2, 3, optional, default value is 1):
    This option is for PET or VBM data (not second level fMRI). Global normalisation:

        | - 1: None
        | - 2: Proportional
        | - 3: ANCOVA

**Optional inputs:**

- *covariates_names* <=> Covariate.Name* [#label]_ (a list of string, optional):
    Names of the covariates.
    The covariate option allows to specified covariates or nuisance variables.
    The aims is to examine whether there are any correlations between this additionl data and brain activation.

    covariates_vectors, covariates_interactions and covariates_centerings should be also filled.

    ::

        ex. ['Reaction time']

- *covariates_vectors* <=> Covariate.Vector* [#label]_ (a list of list of float, optional):
    Vectors of the covariates values.
    For each covariate, the values should be entered "per subject" (i.e all for subject 1, then all for subject 2,...).

    covariates_names, covariates_interactions and covariates_centerings should be also filled.

    ::

        ex. [[0.0, 1.0, 2.0]]


- *covariates_interactions* <=> Covariate.Interactions* [#label]_ (a list of int among 1, 2, 3, 4 optional):
    For each covariate you have defined, there is an opportunity to create an additional regressor that
    is the interaction between the covariate and a chosen experimental facto:

        | - 1: None
        | - 2: With Factor 1
        | - 3: With Factor 2
        | - 4: With Factor 3

    covariates_vectors, covariates_names and covariates_centerings should be also filled.

    ::

        ex. [1]

- *covariates_centerings* <=> Covariate.Centering* [#label]_ (a list of int among 1, 2, 3, 4, 5, 6, 7, 8, optional):
    Covariates centering values:

        | - 1: Overall mean
        | - 2: Factor 1 mean
        | - 3: Factor 2 mean
        | - 4: Factor 3 mean
        | - 5: No centering
        | - 6: User specified value
        | - 7: As implied by ANCOVA
        | - 8: GM

    covariates_vectors, covariates_interactions and covariates_names should be also filled.

    ::

        ex. [1]

- *threshold_mask_value* <=> Threshold masking.Threshold* [#label]_ (a float, optional):
    Value of the threshold.
    If "Absolute" is used for threshold_masking parameter, enter the absolute value of the threshold.
    If "Relative" is used for threshold_masking parameter, enter the threshold as a proportion of the global value.

    ::

        ex. 0.8

- *explicit_mask_file* <=> Explicit mask* [#label]_ (a string that representing an exiting file, optional):
    Use an implicit mask file to threshold.
    Explicit masks are other images containing (implicit) masks that are to be applied to the current analysis.
    All voxels with value NaN (for image data-types with a representation of NaN), or zero (for other data types) are excluded from the analysis.

    Explicit mask images can have any orientation and voxel/image size. Nearest neighbour interpolation of a mask image is used if the voxel centers of the input images do not coincide with that of the mask image.

    ::

        ex. '/home/username/data/derived_data/mask.nii'

- *global_calc_values* <=> Global Calculation.Global values* [#label]_ (a list of float, optional)
    Vector of global values used for global calculation (to be filled only if "User" choose for global_calc parameter )

    ::

        ex. [0.0, 1.1, 0.2]

**Outputs parameters:**

- *out_spm_mat_file* (a pathlike object or string representing a file):
    The SPM.mat file containing specification of the design and estimated model parameters.
    Note that the input smp_mat_file will be overwritten.

    ::

        ex. '/home/username/data/derived_data/spm_stat_2ndLevel/SPM.mat'

-------------

.. [#label] Syntax: mia_processes/nipype OneSampleTTestDesign brick <=> SPM12 fMRI factorial design.

Useful links:

`SPM12 fMRI Factorial Design <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=81>`_,
`nipype OneSampleTTestDesign <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.model.html#onesamplettestdesign>`_
