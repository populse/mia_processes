:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=======================
EstimateConstrast brick
=======================

Estimate contrasts of interest (SPM12 fMRI contrast manager)
------------------------------------------------------------

To define a new T-contrast:
    - add the contrast name in 'T_contrast_names' parameter
    - add the name of the conditions used for this contrast in 'T_condition_names' parameter
    - add the contrast weight for this contrast in 'T_contrast_weights' parameters
    - if some session should not be include, specify this in 'session' parameter

To define a new F-contrast:
    - define the T-contrasts that should be used for estimate this F-contrast
    - add the contrast name in 'F_contrast_names' parameter
    - add the names of T-contrasts that should be used in 'F_contrast_T_names' parameter

--------------------------------------

**Mandatory inputs parameters:**

- *spm_mat_file* (a string representing an existing file):
    The SPM.mat file that contains the estimate model.

    ::

        ex. '/home/username/data/derived_data/SPM.mat'

- *out_dir_name*
    Out directory name (a string, default is "stats").
    It will be used only if the SPM.mat file is not already in "project/data/derived_data/subjectName_data"
    (for example, an SPM.mat file added in the database with the import or download button).
    This directory will be created in a the folder "project/data/derived_data/subjectName_data".

    :

      ex. "stats"

- *beta_images* (a list of items which are a pathlike object or string representing an existing file)
    Estimated regression coefficients.

    ::

        ex. ['/home/username/data/derived_data/beta_0001.nii',
             '/home/username/data/derived_data/beta_0002.nii',
             '/home/username/data/derived_data/beta_0003.nii',
             '/home/username/data/derived_data/beta_0004.nii',
             '/home/username/data/derived_data/beta_0005.nii',
             '/home/username/data/derived_data/beta_0006.nii',
             '/home/username/data/derived_data/beta_0007.nii',
             '/home/username/data/derived_data/beta_0008.nii']

- *residual_image* (a string representing an existing file):
    Mean-squared image of the residuals.

    ::

        ex. '/home/username/data/derived_data/ResMS.nii'

**Optional inputs parameters with a default value:**

- *T_contrast_names* (a list of string, optional, default value is ['+']):
    Names of each T-contrast.

    ::

        ex. ['Task vs control', 'Task vs baseline']

- *T_condition_names*  (a list of list of string, optional, default value is [['R1_1']]):
    Names of conditions used for each T-contrast.

    ::

        ex. [['T', 'C'], ['T', 'B']]

- *T_contrast_weights* (a list of list of float, optional, default value is [[1.0]]):
    List of contrast weights for each T-contrast.

    ::

        ex. [[1, -1], [1 -1]]

- *use_derivs* (a boolean, optional, default value is False):
    Use derivatives for estimation. Mutually exclusive group_contrast parameter.

    ::

        ex. True

- *group_contrast* (a boolean, optional, default value is False):
    Higher level contrast. Mutually exclusive use_derivs parameter.

    ::

        ex. True

**Optional inputs:**

- *session* (a list of list of float, optional):
    Session list. Length should be equal to the number of sessions, with 1. for sessions to include and 0. elsewhere
    Default is Undefined.

    ::

        ex. ['0', '1.0']

- *F_contrast_name* (a list of string, optional):
    Names of each F contrast.

    ::

        ex. ['Positif effect of condition']

- *F_contrast_T_names* (a list of list of string, optional):
    Names of T contrast used for each F contrast.

    ::

        ex. [['Task vs control', 'Task vs baseline']]

**Outputs parameters:**

- *out_spm_mat_file* (a pathlike object or string representing a file):
    The SPM.mat file containing specification of the design, estimated model parameters and contrast.
    Note that the input smp_mat_file will be overwritten.

    ::

        ex. '/home/username/data/derived_data/SPM.mat'

- *con_images* (a list of items which are a pathlike object or string representing a file, optional):
    Contrast images from a T-contrast.

    ::

        ex. ['/home/username/data/derived_data/con_0001.nii',
             '/home/username/data/derived_data/con_0002.nii',
             '/home/username/data/derived_data/con_0003.nii',
             '/home/username/data/derived_data/con_0004.nii']

- *spmT_images* (a list of items which are a pathlike object or string representing a file, optional):
    Stat images from a T-contrast.

    ::

        ex. ['/home/username/data/derived_data/spmT_0001.nii',
             '/home/username/data/derived_data/spmT_0002.nii',
             '/home/username/data/derived_data/spmT_0003.nii',
             '/home/username/data/derived_data/spmT_0004.nii']

- *ess_images* (a list of items which are a pathlike object or string representing a file, optional):
    Contrast images from a F-contrast.

    ::

        ex. ['/home/username/data/derived_data/ess_0001.nii',
             '/home/username/data/derived_data/ess_0002.nii',
             '/home/username/data/derived_data/ess_0003.nii',
             '/home/username/data/derived_data/ess_0004.nii']

- *spmF_images* (a list of items which are a pathlike object or string representing a file, optional):
    Stat images from a F-contrast.

    ::

        ex. ['/home/username/data/derived_data/spmF_0001.nii',
             '/home/username/data/derived_data/spmF_0002.nii',
             '/home/username/data/derived_data/spmF_0003.nii',
             '/home/username/data/derived_data/spmF_0004.nii']


-------------

Useful links:

`nipype EstimateContrast <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.model.html#estimatemodel>`_
