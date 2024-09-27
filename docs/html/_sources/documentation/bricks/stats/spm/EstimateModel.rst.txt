:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===================
EstimateModel brick
===================

Model Estimation using classical (ReML - Restricted Maximum Likelihood) procedure (SPM12 fMRI model estimation)
---------------------------------------------------------------------------------------------------------------

This step will estimate on fitted data the parameters of the model created during the model specification stage,
i.e. estimate the value of the β associated with each conditions.

Bayesian methods are not currently implemented in this brick. If you need to use Bayesian procedures, please open a `ticket <https://github.com/populse/mia_processes/issues>`_.

.. warning::
    In some cases, this brick is not compatible with the Mia V2 controller (risk of Mia crash).
    It is strongly recommended to use the Mia V1 controller with this brick (see: File > Mia Preferences, to change the controller version).

--------------------------------------

**Mandatory inputs parameters:**

- *spm_mat_file* <=> spmmat (a string representing an existing file):
    The SPM.mat file that contains the design specification.

    ::

        ex. '/home/username/data/derived_data/SPM.mat'

**Optional inputs parameters with a default value:**

- *out_dir_name*
    Out directory name (a string, default is "stats").
    It will be used only if the SPM.mat file is not already in "project/data/derived_data/subjectName_data"
    (for example, an SPM.mat file added in the database with the import or download button).
    This directory will be created in a the folder "project/data/derived_data/subjectName_data".

    :

      ex. "stats"

- *estimation_method* <=> method (a dictionary, optional, default value is {‘Classical’: 1}):
    Estimation procedures for fMRI models:
        - Classical: model parameters are estimated using Restricted Maximum Likelihood (ReML). This correlation can be specified using either an AR(1) or an Independent and Identically Distributed (IID) error model (options chosen at the `model specification stage <bricks/stats/spm/Level1Design.html>`_)

        - Bayesian methods are not implemented yet.

    ::

        ex. {'Classical': 1}

- *write_residuals* <=> write_residuals (a boolean, optional, default value is False):
    Write images of residuals for each dynamique.

    ::

        ex. True

- *version* (a string, optional, default value is spm12):
    Version of spm.

    ::

        ex. 'spm12'

**Optional inputs:**

- *tot_reg_num* (an integer, optional):
    The total number of estimated regression coefficients.
    Default is Undefined. In this case, the number of estimated regression coefficients used will be the one defined in the database for the smp_mat_file.

    ::

        ex. 8

- *factor_info* (a list of items which are a dictionary with keys which are ‘name’ and ‘levels’ and with values which are a string or an integer, optional):
    factor_info  parameter used in `Level1Design brick <bricks/stats/spm/Level1Design.html>`_ to specified a factorial design. In this case, SPM will automatically generate the contrasts necessary to test for the main effects and interactions.
    This parameter should be used only if factorial design is specified in Level1Design brick.

    ::

        ex. [{"name": "Factor1", "levels": 2}, {"name": "Factor2", "levels": 2}]

- *bases* (a dictionary with keys which are 'hrf' or 'fourier' or 'fourier_han' or 'gamma' or 'fir' and with values which are a dictionary with keys which are 'derivs' or 'length' or 'order' and with values which are a list or a float or an integer, optional):
    Information used in `Level1Design brick <bricks/stats/spm/Level1Design.html>`_ brick to define basic functions for modeling hemodynamic response.
    Default is Undefined.
    This parameter should be used only if factorial design is specified in Level1Design brick.

    ::

        ex. {"hrf": {"derivs": [1, 1]}}

**Outputs parameters:**

- *beta_images* (a list of items which are a pathlike object or string representing a file):
    Images of estimated regression coefficients "beta_000k" where k indexes the kth regression coefficient.

    ::

        ex. ['/home/username/data/derived_data/beta_0001.nii',
             '/home/username/data/derived_data/beta_0002.nii',
             '/home/username/data/derived_data/beta_0003.nii',
             '/home/username/data/derived_data/beta_0004.nii',
             '/home/username/data/derived_data/beta_0005.nii',
             '/home/username/data/derived_data/beta_0006.nii',
             '/home/username/data/derived_data/beta_0007.nii',
             '/home/username/data/derived_data/beta_0008.nii']

- *mask_image* (a pathlike object or string representing a file):
    The mask image indicating which voxels were included in the analysis.

    ::

        ex. '/home/username/data/derived_data/mask.nii'

- *out_spm_mat_file* (a pathlike object or string representing a file):
    The SPM.mat file containing specification of the design and estimated model parameters.
    Note that the input smp_mat_file will be overwritten.

    ::

        ex. '/home/username/data/derived_data/SPM.mat'

- *residual_image* (a pathlike object or string representing a file):
    The image of the variance of the error.

    ::

        ex. '/home/username/data/derived_data/ResMS.nii'

- *residual_images* (a list of items which are a pathlike object or string representing a file, optional):
    The individual error "Res_000k" images where k indexes the kth dynamic (fourth dimensional points of the functional).
    These images are generated only if write_residuals is True.

    ::

        ex. ['/home/username/data/derived_data/Res_0001.nii',
             '/home/username/data/derived_data/Res_0002.nii',
             '/home/username/data/derived_data/Res_0003.nii',
            ...,
             '/home/username/data/derived_data/Res_0238.nii',
             '/home/username/data/derived_data/Res_0239.nii',
             '/home/username/data/derived_data/Res_0240.nii']

- *RPVimage* (a pathlike object or string representing a file):
    The image of the estimated resolution elements per voxel.

    ::

        ex. '/home/username/data/derived_data/RPV.nii'

- *con_images* (a list of items which are a pathlike object or string representing a file, optional):
    Contrast images from a t-contrast. Only created if factor_info used in Level1Design brick.

    ::

        ex. ['/home/username/data/derived_data/con_0005.nii',
             '/home/username/data/derived_data/con_0006.nii',
             '/home/username/data/derived_data/con_0007.nii',
             '/home/username/data/derived_data/con_0008.nii']

- *spmT_images* (a list of items which are a pathlike object or string representing a file):
    Stat images from a t-contrast. Only created if factor_info used in Level1Design brick.

    ::

        ex. ['/home/username/data/derived_data/spmT_0005.nii',
             '/home/username/data/derived_data/spmT_0006.nii',
             '/home/username/data/derived_data/spmT_0007.nii',
             '/home/username/data/derived_data/spmT_0008.nii']

- *ess_images* (a list of items which are a pathlike object or string representing a file, optional):
    Contrast images from a f-contrast. Only created if factor_info used in Level1Design brick.

    ::

        ex. ['/home/username/data/derived_data/ess_0001.nii',
             '/home/username/data/derived_data/ess_0002.nii',
             '/home/username/data/derived_data/ess_0003.nii',
             '/home/username/data/derived_data/ess_0004.nii']

- *spmF_images* (a list of items which are a pathlike object or string representing a file):
    Stat images from a f-contrast. Only created if factor_info used in Level1Design brick.

    ::

        ex. ['/home/username/data/derived_data/spmF_0001.nii',
             '/home/username/data/derived_data/spmF_0002.nii',
             '/home/username/data/derived_data/spmF_0003.nii',
             '/home/username/data/derived_data/spmF_0004.nii']



-------------

Useful links:

`SPM12 fMRI model estimation <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=73>`_,
`nipype EstimateModel <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.model.html#estimatemodel>`_
