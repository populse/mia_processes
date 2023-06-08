:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===================
EstimateModel brick
===================

Model Estimation using classical (ReML - Restricted Maximum Likelihood) algorithms (SPM12 fMRI model estimation).

The beta weights for each condition will be estimated.

Bayesian methods are not implemented in this brick. If you need to use Bayesian method, please open a `ticket <https://github.com/populse/mia_processes/issues>`_

--------------------------------------

**Mandatory inputs parameters:**

- *spm_mat_file* <=> spmmat (a string representing an existing file):
    The SPM.mat file that contains the design specification.

    ::
        ex. '/home/username/data/derived_data/SPM.mat'

**Optional inputs parameters with a default value:**

- *estimation_method* <=> method (a dictionary, optional, default value is {‘Classical’: 1}):
    Estimation procedures for fMRI models:
        - Classical: model parameters are estimated using Restricted Maximum Likelihood (ReML). This correlation can be specified using either an AR(1) or an Independent and Identically Distributed (IID) error model (options chosen at the model specification stage)

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

**Outputs parameters:**

- *beta_images* (a list of items which are a pathlike object or string representing a file):
    Images of estimated regression coefficients "beta_000k" where k indexes the kth regression coefficient .

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
    Be carefull the input smp_mat_file will be overwritten.

    ::
        ex. '/home/username/data/derived_data/SPM.mat'

- *residual_image* (a pathlike object or string representing a file):
    The image of the variance of the error.

    ::
        ex. '/home/username/data/derived_data/ResMS.nii'

- *residual_images* (a list of items which are a pathlike object or string representing a file, optional):
    The individual error "Res_000k" images where k indexes the kth dynamic (fourth dimensional points of the fuctional).
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


-------------

.. [#label] Syntax: mia_processes/nipype EstimateModel brick <=> SPM12 fMRI model estimation.

	    Usefull links:
	    `SPM12 fMRI model estimation <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=73>`_,
	    `nipype EstimateModel <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.model.html#estimatemodel>`_
