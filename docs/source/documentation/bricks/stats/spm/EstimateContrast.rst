:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=======================
EstimateConstrast brick
=======================


Estimate contrasts of interest (SPM12 fMRI contrast manager).

--------------------------------------

**Mandatory inputs parameters:**

- *spm_mat_file* <=> spmmat (a string representing an existing file):
    The SPM.mat file that contains the estimate model.

    ::
        ex. '/home/username/data/derived_data/SPM.mat'

**Optional inputs parameters with a default value:**

- *session_type* <=> ?? (a string among "tcon", "fcon", "tconsess")


    ::
        ex. 'tcon'

- *contrast_name* <=> ?? (a list of string, optional, default value is ["+"]):
    Names of contrast.

    ::
        ex. ["Task vs control", "Task vs baseline"]

- *condition_name* <=> ?? (a list of list of string, optional, default value is [["R1_1"]]):
    Names of conditions used for each contrast.

    ::
        ex. [[]]

- *contrast_weight* (a list of list of float, default value is [[1.0]])
    List of contrast weights.

    ::
        ex. [[1.0]]

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

**Optional inputs:**

- *session* (a list of list of float, optional):
    Session list. Length should be equal to the number of sessions, with 1. for sessions to include and 0. elsewhere
    Default is Undefined.

        ::
            ex. ["0", "1.0"]

**Outputs parameters:**

- *out_spm_mat_file* (a pathlike object or string representing a file):
    The SPM.mat file containing specification of the design, estimated model parameters and contrast.
    Be carefull the input smp_mat_file will be overwritten.

    ::
        ex. '/home/username/data/derived_data/SPM.mat'

- *con_images* (a list of items which are a pathlike object or string representing a file):
    Contrast images from a t-contrast.

    ::
        ex. ['/home/username/data/derived_data/con_0001.nii',
            '/home/username/data/derived_data/con_0002.nii',
            '/home/username/data/derived_data/con_0003.nii',
            '/home/username/data/derived_data/con_0004.nii']

- *spmT_images* (a list of items which are a pathlike object or string representing a file):
    Stat images from a t-contrast.

    ::
        ex. ['/home/username/data/derived_data/spmT_0001.nii',
            '/home/username/data/derived_data/spmT_0002.nii',
            '/home/username/data/derived_data/spmT_0003.nii',
            '/home/username/data/derived_data/spmT_0004.nii']


-------------

.. [#label] Syntax: mia_processes/nipype EstimateConstrast brick <=> SPM12 fMRI contrast manager.

	    Usefull links:
	    `nipype EstimateContrast <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.model.html#estimatemodel>`_
