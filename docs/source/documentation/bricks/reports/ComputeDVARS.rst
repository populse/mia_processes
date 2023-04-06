:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

====================
ComputeDVARS brick
====================

Computes DVARS (derivative of variance). 

The average change in mean intensity between each pair of fMRI volumes in a series.
Higher values indicate more dramatic changes (e.g., due to motion or spiking).

Adapted from `nipype Cofunds <https://github.com/nipy/nipype/blob/f662acfce8def4717e0c3414618f3a5de5913b31/nipype/algorithms/confounds.py#L100>`_.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input functional image  after head motion correction (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/reg_func_valid.nii'

- *in_mask* (a string representing an existing file)
    Brain mask (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/automask_mean_reg_func_valid.nii'

**Optional inputs with default value parameters:**

- *intensity_normalization* (a float, optional, default value is 1000)
    Divide value in each voxel at each' timepoint by the median calculated across all voxels and timepoints
    within the mask (if specified) and then multiply by the value specified by this parameter.
    
    By using the default (1000) output DVARS will be expressed in x10 % BOLD units compatible with Power et al.2012.'
    Set this to 0 to disable intensity normalization altogether.

    ::

      ex. 1000

- *out_prefix* (a string, optional, default value is 'dvars')
    Specify the string to be prepended to the filename of the output file.

    ::

      ex. 'dvars_'

- *remove_zero_variance* (a boolean, optional, default value is True)
    Remove voxels with zero variance.

    ::

      ex. True

- *variance_tol* (a float, optionnal, default value is 0.0000001000)
    Maximum variance to consider "close to" zero for the purposes of removal.

    ::

      ex. 0.0000001000


**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out file with DVARS information.
    
    ::

      ex. '/home/username/data/derived_data/dvars_reg_func_valid.out'

-------------

Usefull links:

`nipype Cofunds <https://github.com/nipy/nipype/blob/f662acfce8def4717e0c3414618f3a5de5913b31/nipype/algorithms/confounds.py#L100>`_ 
`mriqc IQMS <https://mriqc.readthedocs.io/en/22.0.6/iqms/bold.html>`_
