:orphan:

.. toctree::

+--------------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+-------------------------------------------+----------------------------------------------------+

===========================
FramewiseDisplacement brick
===========================

Calculate the FD (framewise displacement) as in `[Power2012] <https://doi.org/10.1016/j.neuroimage.2011.10.018>`_
-----------------------------------------------------------------------------------------------------------------

This implementation reproduces the calculation in fsl_motion_outliers.


Adapted from `nipype Cofunds <https://github.com/nipy/nipype/blob/f662acfce8def4717e0c3414618f3a5de5913b31/nipype/algorithms/confounds.py#L298>`_

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Motion parameters file.

    ::

      ex. '/home/username/data/derived_data/reg_func_valid_oned.txt'

**Optional inputs with default value parameters:**

- *normalize* (a boolean, optional, default value is False)
    Calculate FD in mm/s.

    ::

      ex. False

- *out_prefix* (a string, optional, default value is 'fd')
    Specify the string to be prepended to the filename of the output file.

    ::

      ex. 'fd_'

- *parameter_source* (a string which is FSL or AFNI or SPM or FSFAST or NIPY', optional, default value is FSL)
    Source of movement parameters.

    ::

      ex. FSL

- *radius* (a float, optional, default value is 50.0)
    Radius in mm to calculate angular FDs.

    ::

      ex. 50.0


**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out file.

    ::

      ex. '/home/username/data/derived_data/fd_reg_func_valid_oned.out'

-------------

Useful links:

`nipype Cofunds <https://github.com/nipy/nipype/blob/f662acfce8def4717e0c3414618f3a5de5913b31/nipype/algorithms/confounds.py#L298>`_

`mriqc IQMS <https://mriqc.readthedocs.io/en/22.0.6/iqms/bold.html>`_
