:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==================
OutlierCount brick
==================

Using AFNI 3dToutcount, computes outliers for all sub-bricks (3D volumes for each TR) in the input dataset
----------------------------------------------------------------------------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image.

    ::

      ex. '/home/username/data/derived_data/reg_func_valid.nii'

**Optional inputs with default value parameters:**

- *autoclip* (a boolean, optional, default value is False)
    Clip off small voxels. Mutually exclusive with mask.

    ::

      ex. False

- *automask* (a boolean, optional, default value is False)
    Clip off small voxels. Mutually exclusive with mask.

    ::

      ex. False

- *fraction* (a boolean, optional, default value is True)
    Combine the final measurements along each axis.

    ::

      ex. True

- *interval* (a boolean or an integer, optional, default value is False)
    Write out the median + 3.5 MAD of outlier count with each timepoint.

    ::

      ex. False

- *legendre* (a boolean or an integer, optional, default value is False)
    Use Legendre polynomials.

    ::

      ex. False

- *out_prefix* (a string, optional, default value is 'outliers')
    Specify the string to be prepended to the filename of the output file.

    ::

      ex. 'outliers_'

- *qthr* (a float between 0.0 and 1.0, default value is 0.001)
    Indicate a value for q to compute alpha.

    ::

      ex. 0.001

**Optional inputs:**

- *mask_file* (a string representing an existing file, optional)
    Mask image. Compute correlation only across masked voxels. Mutually exclusive with automask and autoclip.

    ::

      ex. '/home/username/data/derived_data/automask_mean_reg_func_valid.nii'

- *polort* (an integer, optional)
    Detrend each voxel timeseries with polynomials.

    ::

      ex. 3

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out file.

    ::

      ex. '/home/username/data/derived_data/outliers_reg_func_valid.out'

-------------

Useful links:

`AFNI 3dToutcount <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dToutcount.html>`_

`AFNI OutlierCount - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.preprocess.html#outliercount>`_
