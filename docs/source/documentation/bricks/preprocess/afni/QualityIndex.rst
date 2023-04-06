
+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===================
QualityIndex brick
===================

Computes a quality index for each sub-brick in a 3D+time dataset using AFNI 3dTqual command.
The output is a 1D time series with the index for each sub-brick 

--------------------------------------

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

- *automask* (a boolean, optional, default value is True)
    Clip off small voxels. Mutually exclusive with mask

    ::

      ex. True

- *interval* (a boolean or an integer, optional, default value is False)
    Write out the median + 3.5 MAD of outlier count with each timepoint.

    ::

      ex. False

- *out_prefix* (a string, optional, default value is 'QI')
    Specify the string to be prepended to the filename of the output file.

    ::

      ex. 'QI_'

- *quadrant* (a boolean or an integer, optional, default value is False)
    Similar to spearman parameter, but using 1 minus the quadrant correlation coefficient as the quality index.

    ::

      ex. False

- *spearman* (a boolean or an integer, optional, default value is False)
    Quality index is 1 minus the Spearman (rank) correlation coefficient of each sub-brick with the median sub-brick.

    ::

      ex. False


**Optional inputs:**

- *mask_file* (a string representing an existing file, optional)
    Mask image. Compute correlation only across masked voxels. Mutually exclusive with automask and autoclip.

    ::

      ex. '/home/username/data/derived_data/automask_mean_reg_func_valid.nii'

- *polort* (an integer, optional)
    Detrend each voxel timeseries with polynomials.
    Default value is Undefined (i.e parameter not used)

    ::

      ex. 3

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out file.
    
    ::

      ex. '/home/username/data/derived_data/outliers_reg_func_valid.out'

-------------

Usefull links:

`AFNI 3dTqual <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTqual.html>`_ 
`AFNI QualityIndex - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.preprocess.html#qualityindex>`__
