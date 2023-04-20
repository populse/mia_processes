:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============================
FWHMx brick
============================

Computes FWHMs using AFNI FWHMx command.
This program computes FWHMs for all sub-bricks in the input dataset, each one separately

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image.

    ::

      ex. '/home/username/data/derived_data/mean_reg_func_valid.nii'

**Optional inputs with default value parameters:**

- *combine* (a boolean, optional, default value is True)
    Combine the final measurements along each axis.

    ::

      ex. True

- *detrend* (a boolean or an integer, optional, default value is True)
    Detrend to the specified order.

    ::

      ex. True

- *out_prefix* (a string, optional, default value is 'fwhm')
    Specify the string to be prepended to the filename of the output file.

    ::

      ex. 'fwhm_'

**Optional inputs:**

- *args* (a string representing an existing file, optional)
    Add an option as '-ShowMeClassicFWHM' (use old classic method to compute FWHM)

    ::

      ex. '-ShowMeClassicFWHM'

- *mask_file* (a string representing an existing file, optional)
    Mask image.

    ::

      ex. '/home/username/data/derived_data/automask_mean_reg_func_valid.nii'


**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out file.

    ::

      ex. '/home/username/data/derived_data/fwhm_mean_reg_func_valid.out'

-------------

Usefull links:

`AFNI FWHMx <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dFWHMx.html>`_
`AFNI FWHMx - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.utils.html#fwhmx>`_
