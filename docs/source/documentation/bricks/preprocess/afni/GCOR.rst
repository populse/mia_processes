:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==========
GCOR brick
==========

Using AFNI GCOR, computes the average correlation between every voxel and every other voxel, over any given mask
----------------------------------------------------------------------------------------------------------------

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

- *no_demean* (a boolean or an integer, optional, default value is False)
    Do not (need to) demean as first step.

    ::

      ex. False


**Optional inputs:**

- *mask_file* (a string representing an existing file, optional)
    Mask image.

    ::

      ex. '/home/username/data/derived_data/automask_mean_reg_func_valid.nii'

- *nfirst* (an integer, optional)
    Number of initial TRs to ignore. Default is Undefined (ie parameter not used)

    ::

      ex. 3



**Outputs parameters:**

- *out* (a float)
    Global correlation.

    ::

      ex. 0.06123

-------------

Usefull links:

`AFNI GCOR <https://afni.nimh.nih.gov/pub/dist/doc/program_help/@compute_gcor.html>`_

`AFNI GCOR - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.utils.html#gcor>`_
