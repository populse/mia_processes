:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===============
TStatMean brick
===============

Mean of bold images (using mean option of the AFNI 3dTstat)
-----------------------------------------------------------

Compute mean of input voxels for a 3D+time dataset image.

---------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input bold file to be averaged (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

**Optional inputs with default value parameters:**

- *output_type* (NIFTI or NIFTI_GZ, optional, default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   - NIFTI: \*.nii
    |   - NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *out_prefix* (a string, optional, default value is 'mean')
    Prefix of the output image.

    ::

        ex. 'mean_'

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out image (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/mean_func.nii'

-------------

Useful links:

`AFNI 3dTstat <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTstat.html>`_

`AFNI TStat - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.utils.html#tstat>`_
