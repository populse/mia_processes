:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Despike brick
============

Removes spikes from the 3D+time input dataset using AFNI 3DDespike command. 
The spike values replaced by something more pleasing to the eye.
The output dataset will always be stored in floats.

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import Despike

>>> Despike.help()

**Inputs parameters:**

- *in_file* (a string representing an existing file)
    Input file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

- *despike* (a boolean, optional)
    Despike dataset. 
    Default is True. 

    ::

      ex. True

- *output_type* (NIFTI or NIFTI_GZ, optional)
    | Format of the output image (one of NIFTI, NIFTI_GZ, AFNI).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *out_prefix* (a string, optional)
    Prefix of the output image. Default is 'd_'.
    
    ::

      ex. 'autpmd_ask_'


**Outputs parameters:**

- *out_file* (a strings representing a file)
    Despike file (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/d_func.nii'

-------------

Usefull links:

`AFNI 3dDespike <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html>`_
`AFNI Despike - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.preprocess.html#despike>`_
