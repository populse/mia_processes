:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
RefitDeoblique brick
============

Deoblique dataset (ie transform dataset from oblique to cardinal) using AFNI 3drefit command. 

Output file name is the same as input file name. 

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import RefitDeoblique

>>> RefitDeoblique.help()

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

**Optional inputs with default value parameters:**

- *deoblique* (a boolean, optional)
    Deoblique dataset. 

    ::

      default value. True


**Outputs parameters:**

- *out_file* (a strings representing a file)
    Deoblique file (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/func.nii'

-------------

Usefull links:

`AFNI 3drefit <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3drefit.html>`_
`AFNI Refit - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.utils.html#refit>`_
