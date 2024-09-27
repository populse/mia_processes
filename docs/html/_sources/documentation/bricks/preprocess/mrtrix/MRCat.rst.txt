:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===========
MRCat brick
===========

Concatenate several images into one
------------------------------------

(mrtrix mrcat command)

*Please note that, in mia_processes, MRtrix bricks required diffusion data in MRtrix .mif format.
The MRConvert brick could be used to convert diffusion data in NIfTI format into MRtrix .mif format.
In this case, bvec and bval file should be specified.*

--------------------------------------

**Mandatory inputs parameters:**

- *in_files* (a list of items which are a pathlike object or a string representing an existing file)
    Input images to concatenate (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. ['/home/username/data/raw_data/DWI1.mif', '/home/username/data/raw_data/DWI2.mif']

**Optional inputs with default value parameters:**

- *axis* (an integer, default value is 3, optional)
    Specify axis along which concatenation should be performed

    ::

      ex. 3


- *out_file_name* (a string, default value is 'concatenated', optional)
    Out file name

    ::

      ex. 'concatenated'



**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    The output concatenated image.

    ::

      ex. '/home/username/data/derived_data/concatenated.mif'

-------------

Useful links:

`mrtrix mrcat <https://mrtrix.readthedocs.io/en/latest/reference/commands/mrcat.html>`_

`mrtrix mrcat - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.utils.html#mrcat>`_
