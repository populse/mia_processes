:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
DWICat brick
============

Concatenating multiple DWI series
---------------------------------

(mrtrix dwicat command)

Concatenate 4D DWI images with an intensity scaling.

*Please note that, in mia_processes, MRtrix bricks required diffusion data in MRtrix .mif format.
The MRConvert brick could be used to convert diffusion data in NIfTI format into MRtrix .mif format.
In this case, bvec and bval file should be specified.*

--------------------------------------

**Mandatory inputs parameters:**

- *in_files* (a list of items which are a pathlike object a string representing an existing file)
    Inputs DWI files (valid extensions: [.mif]).

    ::

      ex. ['/home/username/data/derived_data/dwi1.mif', '/home/username/data/derived_data/dwi2.mif']


**Optional inputs parameters:**

- *in_mask* (a string representing an existing file, optional)
    A binary mask within which image intensities will be matched

    ::

      ex. '/home/username/data/derived_data/brainmask.mif'

- *out_name* (a string, optional)
    Out file name.
    If no out file name provide, the output file name will be:
    "concat_filename1_filename...""

    ::

      ex. 'new_concate_image'

**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    Output image (all DWI concatenated)

    ::

      ex. '/home/username/data/derived_data/concat_dwi1_dwi2.mif'


-------------

Usefull links:

`mrtrix dwicat <https://mrtrix.readthedocs.io/en/latest/reference/commands/dwicat.html>`_
