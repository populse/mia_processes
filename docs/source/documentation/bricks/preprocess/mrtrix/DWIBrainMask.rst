:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==================
DWIBrainMask brick
==================

Generates a whole brain mask from a DWI image.
----------------------------------------------

All diffusion weighted and b=0 volumes are used to obtain a mask that
includes both brain tissue and CSF.

(mrtrix dwi2mask command)

**Mandatory nputs parameters:**

- *in_file* (a pathlike object or a string representing an existing file)
    Input DWI image (valid extensions: [.mif]).

    ::

      ex. ['/home/username/data/raw_data/DWI.mif']


**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    Output brain mask.

    ::

      ex. '/home/username/data/derived_data/DWI_brainmask.mif'

-------------

Usefull links:

`mrtrix dwi2mask <https://mrtrix.readthedocs.io/en/latest/reference/commands/dwi2mask.html>`_
