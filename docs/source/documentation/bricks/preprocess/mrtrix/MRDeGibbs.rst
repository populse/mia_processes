:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===============
MRDeGibbs brick
===============

Remove Gibbs ringing artifacts
-------------------------------

Remove Gibbs ringing artefacts from MRI images using the method of local subvoxel-shifts proposed by Kellner et al.

This brick should be run before any interpolation and before motion correction (before brick DWIPreproc)

(mrtrix mrdegibbs command)

*Please note that, in mia_processes, MRtrix bricks required diffusion data in MRtrix .mif format.
The MRConvert brick could be used to convert diffusion data in NIfTI format into MRtrix .mif format.
In this case, bvec and bval file should be specified.*

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a pathlike object or a string representing an existing file)
    Input image (valid extensions: [.mif]).

    ::

      ex. '/home/username/data/raw_data/DWI.mif

**Optional inputs with default value parameters:**

- *axes* (a list of items which are an integer, default value is [0, 1], optional)
    | Indicate the plane in which the data was acquire:
    |   - [0, 1]: axial
    |   - [0, 2] :coronal
    |   - [1, 2]: sagittal

    ::

      ex. [0, 1]


- *maxW* (an integer, default value is 3, optional)
    Right border of window used for total variation (TV) computation

    ::

      ex. 3

- *minW* (an integer, default value is 1, optional)
    Left  border of window used for total variation (TV) computation

    ::

      ex. 1

- *nshifts* (an integer, default value is 20, optional)
    Discretization of subpixel spacing

    ::

      ex. 20


**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    The output unringed DWI image

    ::

      ex. '/home/username/data/derived_data/DWI_unringed.mif'

-------------

Useful links:

`mrtrix mrdegibbs <https://mrtrix.readthedocs.io/en/latest/reference/commands/mrdegibbs.html>`_

`mrtrix mrdegibbs - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.preprocess.html#mrdegibbs>`_

`E.Kellner - Gibbs-ringing artifact removal based on local subvoxel-shifts <https://pubmed.ncbi.nlm.nih.gov/26745823/>`_
