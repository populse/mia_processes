:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==================
MRGridRegrid brick
==================

Perform changes of the voxel grid that require interpolation of the image
--------------------------------------------------------------------------

Performs changes of the voxel grid such as changing the resolution or location and orientation of the voxel grid.
If the image is down-sampled, the appropriate smoothing is automatically applied using Gaussian smoothing
unless nearest neighbour interpolation is selected or oversample is changed explicitly.

The resolution can only be changed for spatial dimensions.

(mrtrix mrgrid regrid command)

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

- *interp* (a string among cubic, nearesr, linear, sinc, defaunt value is cubic, optional)
    Interpolation method to use when reslicing.

    ::

      ex. 'cubic'


**Optional inputs:**

- *suffix* (a string, optional)
    Suffix of output image

    ::

      ex. 'regrid'

- *voxel* (a float or an integer or a list of 3 integer, optional)
    The new voxel size for the output image. his can be specified either as a single value (int or float)
    to be used for all spatial dimensions or as list of three integer for each voxel dimension.

    ::

      ex. 3.0

- *size* (a list of 3 integer, optional)
    The new size (number of voxels) in each spatial dimension for the output image.

    ::

      ex. [120, 120, 120]

- *scale* (an integer or a float or a list of 3 integer or float, optional)
    Scale the image resolution by the supplied factor. This can be specified either as a single value
    to be used for all dimensions or as list of three integer or float for each voxel dimension.

    ::

      ex. 2.3

- *oversample* (an integer or a list of 3 integer, optional)
    Amount of over-sampling (in the target space) to perform when regridding.
    This can consist of a single integer, or list of 3 integers
    if different oversampling factors are desired along the different axes.
    Default is determined from ratio of voxel dimensions (disabled for nearest-neighbour interpolation).

    ::

      ex. 2

- *template* (an existing file, optional)
    A reference image, the output image will match the reference image grid (voxel spacing, image size, header transformation).
    ::

      ex. '/home/username/data/downloaded_data_data/template.mif


**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    The output regrid image

    ::

      ex. '/home/username/data/derived_data/DWI_regrid.mif'

-------------

Useful links:

`mrtrix mrgrid <https://mrtrix.readthedocs.io/en/latest/reference/commands/mrgrid.html>`_
