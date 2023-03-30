:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Smooth brick
============

3D Gaussian smoothing of image volumes
--------------------------------------

>>> from mia_processes.bricks.preprocess.fsl import Smooth
>>> Smooth.help()

**Mandatory nputs parameters:**

- *in_file*
    An image to be smoothed. An item that is an existing, uncompressed file (valid extensions: [.img, .nii, .hdr]).

    ::

      ex. ['/home/username/data/raw_data/Func.nii']

**Optional inputs with default value parameters:**

- *fwhm*
    Specify the full-width at half maximum of the Gaussian smoothing kernel (a float) in mm (not voxels).
    Will be converted internally to sigma in mm (not voxels), using approximately `2.355 * sigma = fwhm <https://brainder.org/2011/08/20/gaussian-kernels-convert-fwhm-to-sigma/>`_. Mutually exclusive with sigma.

    ::

      default value. 6.0

- *out_prefix*
    Specify the string to be prepended to the filename of the smoothed image file (a string).

    ::

      default value. s_

- *output_type*
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz

    ::

      default value. NIFTI

**Optional inputs parameters:**

- *sigma*
    Specify the standard deviation of the filter (a float) in mm (not voxels). Mutually exclusive with fwhm.

    ::

      ex. 2.55

**Outputs parameters:**

- *out_file*
    Smoothed file (an item which are a file name).

    ::

      ex. /home/username/data/derived_data/s_Func.nii

-------------

Usefull links:
`FSL Smooth <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.fsl.utils.html#smooth>`_