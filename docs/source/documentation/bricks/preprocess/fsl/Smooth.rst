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

**Mandatory nputs parameters:**

- *in_file*
    An image to be smoothed. An item that is an existing, uncompressed file (valid extensions: [.img, .nii, .hdr]).

    ::

      ex. ['/home/username/data/raw_data/Func.nii']

**Optional inputs with default value parameters:**

- *fwhm* (a float, optional, default value is 6.0)
    Specify the full-width at half maximum of the Gaussian smoothing kernel in mm (not voxels).
    Will be converted internally to sigma in mm (not voxels), using approximately `2.355 * sigma = fwhm <https://brainder.org/2011/08/20/gaussian-kernels-convert-fwhm-to-sigma/>`_. Mutually exclusive with sigma.

    ::

      ex. 6.0

- *out_prefix* (a string, optional, default value is 's')
    Specify the string to be prepended to the filename of the smoothed image file.

    ::

      ex. 's_'

- *output_type* ('NIFTI' or 'NIFTI_GZ', optional, default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   - NIFTI: \*.nii
    |  -  NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

**Optional inputs parameters:**

- *sigma* (a float, optionnal)
    Specify the standard deviation of the filter in mm (not voxels). Mutually exclusive with fwhm.

    ::

      ex. 2.55

**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    Smoothed file.

    ::

      ex. '/home/username/data/derived_data/s_Func.nii'

-------------

Usefull links:

`FSL Smooth <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.fsl.utils.html#smooth>`_
