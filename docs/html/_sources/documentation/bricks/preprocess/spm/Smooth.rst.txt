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

>>> from mia_processes.bricks.preprocess.spm import Smooth
>>> Smooth.help()

**Inputs parameters:**

- *in_files <=> data* [#label]_
    List of files to smooth. A list of items which are an existing, uncompressed file (valid extensions: [.img, .nii, .hdr]).

    ::

      ex. ['/home/ArthurBlair/data/raw_data/Func.nii']

- *fwhm <=> fwhm* [#label]_
    Specify the full-width at half maximum (FWHM) of the Gaussian smoothing kernel in mm (a float or a list of 3 items which are a float).
    Three values can be entered, indicating the FWHM in the x, y and z directions, or one value only, indicating the same FWHM in all three
    directions.

    ::

      ex. [6, 6, 6]

- *data_type <=> dtype* [#label]_
    | Data type of the output images (an integer [int or long]).
    |   0: same as the original images
    |   2: UINT8 (unsigned char)
    |   4: INT16 (signed short)
    |   6: INT32 (signed int)
    |   8: FLOAT32 (single prec. float)
    |   10: FLOAT64 (double prec. float)

    ::

      ex. 0

- *implicit_masking <=> im* [#label]_
    A mask implied by a particular voxel value (a boolean). If set to True, the implicit masking of the input image is preserved in the
    smoothed image.

    ::

      ex. False

- *out_prefix <=> prefix* [#label]_
    Specify the string to be prepended to the filenames of the smoothed image file(s) (a string).

    ::

      ex. s

**Outputs parameters:**

- *smoothed_files*
    Smoothed files (a list of items which are an existing file name).

    ::

      ex. /home/ArthurBlair/data/derived_data/sFunc.nii

-------------

.. [#label] Syntax: mia_processes/nipype Smooth <=> SPM12 Smooth.

	    Usefull links:
	    `SPM12 Smooth <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=57>`_, 
	    `nipype Smooth <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.preprocess.html#smooth>`_
