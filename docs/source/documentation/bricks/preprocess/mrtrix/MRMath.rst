:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
MRMath brick
============

Compute summary statistic on image intensities along a specified axis of a single image
---------------------------------------------------------------------------------------

(mrtrix mrmath command)

*Please note that, in mia_processes, MRtrix bricks required diffusion data in MRtrix .mif format.
The MRConvert brick could be used to convert diffusion data in NIfTI format into MRtrix .mif format.
In this case, bvec and bval file should be specified.*

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a pathlike object or a string representing an existing file)
    Input image (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/DWI.mif

**Optional inputs with default value parameters:**

- *operation* (mean or median or sum or product or rms or norm or var or std or min or max or absmax or magmax, default value is mean, optional)
    | Operation to computer along a specified axis:
    |   - mean
    |   - median
    |   - sum
    |   - product
    |   - rms
    |   - norm
    |   - var
    |   - sdt
    |   - min
    |   - max
    |   - absmax
    |   - absmin
    |   - magmax

    ::

      ex. mean


- *axis* (an integer, default value is 3, optional)
    Specify axis along which concatenation should be performed.
    Note that axis indices start from 0; thus, axes 0, 1 & 2 are the three spatial axes, and axis 3 operates across volumes.

    ::

      ex. 3

**Optional inputs:**

- *out_file_name* (a string, optional)
    Out file name. If not provided, name of the oeration is used.

    ::

      ex. "mean_along_axis3"


**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    The output image

    ::

      ex. '/home/username/data/derived_data/DWI_mean.mif'

-------------

Useful links:

`mrtrix mrmath <https://mrtrix.readthedocs.io/en/latest/reference/commands/mrmath.html>`_

`mrtrix mrmath - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.utils.html#mrmath>`_
