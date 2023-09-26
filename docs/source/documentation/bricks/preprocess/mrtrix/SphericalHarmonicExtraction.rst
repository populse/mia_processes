:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=================================
SphericalHarmonicExtraction brick
=================================

Extract the peaks of a spherical harmonic function in each voxel
----------------------------------------------------------------

(mrtrix sh2peaks command)

*Please note that, in mia_processes, MRtrix bricks required diffusion data in MRtrix .mif format.
The MRConvert brick could be used to convert diffusion data in NIfTI format into MRtrix .mif format.
In this case, bvec and bval file should be specified.*

--------------------------------------

**Mandatory inputs parameters:**

- *in_SH_coeff* (a pathlike object a string representing an existing file)
    The input image of SH coefficients (valid extensions: [.mif]).

    ::

      ex. '/home/username/data/raw_data/wm_fod.mif'


**Optional inputs with default value parameters:**

- *suffix* (a string, default value is peaks, optional)
    Output file suffix

    ::

      ex. 'peaks'


- *num* (an integer, default value is 3, optional)
    The number of peaks to extract

    ::

      ex. 3

- *fast* (a boolean, default value is False, optional)
    Use lookup table to compute associated Legendre polynomials (faster, but approximate)

    ::

      ex. False

**Optional inputs parameters:**

- *direction* (a tuple of the form (a Float, a Float))
    The direction of a peak to estimate (phi, theta)

    ::

      ex. (30, 60)

- *peaks_image* (string representing an existing file, optional)
    The program will try to find the peaks that most closely match those in the image provided.

    ::

      ex. '/home/username/data/derived_data/peaks.mif'

- *thresh_value* (an float, optional)
    Only peak amplitudes greater than the threshold will be considered.

    ::

      ex. 2.3

- *seeds_file* (string representing an existing file, optional)
    Specify a set of directions from which to start the multiple restarts of the optimisation

    ::

      ex. '/home/username/data/derived_data/seeds.mif'


- *seeds_file* (string representing an existing file, optional)
    Only perform computation within the specified binary brain mask image

    ::

      ex. '/home/username/data/derived_data/brain_mask.mif'


**Outputs parameters:**

- *output_image* (a pathlike object or string representing a file)
    The output image. Each volume corresponds to the x, y & z component of each peak direction vector in turn.

    ::

      ex. '/home/username/data/derived_data/wm_fod_peaks.mif'

-------------

Usefull links:

`mrtrix sh2peaks <https://mrtrix.readthedocs.io/en/latest/reference/commands/sh2peaks.html>`_
