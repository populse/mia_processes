:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
N4BiasFieldCorrection brick
============

Bias field correction using ANTs N4BiasFieldCorrection

N4 is a variant of the popular N3 (nonparameteric nonuniform normalization) retrospective bias correction algorithm. 
Based on the assumption that the corruption of the low frequency bias field can be modeled as a convolution of the intensity histogram by a Gaussian, 
the basic algorithmic protocol is to iterate between deconvolving the intensity histogram by a Gaussian, remapping the intensities, 
and then spatially smoothing this result by a B-spline modeling of the bias field itself

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import N4BiasFieldCorrection

>>> N4BiasFieldCorrection.help()

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    A image file (2D, 3D or 4D) to correct (valid extensions: [.nii, .nii.gz]).
    If thera are negative values or values close to zero, they are processed prior to correction.

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

**Optional inputs parameters with a default value:**

- *copy_header* (a boolean, optionnal)
    Copy headers of the original image into the output (corrected) file.

    ::

      default value. True

- *dimension* (2 or 3 or 4, optional)
    Image dimension (2 or 3 or 4).

    ::

      default value. 3

- *out_prefix* (a string, optional)
    Prefix of the output image.
    
    ::

        default value. 'n4c_'

- *rescale_intensities* (a boolean, optionnal)
    Rescales to the [min,max] range of theoriginal image intensities within the user-specified mask.

    ::

      default value. True

- *save_bias* (a boolean, optionnal)
    Save the estimated bias in a file .

    ::

      default value. True

**Optional inputs parameters:**

- *bspline_fitting_distance* (a float, optional)
    Set bspline fitting distance.
    Default is Undefined (ie parameter not used).

    ::

      ex. 300.0

- *bspline_order* (an integer, optional)
    Set bspline order. Require bspline_fitting_distance parameters.
    Default is Undefined (ie parameter not used).

    ::

      ex. 5

- *convergence_threshold* (a float, optional)
    Set convergence threshold.
    Require n_iterations parameters. 
    Default is Undefined (ie parameter not used).

    ::

      ex. 1e-6

- *histogram_sharpening* (a tuple of the form: (a float, a float, an integer), optional)
    Three-values tuple of histogram sharpening parameters (FWHM, wienerNose, numberOfHistogramBins). 
    These options describe the histogram sharpening parameters, i.e. the deconvolution step parameters described in the original N3 algorithm.
    Note that the default values have been shown to work fairly well.
    Default is Undefined (ie parameter not used).

    ::

      ex. (0.12, 0.02, 200)

- *mask_image* (a string representing an existing file, optional)
    Image to specify region to perform final bias correction in.
    Default is Undefined (ie parameter not used).

    ::

      ex. '/home/username/data/derived_data/T1w_mask.nii'

- *n_iterations* (a list of integer, optional)
    Set number of iterations.
    Default is Undefined (ie parameter not used).

    ::

      ex. [50, 50, 30, 20]

- *shrink_factor* (an integer, optional)
    Shrink factor.
    Default is Undefined (ie parameter not used).

    ::

      ex. 3

- *weight_image* (a string representing an existing file, optional)
    Image for relative weighting (e.g. probability map of the white matter) of voxels during the B-spline fitting.
    Default is Undefined (ie parameter not used).

    ::

      ex. '/home/username/data/derived_data/T1w_pve_0.nii'

**Outputs parameters:**

- *bias_image* (a strings representing a file)
    Estimated bias (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/T1w_bias.nii'
  
- *negative_values* (a strings representing a file)
    True if negative values are present in the input image. False otherwise.
    If thera are negative values or values close to zero, they are processed prior to correction.
    
    ::

      ex. False

- *out_file* (a strings representing a file)
    Out image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/n4c_T1w.nii'

-------------

Usefull links:

`ANTS N4BiasFieldCorrection - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.ants.html#n4biasfieldcorrection>`_

