:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==================
Registration brick
==================

Image registration using ANTs Registration command. 

This brick registers a moving_image to a fixed_image, using a predefined (or a sequence of) cost function(s) 
and transformation operations. The cost function is defined using one or more ‘metrics’. 

Both linear and non-linear transformations can be used. Usually, registration is done in multiple stages. 
For example first an Affine, then a Rigid, and ultimately a non-linear (Syn)-transformation (see 'transform' and 'transform_parameters' parameter).

It is possible to initilize the registration by using one or more transforms from moving_image to fixed_image with the 
'initial_moving_transform' parameter. 
For example, when you already have a warpfield that corrects for geometrical distortions in an EPI (functional) image, that you want to apply before an Affine registration to a structural image. 
You could put this transform into 'intial_moving_transform'.

Here can be found some of classical registration parameters(used in fMRIPrep and MRIQC):
`niworkflows parameters examples <https://github.com/nipreps/niworkflows/tree/master/niworkflows/data>`_ 

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import Registration

>>> Registration.help()

**Mandatory inputs parameters:**

- *fixed_image* (a string representing an existing file)
    Image to which the moving image should be transformed (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/template.nii'

- *moving_image* (a string representing an existing file)
    Image that will be registered to the space of the fixed image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'


**Optional inputs parameters with a default value:**

- *collapse_output_transforms* (a boolean, optional)
    Collapse output transforms.
    
    ::

        default value. True

- *dimension* (2 or 3 or 4, optional)
    Image dimension (2 or 3).

    ::

        default value. 3

- *float* (a boolean, optional)
    Use float instead of double for computations.
    
    ::

        default value. False

- *initialize_transforms_per_stage* (a boolean, optional)
    Initialize linear transforms from the previous stage.

    ::

        default value. False

- *interpolation* (Linear or NearestNeighbor or CosineWindowedSinc or WelchWindowedSinc or HammingWindowedSinc or LanczosWindowedSinc or BSpline or MultiLabel or Gaussian or GenericLabel, optional)
    Interpolation model.
    
    ::

        default value. 'Linear'


- *metrics* (a list of items which are CC or MeanSquares or Demon or GC or MI or mattes, optional)
    | The metric to use for each stage.
    | Metrics available:
    |   cross-correlation (CC)
    |   Mean Squares (MeanSquares)
    |   Demons (Demons)
    |   global correlation (GC)
    |   Mutual Information (Mattes or MI)
   
    ::

        default value. ['Mattes', 'Mattes', 'Mattes']

- *metric_weight* (a list of items which are a float, optional)
    The metric weight(s) for each stage. The weights must sum to 1 per stage. 
    
    ::

        default value. [1.0, 1.0, 1.0]

- *output_inverse_warped_image* (a boolean, optional)
    Get output inverse_warped_image.
    Default is False.
    
    ::

        default value. False

- *out_prefix* (a string, optional)
    Prefix of the warped output image.
    
    ::

        default value. 'w_'

- *radius_bins_item_trait* (an integer, optional)
    Radius bins item.
    
    ::

        default value. 5

- *signa_units* (a list of items which are vox or mm, optional)
    Units for smoothing sigmas (mm or vox).
    
    ::

        default value. ['vox'] * len(metric)

- *smoothing_sigmas* (a list of items which are a list of items which are a float, optional)
    Smoothing sigmas.
    
    ::

        default value. [[4.0], [4.0, 2.0, 0.0], [1.0, 0.0]]

- *shrink_factors* (a list of items which are a list of items which are an integer, optional)
    Shrink factor.
    
    ::

        default value. [[4], [4, 2, 1], [2, 1]]

- *transforms* (a list of items which are Rigid or Affine or CompositeAffine or Similarity or Translation or BSpline or GaussianDisplacementField or TimeVaryingVelocityField or TimeVaryingBSplineVelocityField or SyN or BSplineSyN or Exponential or BSplineExponential, optional)
    | List of the transformations.
    |   Linear transformations:
    |     Translation
    |     Rigid
    |     Affine
    |     CompositeAffine
    |     Translation
    |   Non-linear transformations:
    |     BSpline
    |     GaussianDisplacementField
    |     TimeVaryingVelocityField
    |     TimeVaryingBSplineVelocityField
    |     SyN
    |     BSplineSyN
    |     Exponential
    |     BSplineExponential
    
    ::

        default value. ['Rigid', 'Affine', 'SyN']

- *transform_parameters* (a list of items which are a tuple of the form: (a float) or a tuple of the form: (a float, a float, a float) or a tuple of the form: (a float, an integer, an integer, an integer) or a tuple of the form: (a float, an integer, a float, a float, a float, a float) or a tuple of the form: (a float, a float, a float, an integer) or a tuple of the form: (a float, an integer, an integer, an integer, an integer), optional)
    Transforms parameters.
    
    ::

        default value. [(0.01,), (0.08,), (0.1, 3.0, 0.0)]

- *use_estimate_learning_rate_once* (a list of boolean, optional)
    Use estimate learning rate.
    
    ::

        default value. [True] * len(metric)

- *use_histogram_matching* (a list of boolean, optional)
    Use histogram matching.
    
    ::

        default value. [True] * len(metric)

- *winsorize_lower_quantile* (a float between 0.0 and 1.0, optional)
    The lower quantile to clip image ranges.
    
    ::

        default value. 0.005

- *winsorize_upper_quantile* (a float between 0.0 and 1.0, optional)
    The upper quantile to clip image ranges.
    
    ::

        default value. 0.995

- *write_composite_transform* (a boolean, optional)
    Write composite transform.
    
    ::

        default value. True

**Optional inputs parameters:**

- *convergence_threshold* (a list of items which are a float, optional)
    Convergence threshold.
    Requires number_of_iterations parameter.
    Default is Undefined (ie parameter not used).
    
    ::

        ex. [1e-07, 1e-08, 1e-8]

- *convergence_window_size* (a list of items which are an integer, optional)
    Convergence window size.
    Default is Undefined (ie parameter not used).
    
    ::

        ex. [15, 5, 3]

- *fixed_image_masks* (a string representing an existing file or a list of strings representing existing files or NULL, optional)
    Mask used to limit metric sampling region of the fixed image defined per registration stage (valid mask extensions: [.nii, .nii.gz]).
    If a list of items is used, use “NULL” to omit a mask at a given stage. 
    Default is NULL (ie option not used).

    ::

      ex. '/home/username/data/raw_data/template_desc-brain_mask.nii'

- *initial_moving_transform* (multiples string representing existing files, optional)
    A transform or a list of transform that should be applied before the. Mutually exclusive with initial_moving_transform_com.
    
    ::

      ex. ['trans.mat']

- *initial_moving_transform_com* (0 or 1 or 2, optional)
    Align the moving_image and fixed_image before registration using the geometric center of the images (=0), 
    the image intensities (=1), or the origin of the images (=2). Mutually exclusive with initial_moving_transform.
    
    ::

      ex. 0

- *interpolation_parameters* (a tuple of the form: (an integer) or a tuple of the form: (a float, a float) or a tuple of the form: (a string), optional)
    Interpolation parameters. For example for BSpline order or for sigma/alphaparameters for Gaussian/Multilabel  
    Default is Undefined (ie parameter not used).

    ::

      ex. (5,) (if interpolation = 'BSpline')

- *invert_initial_moving_transform* (a list of boolean)
    A list of booleans that indicatewhether the inverse(s) of the transform(s) defined in initial_moving_transform should be used. 
    Requires initial_moving_transform parameter.
    Default is Undefined (ie parameter not used).
   
    ::

      ex. [False]

- *moving_image_mask* (a string representing an existing file or a list of strings representing existing files or NULL, optional)
    Masks used to limit metric sampling region of the moving image, defined per registration stage (valid extensions: [.nii, .nii.gz]).
    If a list of items is used, use “NULL” to omit a mask at a given stage. 
    Default is NULL (ie option not used)

    ::

      ex. '/home/username/data/raw_data/T1w_desc-brain_mask.nii'

- *number_of_iterations* (a list of items which are a list of item which are an integer, optional)
    Number of iterations.
    Default is Undefined (ie parameter not used).
    
    ::

        ex. [[10000, 1000, 100], [500, 250, 100], [100, 30, 20]]

- *radius_or_number_of_bins* (a list of items which are an integer, optional)
    The number of bins in each stage for the MI and Mattes metric, the radius for other metrics.
    Default is Undefined (ie parameter not used).
    
    ::

        ex. [56, 56, 56]

- *random_seed* (an integer, optional)
    Fixed seed for random number generation.
    Default is Undefined (ie parameter not used).
    
    ::

        ex. 5

- *sampling_percentage* (a list of items which are a float between 0.0 and 1.0, optional)
    The metric sampling percentages to use for each stage. 
    Requires sampling strategy parameter. 
    Default is Undefined (ie parameter not used).
    
    ::

        ex.[0.2, 0.1, 0.1]

- *sampling_strategy* (a list of items which are Random, Regular or None, optional)
    The metric sampling strategies for each stage.
    Default is Undefined (ie parameter not used).
    
    ::

        ex.['Random', 'Random', 'Random']


**Outputs parameters:**


- *composite_transform* (a strings representing a file)
    Composite transform file (extensions: [.h5]).
    
    ::

      ex. '/home/username/data/derived_data/T1w_Composite.h5'

- *inverse_composite_transform* (a strings representing a file)
    Inverse composite transform file (extensions: [.h5]).
    
    ::

      ex. '/home/username/data/derived_data/T1w_InverseComposite.h5'

- *inverse_warped_image* (a strings representing a file)
    Inverse warped image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/w_inverse_T1w.nii'

- *warped_image* (a strings representing a file)
    Warped image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/w_T1w.nii'

-------------

Usefull links:

`ANTs Registration - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.ants.html#registration>`_

