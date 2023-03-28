:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
AffineInitializer brick
============

A multi-start optimizer for affine registration.

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import ApplyTransforms

>>> ApplyTransforms.help()

**Inputs parameters:**

- *input_image* (a string representing an existing file)
    Input file to apply transformation to (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *reference_image* (a string representing an existing file)
    Input file to apply transformation to (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *transforms* (multiples string representing existing files)
    Transform files that will be applied in reverse order.

    ::

      ex. ['identity', 'ants_Warp.nii.gz', 'trans.mat']

- *interpolation*(Linear or NearestNeighbor or CosineWindowedSinc or WelchWindowedSinc or HammingWindowedSinc or LanczosWindowedSinc or MultiLabel or Gaussian or BSpline, optional)
    Interpolation model.
    Default is Linear.
    ::

      ex. 'Linear'

- *default_value* (a float, optional)
    Default value. Default is 0.0.
    ::

      ex. 0.0

- *dimension* (2 or 3 or 4, optional)
    This option forces the image to be treated as a specified-dimensional image (2 or 3 or 4).
    Default is 3.
    ::

      ex. 3

- *float* (a boolean, optional)
    Use float instead of double for computations.
    Default is False.
    ::

      ex. False

- *input_image_type* (0 or 1 or 2 or 3, optional)
    Option specifying the input image type of scalar (0), vector (1), tensor(2), or time series(3). 
    Default is 0 (scalar).
    ::

      ex. 0

- *interpolation_parameters* (a tuple of the form: (an integer) or a tuple of the form: (a float, a float), optional)
    Interpolation parameters. For example for BSpline order or for sigma/alphaparameters for Gaussian/Multilabel  
    Default is Undefined (ie parameter not used).
    ::

      ex. (5,) (if interpolation = 'BSpline')

- *invert_transform_flags* (a list of boolean, optional)
    Invert transform flags. 
    For example if transforms parameter is equal to ['ants_Warp.nii.gz', 'trans.mat'] and invert_transform_flags is [False, True], the follinw trnasform will be applied: 
    `--transform ants_Warp.nii.gz --transform [ trans.mat, 1 ]`
    Default is Undefined (ie parameter not used).
    ::

      ex. [False, True]

- *print_out_composite_warp_file* (a boolean, optional)
    Output a composite warp file instead of a transformed image.
    Default is False.
    ::

      ex. False

- *out_prefix* (a string, optional)
    Prefix of the output image. Default is 't_'.
    
    ::

        ex. 't_'

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/denoise_func.nii'

-------------

Usefull links:

`ANTS ApplyTransforms - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.ants.html#applytransforms>`_

