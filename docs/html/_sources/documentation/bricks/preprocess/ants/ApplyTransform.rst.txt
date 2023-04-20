:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
ApplyTransforms brick
============

ApplyTransforms, applied to an input image, transforms it according to a reference image and a transform (or a set of transforms)
using ANTs ApplyTransforms command.

--------------------------------------

**Mandatory inputs parameters:**

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

**Optional inputs parameters with a default value:**

- *default_value* (a float, optional, default value is 0.0)
    ex.
    ::

      ex. 0.0

- *dimension* (2 or 3 or 4, optional, default value is 3)
    This option forces the image to be treated as a specified-dimensional image (2 or 3 or 4).

    ::

      ex. 3

- *float* (a boolean, optional, default value is False)
    Use float instead of double for computations.

    ::

      ex. False

- *input_image_type* (0 or 1 or 2 or 3, optional, default value is 0)
    Option specifying the input image type of scalar (0), vector (1), tensor(2), or time series(3).

    ::

      ex. 0

- *interpolation*(Linear or NearestNeighbor or CosineWindowedSinc or WelchWindowedSinc or HammingWindowedSinc or LanczosWindowedSinc or MultiLabel or Gaussian or BSpline, optional, default value is Linear)
    Interpolation model.

    ::

      ex. 'Linear'

- *out_prefix* (a string, optional, default value is 't')
    Prefix of the output image.

    ::

        ex. 't_'

- *print_out_composite_warp_file* (a boolean, optional, default value is False)
    Output a composite warp file instead of a transformed image.
    Default is False.
    ::

      ex. False

**Optional inputs parameters:**

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

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out image (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/denoise_func.nii'

-------------

Usefull links:

`ANTS ApplyTransforms - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.ants.html#applytransforms>`_
