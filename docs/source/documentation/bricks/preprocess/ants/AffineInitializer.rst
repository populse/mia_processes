:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=======================
AffineInitializer brick
=======================

Initialize an affine transform using ANTs AffineInitizer command.

--------------------------------------

**Mandatory inputs parameters:**

- *moving_image* (a string representing an existing file)
    The moving image to be mapped to the fixed space (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *fixed_image* (a string representing an existing file)
    The fixed reference image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w_template.nii'

**Optional inputs parameters with a default value:**

- *dimension* (3 or 2, optional, default value is 3)
    Image dimension.

    ::

        ex. 3

- *local_search* (an integer, optional, default value is 10)
    Determines if a local optimization is run at each search point for the set number of iterations.

    ::

        ex. 10

- *out_prefix* (a string, optional, default value is 'AffineTransform')
    Prefix of the out transform file.

    ::

        ex. 'AffineTransform_'

- *principal_axes* (a boolean, optional, default value is False)
    Whether the rotation is searched around an initial principal axis alignment.

    ::

        ex. False

- *radian_fraction* (a float between 0.0 and 1.0, optional, default value is 0.1)
    Search this arc +/- principal axes.

    ::

        ex. 0.1

- *search_factor* (a float, optional, default value is 15.0)
    Increments (degrees) for affine search.

    ::

        ex. 15.0

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Output transform file  (extensions: [.mat]).

    ::

      ex. '/home/username/data/derived_data/AffineTransform_T1w.mat'

-------------

Usefull links:

`ANTS AffineInitializer - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.ants.html#affineinitializer>`_
