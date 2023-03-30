:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
AffineInitializer brick
============

Initialize an affine transform using ANTs AffineInitizer command.

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import AffineInitializer

>>> AffineInitializer.help()

**Mandatory inputs parameters:**

- *moving_image* (a string representing an existing file)
    The moving image to be mapped to the fixed space(valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *fixed_image* (a string representing an existing file)
    The fixed reference image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w_template.nii'

**Optional inputs parameters with a default value:**

- *dimension* (3 or 2, optional)
    Image dimension.
    
    ::

        default value. 3

- *local_search* (an integer, optional)
    Determines if a local optimization is run at each search point for the set number of iterations.
    
    ::

        default value. 10

- *out_prefix* (a string, optional)
    Prefix of the out transform file. Default is 'AffineTransform_'.
    
    ::
        default value. 'AffineTransform_'

- *principal_axes* (a boolean, optional)
    Whether the rotation is searched around an initial principal axis alignment.
    
    ::

        default value. False

- *radian_fraction* (a float between 0.0 and 1.0, optional)
    Search this arc +/- principal axes. 
    
    ::

        default value. 0.1

- *search_factor* (a float, optional)
    Increments (degrees) for affine search.
    
    ::

        default value. 15.0

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Output transform file  (extensions: [.mat]).
    
    ::

      ex. '/home/username/data/derived_data/AffineTransform_T1w.mat'

-------------

Usefull links:

`ANTS AffineInitializer - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.ants.html#affineinitializer>`_


