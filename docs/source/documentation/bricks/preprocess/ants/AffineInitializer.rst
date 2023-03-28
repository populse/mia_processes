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

>>> from mia_processes.bricks.preprocess.others import AffineInitializer

>>> AffineInitializer.help()

**Inputs parameters:**

- *moving_image* (a string representing an existing file)
    The moving image to be mapped to the fixed space(valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *reference_image* (a string representing an existing file)
    The fixed reference image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w_template.nii'

- *dimension* (3 or 2, optional)
    Image dimension. Default is 3.
    
    ::

        ex. 3

- *local_search* (an integer, optional)
    Determines if a local optimization is run at each search point for the set number of iterations. Default is 10.
    
    ::

        ex. 10

- *search_factor* (a float, optional)
    Increments (degrees) for affine search. Default is 15.0.
    
    ::

        ex. 15.0


- *search_factor* (a float between 0.0 and 1.0, optional)
    Search this arc +/- principal axes. Default is 0.1.
    
    ::

        ex. 0.1

- *principal_axes* (a boolean, optional)
    Whether the rotation is searched around an initial principal axis alignment. Default is False.
    
    ::

        ex. False

- *out_prefix* (a string, optional)
    Prefix of the out transform file. Default is 'AffineTransform_'.
    
    ::

        ex. 'AffineTransform_'

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Output transform file  (extensions: [.mat]).
    
    ::

      ex. '/home/username/data/derived_data/AffineTransform_T1w.mat'

-------------

Usefull links:

`ANTS AffineInitializer - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.ants.html#affineinitializer>`_


