:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=================
MRTransform brick
=================

Apply spatial transformations or reslice images
-----------------------------------------------

If a linear transform is applied without a template image, the image header transform matrix will be modified.

Fibre orientation distribution (FOD) reorientation (with apodised point spread functions) can be performed if the number of volumes in the 4th dimension equals the number of coefficients in an antipodally symmetric spherical harmonic series (e.g. 6, 15, 28 etc).
The fod_reorient should be specified.


(mrtrix mrtransform command)

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a pathlike object or a string representing an existing file)
    Input images to be transformed (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/DWI.mif

**Optional inputs with default value parameters:**

- *interpolation* (cubic, nearest, linear, sinc, default value is cubic, optional)
    Set the interpolation method to use when reslicing

    ::

      ex. cubic

- *inverse* (a boolean, default value is False, optional)
    Invert the specified transform before using it

    ::

      ex. False

- *half* (a boolean, default value is False, optional)
    Apply the matrix square root of the transformation

    ::

      ex. False

- *identity* (a boolean, default value is False, optional)
    Set the header transform of the image to the identity matrix

    ::

      ex. False

- *midway_space* (a boolean, default value is False, optional)
    Reslice the input image to the midway space. Requires either the template_image or warp_image option

    ::

      ex. False

- *fod_reorient* (a boolean, default value is False, optional)
    Specify whether to perform FOD reorientation

    ::

      ex. False


**Optional inputs:**

- *linear_transform* (a pathlike object or a string representing an existing file, optional)
    Specify a linear transform to apply. It should be a 3x4 or 4x4 ascii file.
    Note the standard reverse convention is used, where the transform maps points in the template image to the moving image.

    ::

      ex. '/home/username/data/derived_data/transform.txt

- *flix_axes* (a pathlike object or a string representing an existing file, optional)
    Flip the specified axes (a list of int with 0:x, 1:y and 2:z)

    ::

      ex. [0, 1, 2]

- *replace_file* (a pathlike object or a string representing an existing file, optional)
    Replace the linear transform of the original image by that specified, rather than applying it to the original image.
    The specified transform can be either a template image, or a 3x4 or 4x4 ascii file.

    ::

      ex. '/home/username/data/derived_data/template.nii'

- *template_image* (a pathlike object or a string representing an existing file, optional)
    Reslice the input image to match the specified template image.

    ::

      ex. '/home/username/data/derived_data/template.nii'

- *template_image* (a pathlike object or a string representing an existing file, optional)
    Reslice the input image to match the specified template image.

    ::

      ex. '/home/username/data/derived_data/template.nii'

- *oversample_factor* (an integer or a list of three integers, optional)
    Set the amount of over-sampling (in the target space) to perform

    ::

      ex. 3

- *warp_image* (a pathlike object or a string representing an existing file, optional)
    Apply a non-linear 4D deformation field to warp the input image

    ::

      ex. '/home/username/data/derived_data/deformation_field.nii'

- *warp_full_image* (a pathlike object or a string representing an existing file, optional)
    Warp the input image using a 5D warp file output from mrregister

    ::

      ex. '/home/username/data/derived_data/deformation_field_5D.nii'

- *fod_modulation* (fod or jac, optional)
    | Intensity modulation method for Fibre orientation distribution (fod):
    |   - fod: modulate FODs during reorientation to preserve the apparent fibre density across fibre bundle widths before and after the transformation.
    |   - jac: modulate the image intensity with the determinant of the Jacobian of the warp of linear transformation to preserve the total intensity before and after the transformation.

    ::

      ex. fod

- *fod_direction_file* (a pathlike object or a string representing an existing file, optional)
    Directions defining the number and orientation of the apodised point spread functions used in FOD reorientation.
    If not used, 300 directions are used

    ::

      ex. '/home/username/data/derived_data/direction.txt'


**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    The output image of the transformation

    ::

      ex. '/home/username/data/derived_data/DWI_transformed.mif'

-------------

Usefull links:

`mrtrix mrtransform <https://mrtrix.readthedocs.io/en/latest/reference/commands/mrtransform.html>`_

`mrtrix mrtransform - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.utils.html#mrtransform>`_
