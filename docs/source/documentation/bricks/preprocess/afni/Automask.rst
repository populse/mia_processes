:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==============
Automask brick
==============

Create a brain-only mask of the image using AFNI 3dAutomask
-----------------------------------------------------------

AFNI 3dClipLevel algorithm is used to find clipping level and after an erosion/dilation step, only the largest connected component of the supra-threshold voxels are kept.

Be carreful this program by itself does NOT do 'skull-stripping'.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz]).
    Input dataset can be 4DEPI or a skull-stripped anatomical.

    ::

      ex. '/home/username/data/raw_data/func.nii'

**Optional inputs with default value parameters:**

- *clfrac* (a float, optional, default value is 0.5)
    Sets the clip level fraction (must be 0.1-0.9). A small value will tend to make the mask larger.

    ::

      ex. 0.5

- *out_brain_suffix* (a string, optional, default value is '_masked')
    Suffix of the brain masked image.

    ::

      ex. '_masked'

- *output_type* (NIFTI or NIFTI_GZ, optional, default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *out_prefix* (a string, optional, default value is 'automask')
    Prefix of the output image.

    ::

      ex. 'automask_'

**Optional inputs parameters:**

- *dilate* (an integer, optional)
    Dilate the mask outwards. Default is Undefined (ie parameter not used).

    ::

      ex. 1

- *erode* (an integer, optional)
    Erode the mask inwards. Default is Undefined (ie parameter not used).

    ::

      ex. 1


**Outputs parameters:**

- *out_brain* (a strings representing a file)
    Masked image (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/func_masked.nii'

- *out_file* (a strings representing a file)
    Brain mask image (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/automask_func.nii'

-------------

Usefull links:

`AFNI 3dAutomask <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dAutomask.html>`_

`AFNI 3dAutomask - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.preprocess.html#automask>`_
