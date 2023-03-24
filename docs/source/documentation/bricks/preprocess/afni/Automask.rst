:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Automask brick
============

Create a brain-only mask of the image using AFNI 3dAutomask command.

AFNI 3dClipLevel algorithm is used to find clipping level and after an erosion/dilation step, only the largest connected component of the supra-threshold voxels are kept.  


Be carreful this program by itself does NOT do 'skull-stripping'.

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import Automask

>>> Automask.help()

**Inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz]).
    Input dataset can be 4DEPI or a skull-stripped anatomical.

    ::

      ex. '/home/username/data/raw_data/func.nii'

- *output_type* (NIFTI or NIFTI_GZ, optional)
    | Format of the output image (one of NIFTI, NIFTI_GZ, AFNI).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *out_prefix* (a string, optional)
    Prefix of the output image. Default is 'automask_'.
    
    ::

      ex. 'autpmask_'

- *out_brain_suffix* (a string, optional)
   Suffix of the brain masked image. Default is '_masked'.
    
    ::

      ex. '_masked'

- *clfrac* (a float, optional)
    Sets the clip level fraction (must be 0.1-0.9). A small value will tend to make the mask larger. 
    Default is 0.5.
    
    ::

      ex. '0.5'

- *dilate* (an integer, optional)
    Dilate the mask outwards.
    Default is Undefined (ie parameter not used).
    
    ::

      ex. '1'  

- *erode* (an integer, optional)
    Erode the mask inwards.
    Default is Undefined (ie parameter not used).
    
    ::

      ex. '1'    


**Outputs parameters:**

- *out_file* (a strings representing a file)
    Brain mask image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/automask_func.nii'

- *out_brain* (a strings representing a file)
    Masked image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/func_masked.nii'

-------------

Usefull links:

`AFNI 3dAutomask <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dAutomask.html>`_
`AFNI 3dAutomask - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.preprocess.html#automask>`_
