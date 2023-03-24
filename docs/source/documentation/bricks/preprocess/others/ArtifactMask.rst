:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
ArtifactMask brick
============

Computes the artifact mask using the method described in the step 1 (Background Region Segmentation) of:
Mortamet B et al.
Automatic quality assessment in structural brain magnetic resonance imaging,
Mag Res Med 62(2):365-372, 2009. doi:10.1002/mrm.21992.

Adapted from https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/anatomical.py#L301

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import ArtifactMask
>>> ArtifactMask.help()

**Inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *head_mask* (a string representing an existing file)
    Head mask image (valid extensions: [.nii, .nii.gz])
    
    ::

      ex. '/home/username/data/raw_data/T1w_outskin_mask.nii'

- *rot_mask* (a string representing an existing file, optional)
    Rotation mask image (valid extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/raw_data/T1w_rotmasked.nii'

- *nasion_post_mask* (a string representing an existing file)
    Template mask registred in the in-file space (valid extensions: [.nii, .nii.gz]). 
    
    
    ::

      ex. '/home/username/data/raw_data/t_tpl-MNI152NLin2009cAsym_res-01_desc-head_mask.nii'

- *suffix* (a string, optional)
   Suffix of output images. 
    
    ::

      ex. '_mask'


**Outputs parameters:**

- *out_art_mask* (a pathlike object or string representing a file)
    Out artifact mask image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/art_T1w_mask.nii'

- *out_hat_mask* (a pathlike object or string representing a file)
    Out hat mask image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/hat_T1w_mask.nii'

- *out_air_mask* (a pathlike object or string representing a file)
    Out air mask ("hat" mask without artifacts) image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/air_T1w_mask.nii'


-------------

Usefull links:
`Artifactmask mriqc - nipype <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/anatomical.py#L301>`_
