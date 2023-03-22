:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
SynthStrip brick
============
 
 Skull-stripping using SynthStrip tool (mri_synthstrip in Freesurfer). 
 It is a is a skull-stripping tool that extracts brain signal from a landscape of image types, ranging across imaging modality, contrast, resolution, and subject population. 
 It leverages a deep learning strategy  that synthesizes arbitrary training images from segmentation maps to optimize a robust model agnostic to acquisition specifics.


>>> from mia_processes.bricks.preprocess.freesurfer import SynthStrip
>>> SynthStrip.help()

**Inputs parameters:**

- *in_file* (a string representing an existing file)
    Input file to skull strip. An existing file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. /home/username/data/raw_data/T1w.nii


- *border_mm* (integer)
    Mask border threshold in mm. Controls the boundary distance from the brain.
    Default is 1

    ::

      ex. 1

- *no_csf* (boolean)
    Exclude CSF from brain border. Default is Fasle

    ::

      ex. False

- *output_type* ('NIFTI' or 'NIFTI_GZ' or 'MGZ')
    | Format of the output image (one of NIFTI, NIFTI_GZ, MGZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz
    |   MGZ: \*.mgz

    ::

      ex. NIFTI


**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    Brain-extracted image

    ::

      ex. /home/username/data/raw_data/T1w_desc-brain.nii

- *out_mask* (a pathlike object or string representing a file)
    Binary brain mask 

    ::

      ex. /home/username/data/raw_data/T1w_desc-brain_mask.nii

-------------

Usefull links:
`Freesurfer SynthStrip <https://surfer.nmr.mgh.harvard.edu/docs/synthstrip/>`_