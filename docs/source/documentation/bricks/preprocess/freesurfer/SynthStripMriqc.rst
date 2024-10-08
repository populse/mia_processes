:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=====================
SynthStripMriqc brick
=====================

Skull-stripping using SynthStrip tool (mri_synthstrip in FreeSurfer) as done in `mriqc v22.06 <https://mriqc.readthedocs.io/en/22.0.6/>`_
-----------------------------------------------------------------------------------------------------------------------------------------

Adapted from  `mriqc_1 <https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/synthstrip/cli.py>`_,
`mriqc_2 <https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/synthstrip/model.py>`_ and
`FreeSurfer <https://github.com/freesurfer/freesurfer/blob/2995ded957961a7f3704de57eee88eb6cc30d52d/mri_synthstrip/mri_synthstrip>`_

If you use this bricks in your analysis, please cite: `SynthStrip: Skull-Stripping for Any Brain Image. <https://doi.org/10.1016/j.neuroimage.2022.119474>`_.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input file to skull strip. An existing file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

**Optional inputs with default value parameters:**

- *border_mm* (an integer, optional, default value is 1)
    Mask border threshold in mm. Controls the boundary distance from the brain.

    ::

      ex. 1


- *model* (a string representing an existing file, optional, default value is /freeesurfer_path/models/synthstrip.1.pt)
    Alternative model weights.

    ::

      ex. '/freeesurfer_path/models/synthstrip.1.pt'

- *gpu* (a boolean, optional, default value is False)
    Exclude CSF from brain border.

    ::

      ex. False

- *output_type* ('NIFTI' or 'NIFTI_GZ' or 'MGZ', optional, default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ, MGZ).
    |   - NIFTI: \*.nii
    |   - NIFTI_GZ: \*.nii.gz
    |   - MGZ: \*.mgz

    ::

      ex. NIFTI

**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    Brain-extracted image

    ::

      ex. '/home/username/data/raw_data/T1w_desc-brain.nii'

- *out_mask* (a pathlike object or string representing a file)
    Binary brain mask

    ::

      ex. '/home/username/data/raw_data/T1w_desc-brain_mask.nii'

-------------

Useful links:

`mriqc v22.06 - cli <https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/synthstrip/cli.py>`_

`mriqc v22.06 - model <https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/synthstrip/model.py>`_

`Freesurfer SynthStrip <https://github.com/freesurfer/freesurfer/blob/2995ded957961a7f3704de57eee88eb6cc30d52d/mri_synthstrip/mri_synthstrip>`_
