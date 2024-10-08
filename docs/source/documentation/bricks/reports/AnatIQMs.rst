:orphan:

.. toctree::

+--------------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+-------------------------------------------+----------------------------------------------------+

==============
AnatIQMs brick
==============

Computes the anatomical Image Quality Metrics (IQMs) as defined in mriqc v22.06
-------------------------------------------------------------------------------

A no-reference IQM is a measurement of some aspect of the actual image which cannot be compared to a reference value for
the metric since there is no ground-truth about what this number should be.
All the computed IQMs corresponding to an image are saved in a JSON file.

Adapted from `mriqc <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/anatomical.py#L332>`_.

--------------------------------------

**Mandatory inputs parameters:**

- *in_ras* (a string representing an existing file)
    Input anatomical image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

**Optional inputs with default value parameters:**

- *airmask* (a string representing an existing file)
    Air mask (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/air_T1w_mask.nii'

- *artmask* (a string representing an existing file)
    Artifact mask (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/art_T1w_mask.nii'

- *hatmask* (a string representing an existing file)
    Hat mask (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/hat_T1w_mask.nii'

- *headmask* (a string representing an existing file)
    Head mask (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/n4c_T1w_clipped_brain_outskin_mask.nii'

- *in_fwhm* (a string representing an existing file)
    A file with FWHM information.

    ::

      ex. '/home/username/data/derived_data/fwhm_T1w.out'

- *in_inu* (a string representing an existing file)
    Bias image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/n4c_T1w_clipped_bias.nii'

- *in_noinu* (a string representing an existing file)
    Harmonized image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/n4c_T1w_clipped_harmonized.nii'

- *mni_tpms* (a list with items which are string representing an existing file)
    Templates registered in the subject space (valid extensions: [.nii, .nii.gz]).

    ::

      ex. ['/home/username/data/derived_data/csf_tpl-MNI152NLin2009cAsym_res-01_label-CSF_probseg.nii.gz', '/home/username/data/derived_data/gm_tpl-MNI152NLin2009cAsym_res-01_label-GM_probseg.nii.gz', '/home/username/data/derived_data/wm_tpl-MNI152NLin2009cAsym_res-01_label-WM_probseg.nii.gz']

- *pvms* (a list with items which are string representing an existing file)
    Tissues probality map (valid extensions: [.nii, .nii.gz]).

    ::

      ex. ['/home/username/data/derived_data/ss_n4c_T1w_clippep_pve_0.nii', '/home/username/data/derived_data/ss_n4c_T1w_clippep_pve_1.nii', '/home/username/data/derived_data/ss_n4c_T1w_clippep_pve_2.nii'

- *rotmask* (a string representing an existing file)
    Rotation mask (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/T1w_rotmasked.nii'

- *segmentation* (a string representing an existing file)
    Segmentation mask (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/ss_n4c_T1w_clippep_seg.nii'

**Outputs parameters:**

- *out_file* (a strings representing a file)
    A json file with all the IQMs computed.

    ::

      ex. '/home/username/data/derived_data/T1w_anat_qc.json'

-------------

Useful links:

`mriqc Anat qc  <https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/qc/anatomical.py>`_
