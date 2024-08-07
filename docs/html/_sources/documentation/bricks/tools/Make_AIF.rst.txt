:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

==============
Make_AIF brick
==============

Compute the Arterial Input Function (AIF) from dynamic MRI perfusion data
-------------------------------------------------------------------------

- **Loading Data**: The code begins by loading 4D MRI data (3D spatial data
  over time) and determining its dimensions.

- **Thresholding and Masking**: It applies an intensity threshold to create a
  mask that identifies relevant brain tissue, excluding non-brain areas and
  noisy voxels.

- **Filtering Candidates**: The algorithm selects candidate voxels based on
  their dynamic signal characteristics, such as intensity changes and peak
  widths.

- **Peak Analysis**: It calculates key metrics like the baseline intensity,
  minimum intensity, and the width of the dynamic peak for each voxel.

- **Excluding Artifacts**: Voxels with artifacts or unrealistic dynamic
  behaviors are excluded to ensure accuracy.

- **Scoring and Selection**: Each voxel is scored based on its dynamic
  properties. The top candidate voxels are then selected, and their signals
  are averaged to compute the final AIF.

- **Saving Results**: The resulting AIF and the scoring details of the
  selected voxels are saved in a JSON file for further analysis.

-----------------------------------------------

**Inputs parameters:**

- *func_file*
    T2* functional Magnetic Resonance Imaging (fMRI) experiment recorded
    during gadolinium bolus. Ideally, the data will have been pre-processed
    (realignment, segmentation, etc.). An existing, uncompressed file
    (valid extensions: [.img, .nii, .hdr]).

    ::

      ex. '/home/username/data/raw_data/swrfunc.nii'

**Outputs parameters:**

- *aif_file*
    The output data from the DSC-MRI AIF computation (a file with
    .json format). It includes the computed Arterial Input Function
    (`aif` key) and associated scoring details of selected voxels
    (`scores` key).

    In `scores`, the elements are, in this order:

      - Score value: The more intense and finer the peak, the higher the score.
      - Row index: Row index of the selected voxel.
      - Column index: Column index of the selected voxel.
      - Slice index: slice index of the selected voxel.
      - Number of warnings : The number of warnings when calculating AIF.
      - Pre-bolus baseline test: Indicates whether the pre-bolus baseline is
        too noisy. `None` value indicates no problem detected.
      - Post-bolus baseline test: Indicates whether the post-bolus baseline
        is too noisy. `None` value indicates no problem detected.
      - t0 point test: Indicates whether the voxel value at the time of bolus
        arrival (t0) is greater than 11/10th of the baseline mean. `None`
        value indicates no problem detected.
      - Pre-bolus baseline length test: Indicates whether the pre-bolus
        baseline length is too short (< 8 dynamics). `None` value indicates
        no problem detected.

    ::

      ex. '/home/username/data/derived_data/swrfunc_aif.json'
