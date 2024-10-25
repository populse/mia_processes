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

- *bat_window_size*
    Number of time points (or dynamics) used to determine the bolus arrival
    time (t0) in the MRI signal analysis (an integer). Acts as a sliding
    window that moves in the time dimension of the MRI data, allowing the
    algorithm to analyse the local temporal behaviour at each voxel. A large
    window provides smoother estimates of the signal, but can delay the
    detection of sudden changes, such as the arrival of the bolus. It
    increases stability, but can also *blur* the sharp transition at bolus
    arrival. A smaller window size will be more sensitive to rapid changes,
    but may also pick up more noise or fluctuations in the signal. The default
    value (8) is a moderately sized window to smooth out short-term
    fluctuations while detecting significant changes in the MRI signal value.

    ::

      ex. 8

- *bat_th*
    Threshold multiplier (a float) is used to detect significant changes in
    the signal intensity of MRI data (controls the sensitivity of the bolus
    arrival detection). A high *th* value (e.g., th = 3.0) would make the
    algorithm less sensitive, meaning it will only consider larger, more
    pronounced signal drops as valid bolus arrival points. This reduces false
    positives but may miss subtle or gradual bolus arrivals. A low th value
    (e.g., th = 1.0) would make the algorithm more sensitive, catching even
    smaller signal drops. However, this may lead to false positives, where
    noise or natural fluctuations in the signal are incorrectly identified as
    bolus arrival. The default value (2.0) ensures only signal drops larger
    than 2 standard deviations below the mean are considered significant,
    providing a balance between sensitivity and robustness against noise.

    ::

      ex. 2.0

- *wmaxr*
    Scaling factor (a float) used to define the threshold for acceptable peak
    widths of voxel intensity curves. It is a percentage of the total number
    of dynamics, which represents the maximum width allowed for peaks within
    the time series data of each voxel. Increasing *wmaxr* will relax the
    criterion for peak width, allowing wider peaks to be included in the
    region of interest (ROI). This might allow more voxels with broader peaks
    but could also include noisier or less precise signals. Decreasing *wmaxr*
    will make the criterion stricter, permitting only narrower peaks and
    potentially excluding more noisy voxels. This could increase precision but
    might also exclude valid data.

    ::

      ex. 0.5

- *nb_vox_cand*
    Number of candidate voxels to be evaluated for computing the AIF (an
    integer).

    ::

      ex. 50

- *nb_vox_best_scor*
    Number of best-scoring voxels to be used for computing the AIF after
    evaluating the initial set of candidate voxels (an integer).

    ::

      ex. 5

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
