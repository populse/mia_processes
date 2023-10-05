:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==================
Tractography brick
==================

Performs streamlines tractography
---------------------------------

Perform streamlines tractography  using tckgen mrtrix command.

At least one tractography seeding mechanisms must be provided. The select and seed options allow to choose the number of desired streamlines.

(mrtrix tckgen command)

*Please note that, in mia_processes, MRtrix bricks required diffusion data in MRtrix .mif format.
The MRConvert brick could be used to convert diffusion data in NIfTI format into MRtrix .mif format.
In this case, bvec and bval file should be specified.*

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a pathlike object or a string representing an existing file)
    Input DWI image (valid extensions: [.mif]).

    ::

      ex. '/home/username/data/derived_data/DWI.mif'


**Optional inputs with default value parameters:**

- *algorithm* (iFOD2, FACT, iFOD1, Nulldist, SD_Stream, Tensor_Det or Tensor_Prob default value is iFOD2, optional)
    | Tractography algorithm to be used:
    |   - iFOD2: Second-order Integration over Fiber Orientation Distributions
    |   - iFOD1: First-order Integration over Fiber Orientation Distributions
    |   - FACT: Fiber Assigned by Continuous Tracking
    |   - NullDist: Null Distribution tracking algorithms
    |   - SD_Stream: Streamlines tractography based on Spherical Deconvolution (SD)
    |   - Tensor_Det: A deterministic algorithm that takes as input a 4D diffusion-weighted image (DWI) series.
    |   - Tensor_Prob: A probabilistic algorithm that takes as input a 4D diffusion-weighted image (DWI) series

    ::

      ex. iFOD2

- *cutoff* (a float, default value is 0.1, optional)
    Streamlines tractography options - The FA or FOD amplitude for terminating tracks

    ::

      ex. 0.1

- *nopreselectcompt* (a boolean, default value is False, optional)
    Streamlines tractography options - Do NOT pre-compute legendre polynomial values

    ::

      ex. False

- *stop* (a boolean, default value is False, optional)
    Streamlines tractography options - Stop propagating a streamline once it has traversed all include regions

    ::

      ex. False

- *use_rk4* (a boolean, default value is False, optional)
    Streamlines tractography options - Use 4th-order Runge-Kutta integration

    ::

      ex. False

- *tracto_seed_unidirectional* (a boolean, default value is False, optional)
    Tractography seeding options and parameters - Track from the seed point in one direction only

    ::

      ex. False

- *tracto_get_output_seeds* (a boolean, default value is False, optional)
    Tractography seeding options and parameters - Get the seed location of all successful streamlines into a file

    ::

      ex. False

- *backtrack* (a boolean, default value is False, optional)
    Anatomically-Constrained Tractography options - Allow tracks to be truncated.

    ::

      ex. False

- *crop_at_gmwmi* (a boolean, default value is False, optional)
    Anatomically-Constrained Tractography options - Crop streamline endpoints more precisely as they cross the GM-WM interface

    ::

      ex. False

**Optional inputs:**

- *angle* (a float, optional)
    Streamlines tractography options - Maximum angle between successive steps.

    If not set, default is 60 for deterministic algorithms, 15 for iFOD1 / nulldist1 and 45 for iFOD2 / nulldist2

    ::

      ex. 45

- *downsample_factor* (a float, optional)
    Streamlines tractography options - Downsample the generated streamlines to reduce output file size

    If not set, default is (samples-1) for iFOD2, no downsampling for all other algorithms.

    ::

      ex. 5.0

- *max_lenght* (a float, optional)
    Streamlines tractography options - The max length of any track in mm

    If not set, default value is 100 x voxelsize

    ::

      ex. 200.0


- *min_lenght* (a float, optional)
    Streamlines tractography options - The minimum length of any track in mm

    If not set, default value is without ACT 5 x voxelsize and with ACT 2 x voxelsize

    ::

      ex. 4

- *step_size* (a float, optional)
    | Streamlines tractography options - Step size of the algorithm in mm

    | If not set, default value is:
    | - for first-order algorithms, 0.1 x voxelsize;
    | - if using RK4, 0.25 x voxelsize
    | - for iFOD2: 0.5 x voxelsize

    ::

      ex. 0.32

- *trials* (a integer, optional)
    Streamlines tractography options - Maximum number of sampling trials at each point (only used for probabilistic tracking)
    If not set, default value is 1000.

    ::

      ex. 1000

- *select* (an integer, optional)
    Streamlines tractography options - Desired number of tracks after all selection criteria have been applied (i.e. inclusion/exclusion ROIs, min/max length, etc).
    If this option is not specified, by default, 5000 streamlines are to be selected.

    ::

      ex. 5000


- *seed_dynamic* (a string representing an existing file, optional)
    Tractography seeding mechanisms - Determine seed points dynamically using the SIFT model (must not provide any other seeding mechanism).
    Note that while this seeding mechanism improves the distribution of reconstructed streamlines density,
    it should NOT be used as a substitute for the SIFT method itself.

    ::

      ex. '/home/username/data/derived_data/fod.mif'

- *seed_gmwmi* (a string representing an existing file, optional)
    Tractography seeding mechanisms - Seed from the grey matter - white matter interface (only valid if using ACT framework).

    ::

      ex. '/home/username/data/derived_data/gmwmSeed_coreg.mif'

- *seed_grid_voxel* (a tuple of the form: (a pathlike object or string representing an existing file, an integer), optional)
    Tractography seeding mechanisms - Seed a fixed number of streamlines per voxel in a mask image;
    place seeds on a 3D mesh grid (grid_size argument is per axis; so a grid_size of 3 results in 27 seeds per voxel)

    ::

      ex. ('/home/username/data/derived_data/mesh_grid.mif', 5)

- *seed_image* (a string representing an existing file, optional)
    Tractography seeding mechanisms - Seed streamlines entirely at random within a mask image.

    ::

      ex. '/home/username/data/derived_data/gmwmSeed_coreg.mif'

- *seed_rejection* (a string representing an existing file, optional)
    Tractography seeding mechanisms - Seed from an image using rejection sampling

    ::

      ex. '/home/username/data/derived_data/gmwmSeed_coreg.mif'

- *seed_rnd_voxel* (a tuple of the form: (a pathlike object or string representing an existing file, an integer), optional)
    Tractography seeding mechanisms - Seed a fixed number of streamlines per voxel in a mask image; random placement of seeds in each voxel.

    ::

      ex. ()

- *seed_sphere* (a tuple of the form: (a float, a float, a float, a float), optional)
    Tractography seeding mechanisms - Spherical seed coordinate XYZ position and radius (X, Y, Z, r)

    ::

      ex. (55.6, 95.3, 25.3, 0.9)

- *tracto_seeds_number* (a integer, optional)
    Tractography seeding options and parameters - Number of seeds that the algorithm will attempt to track from.
    If this option is NOT provided, the default number of seeds is set to 1000 × the number of selected streamlines.

    If select option is NOT also specified: the tracking will continue  until this number of seeds has been attempted.
    If select option is also specified, the tracking will stop when the number of seeds attempted reaches the number specified here,
    OR when the number of streamlines selected reaches the number requested with the select option.

    If tracto_seeds_number option is set to zero, the algorithm will keep attempting seeds until the number specified by select option has been reached.

    ::

      ex. 2000

- *tracto_max_attempts_per_seed_number* (a integer, optional)
    Tractography seeding options and parameters - Set the maximum number of times that the tracking algorithm should attempt to find an appropriate tracking direction from a given seed point.
    This should be set high enough to ensure that an actual plausible seed point is not discarded prematurely as being unable to initiate tracking from.

    If tracto_max_attempts_per_seed_number is not set, the maximum attempts is 1000.

    ::

      ex. 2000

- *tracto_seed_cutoff* (a float, optional)
    Tractography seeding options and parameters - The minimum FA or FOD amplitude for seeding tracks.

    If tracto_seed_cutoff is not set, the cut off is set to the normal cut off (see cutoff option)

    ::

      ex. 0.2

- *tracto_seed_direction* (a tuple of the form: (a float, a float, a float), optional)
    Tractography seeding options and parameters - Specify a seeding direction for the tracking.

    ::

      ex. (0.2, 0.9, 1.3)

- *roi_excl* (string representing an existing file or a tuple of the form: (a float, a float, a float, a float), optional)
    Region Of Interest processing options - Specify an exclusion region of interest, streamlines that enter ANY exclude region will be discarded.

    ::

      ex. '/home/username/data/derived_data/roi.mif'

- *roi_incl* (string representing an existing file or a tuple of the form: (a float, a float, a float, a float), optional)
    Region Of Interest processing options - Specify an inclusion region of interest, streamlines must traverse ALL inclusion regions to be accepted.

    ::

      ex. '/home/username/data/derived_data/roi.mif'

- *roi_incl_ordered* (string representing an existing file or a tuple of the form: (a float, a float, a float, a float), optional)
    Region Of Interest processing options - Specify an inclusion region of interest, streamlines must traverse ALL inclusion_ordered regions in the order
    they are specified in order to be accepted.

    ::

      ex. '/home/username/data/derived_data/roi.mif'

- *roi_mask* (string representing an existing file or a tuple of the form: (a float, a float, a float, a float), optional)
    Region Of Interest processing options - Specify a masking region of interest. If defined,streamlines exiting the mask will be truncated.

    ::

      ex. '/home/username/data/derived_data/mask.mif'

- *act_image* (string representing an existing file, optional)
    Anatomically-Constrained Tractography options - Use the Anatomically-Constrained Tractography framework during tracking.
    The provided image must be in the 5TT ie five tissue type format.

    ::

      ex. '/home/username/data/derived_data/T1w_5tt.mif'

- *iFOD_power* (a integer, optional)
    Options specific to the iFOD tracking algorithms - Raise the FOD to the power specified (default is 1/nsamples)

    ::

      ex. 2

- *iFOD2_n_samples* (a integer, optional)
    Options specific to the iFOD2 tracking algorithm - Number of FOD samples to take per step for the 2nd order (iFOD2) method

    ::

      ex. 2

**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    Output file containing tracks

    ::

      ex. '/home/username/data/derived_data/DWI_tracto.tck'

- *out_file* (a pathlike object or string representing a file)
    Out seed location of all successful streamlines

    ::

      ex. '/home/username/data/derived_data/DWI__tracto_out_seeds.mif'


-------------

Usefull links:

`mrtrix tckgen <https://mrtrix.readthedocs.io/en/latest/reference/commands/tckgen.html>`_

`mrtrix tckgen - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.tracking.html#tractography>`_
