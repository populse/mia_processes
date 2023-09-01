:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

====================
FilteringTrack brick
====================

Filter tracking data
--------------------

Filter a whole-brain fibre-tracking data set such that the streamline densities match the FOD lobe integrals

(mrtrix tcksift command)

--------------------------------------

**Mandatory inputs parameters:**

- *in_tracks* (a pathlike object a string representing an existing file)
    Inputs track file (valid extensions: [.tck]).

    ::

      ex. '/home/username/data/raw_data/tracks.tck'

- *in_fod* (a pathlike object a string representing an existing file)
    Input image containing the spherical harmonics of the fibre orientation distributions

    ::

      ex. '/home/username/data/raw_data/wm_fod_norm.mif'

**Optional inputs with default value parameters:**

- *suffix* (a string, default value is sift, optional)
    Output file suffix

    ::

      ex. 'sift'


- *fd_scale_gm* (a boolean, default value is False, optional)
    Provide this option (in conjunction with -act) to heuristically downsize the fibre density
    estimates based on the presence of GM in the vox

    ::

      ex. False

- *no_dilate_lut* (a boolean, default value is False, optional)
    Do NOT dilate FOD lobe lookup tables; only map streamlines to FOD lobes if the precise tangent
    lies within the angular spread of that lobe

    ::

      ex. False

- *make_null_lobes* (a boolean, default value is False, optional)
    Add an additional FOD lobe to each voxel, with zero integral,
    that covers all directions with zero / negative FOD amplitudes

    ::

      ex. False

- *remove_untracked* (a boolean, default value is False, optional)
    Remove FOD lobes that do not have any streamline density attributed to them

    ::

      ex. False

- *get_csv_file* (a boolean, default value is False, optional)
    Output statistics of execution per iteration to a .csv file

    ::

      ex. False

- *get_mu_file* (a boolean, default value is False, optional)
    Output statistics of execution per iteration to a .csv file

    ::

      ex. False

- *get_out_selection_file* (a boolean, default value is False, optional)
    Output a text file containing the binary selection of streamlines

    ::

      ex. False

**Optional inputs parameters:**

- *proc_mask* (string representing an existing file, optional)
    Provide an image containing the processing mask weights for the model.

    ::

      ex. '/home/username/data/derived_data/mask.mif'

- *act_image* (string representing an existing file, optional)
    Provide an ACT five-tissue-type segmented anatomical image to derive the processing mask

    ::

      ex. '/home/username/data/derived_data/5tt_coreg.mif'

- *fd_thresh_value* (an integer, optional)
    Fibre density threshold

    ::

      ex.

- *term_ratio_value* (an integer, optional)
    Termination ratio - defined as the ratio between reduction in
    cost function, and reduction in density of streamlines.

    ::

      ex.

- *term_mu_value* (an integer, optional)
    Terminate filtering once the SIFT proportionality coefficient reaches a given value

    ::

      ex.


**Outputs parameters:**

- *tracks_out* (a pathlike object or string representing a file)
    The output filtered track file

    ::

      ex. '/home/username/data/derived_data/tracks_sift.tck'

- *csv_file_out* (a pathlike object or string representing a file, optional)
    A csv file with the output statistics of execution per iteration

    ::

      ex. '/home/username/data/derived_data/tracks_tcksift_stats.csv

- *mu_file_out* (a pathlike object or string representing a file, optional)
    The final value of SIFT proportionality coefficient mu in a text file

    ::

      ex. '/home/username/data/derived_data/tracks_tcksift_mu.txt'

- *selection_file_out* (a pathlike object or string representing a file, optional)
    A text file containing the binary selection of streamlines

    ::

      ex. '/home/username/data/derived_data/tracks_tcksift_selection.txt'


-------------

Usefull links:

`mrtrix tcksift <https://mrtrix.readthedocs.io/en/latest/reference/commands/tcksift.html>`_
