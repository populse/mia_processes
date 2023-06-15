:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==================
Level1Design brick
==================

fMRI model specification for GLM analysis.
------------------------------------------

**Inputs parameters:**

- *timing_units <=> timing.units* [#label]_
    Units for specification of events or blocks onsets. One of 'scans' or 'secs'.

    ::

      ex. secs

- *interscan_interval <=> timing.RT* [#label]_
    The time between acquiring a plane of one volume and the same plane in the next volume (in secs, a float). If this parameter is undefined at the initialisation time, its value is retrieved from the database

    ::

      ex. 3.0

- *microtime_resolution <=> timing.fmri_t* [#label]_
    The number of time-bins per scan used to build regressors (an integer). Typically, no need to change this parameter, except if a slice-timing correction were performed previously, in order to match the number of slices specified there.

    ::

      ex. 16

- *microtime_onset <=> timing.fmri_t0* [#label]_
    The reference time-bin at which the regressors are resampled to coincide with data acquisition (an integer). If a slice-timing correction has been made, this parameter will have to be modified to match the reference slice specified in it.

    ::

      ex. 8

- *sess_scans <=> sess.scans* [#label]_
    The fMRI scans for the sessions (a list of items which are a pathlike object or string representing an existing file). They must all have the same image dimensions, orientation, voxel size, etc. Each scan will define (will be associated with) one session. In a general way, ``sess_scans`` can be defined as:

    .. code-block:: python

        [sess_1_ScanPath_1, sess_2_scanPath_2, ..., sess_n_scanPath_n]

        ex. ['/home/username/data/raw_data/Func_1.nii', '/home/username/data/raw_data/Func_2.nii']

- *sess_cond_names <=> sess.cond.name* [#label]_
    The name of each condition (list of items which are a list of items which are a string). Each session can have several or no conditions, so any number of condition types (event or time) can be specified for each session. In a general way, ``sess_cond_names`` can be defined as:

    .. code-block:: python

        [
          [sess_1_condName_1, sess_1_condName_2, ..., sess_1_condName_n],
          [sess_2_condName_1, sess_2_condName_2, ..., sess_2_condName_n],
          ...,
          [sess_n_condName_1, sess_n_condName_2, ..., sess_n_condName_n]
        ] # Use None if there is no condition for a session

        ex. [None, ['cond1', 'cond2']]

- *sess_cond_onsets <=> sess.cond.onsets* [#label]_
    The onset times (in seconds or in scans) of the epochs or events within each condition (a list of items which are a list of items which are a list of items which are a float). In a general way, ``sess_cond_onsets`` can be defined as:

    .. code-block:: python

        [
          [
            [sess_1_cond_1_onset_1, sess_1_cond_1_onset_2, ..., sess_1_cond_1_onset_n],
	    [sess_1_cond_2_onset_1, sess_1_cond_2_onset_2, ..., sess_1_cond_2_onset_n],
	    ...,
	    [sess_1_cond_n_onset_1, sess_1_cond_n_onset_2, ..., sess_1_cond_n_onset_n]
          ],
          [
            [sess_2_cond_1_onset_1, sess_2_cond_1_onset_2, ..., sess_2_cond_1_onset_n],
	    [sess_2_cond_2_onset_1, sess_2_cond_2_onset_2, ..., sess_2_cond_2_onset_n],
	    ...,
	    [sess_2_cond_n_onset_1, sess_2_cond_n_onset_2, ..., sess_2_cond_n_onset_n]
          ],
          ...,
          [
            [sess_n_cond_1_onset_1, sess_n_cond_1_onset_2, ..., sess_n_cond_1_onset_n],
	    [sess_n_cond_2_onset_1, sess_n_cond_2_onset_2, ..., sess_n_cond_2_onset_n],
	    ...,
	    [sess_n_cond_n_onset_1, sess_n_cond_n_onset_2, ..., sess_n_cond_n_onset_n]
          ]
        ] # Use None if there is no condition for a session

        ex. [None, [[43.72, 117.27, 152.59, 225.20], [21.23, 98.21, 206.98, 282.41]]]

- *sess_cond_durations <=> sess.cond.durations* [#label]_
    The duration times (in seconds or in scans) of the epochs within each condition (a list of items which are a list of items which are a list of items which are a float). Events are specified with a duration of 0. If a single number is specified for durations, it will be assumed that all trials comply with that duration. In a general way, ``sess_cond_durations`` can be defined as:

    .. code-block:: python

        [
          [
            [sess_1_cond_1_duration_1, sess_1_cond_1_duration_2, ..., sess_1_cond_1_duration_n],
	    [sess_1_cond_2_duration_1, sess_1_cond_2_duration_2, ..., sess_1_cond_2_duration_n],
	    ...,
	    [sess_1_cond_n_duration_1, sess_1_cond_n_duration_2, ..., sess_1_cond_n_duration_n]
          ],
          [
            [sess_2_cond_1_duration_1, sess_2_cond_1_duration_2, ..., sess_2_cond_1_duration_n],
	    [sess_2_cond_2_duration_1, sess_2_cond_2_duration_2, ..., sess_2_cond_2_duration_n],
	    ...,
	    [sess_2_cond_n_duration_1, sess_2_cond_n_duration_2, ..., sess_2_cond_n_duration_n]
          ],
          ...,
          [
            [sess_n_cond_1_duration_1, sess_n_cond_1_duration_2, ..., sess_n_cond_1_duration_n],
	    [sess_n_cond_2_duration_1, sess_n_cond_2_duration_2, ..., sess_n_cond_2_duration_n ],
	    ...,
	    [sess_n_cond_n_duration_1, sess_n_cond_n_duration_2, ..., sess_n_cond_n_duration_n]
          ]
        ] # Use None if there is no condition for a session

        ex. [None, [[0], [10.54, 10.42, 8.59, 9.81]]]

- *sess_cond_tmod <=> sess.cond.tmod* [#label]_
    Allows for the characterisation of linear or nonlinear time effects (a list of items which are a list of items which are 0 or 1 or 2 or 3 or 4 or 5 or 6). One time modulation parameter for each condition must be applied.
      | - 0: No time modulation
      | - 1: 1st order time modulation
      | - 2: 2nd order time modulation
      | - 3: 3rd order time modulation
      | - 4: 4th order time modulation
      | - 5: 5th order time modulation
      | - 6: 6th order time modulation

    In a general way, ``sess_cond_tmod`` can be defined as:

    .. code-block:: python

        [
          [sess_1_cond_1_tmod, sess_1_cond_2_tmod, ..., sess_1_cond_n_tmod],
          [sess_2_cond_1_tmod, sess_2_cond_2_tmod, ..., sess2_cond_n_tmod],
          ...,
          [sess_n_cond_1_tmod, sess_n_cond_2_tmod, ..., sess_n_cond_n_tmod]
        ] # Use None if there is no condition for a session

        ex. [None, [3, 0]]

- *sess_cond_pmod_names <=> sess.cond.pmod.name* [#label]_
    The name of the parametric modulation (a list of items which are a list of items which are a list of items which are a string). Modeling interactions with these user specified parameters, allows model in the design matrix nonlinear effects relating to some other measures. Several (or none) parametric modulations can be applied for each condition. In a general way, ``sess_cond_pmod_names`` can be defined as:

    .. code-block:: python

        [
          [
            [sess_1_cond_1_pmodName_1, sess_1_cond_1_pmodName_2, ..., sess_1_cond_1_pmodName_n],
	    [sess_1_cond_2_pmodName_1, sess_1_cond_2_pmodName_2, ..., sess_1_cond_2_pmodName_n],
	    ...,
	    [sess_1_cond_n_pmodName_1, sess_1_cond_n_pmodName_2, ..., sess_1_cond_n_pmodName_n]
          ],
          [
            [sess_2_cond_1_pmodName_1, sess_2_cond_1_pmodName_2, ..., sess_2_cond_1_pmodName_n],
	    [sess_2_cond_2_pmodName_1, sess_2_cond_2_pmodName_2, ..., sess_2_cond_2_pmodName_n],
	    ...,
	    [sess_2_cond_n_pmodName_1, sess_2_cond_n_pmodName_2, ..., sess_2_cond_n_pmodName_n]
          ],
          ...,
          [
            [sess_n_cond_1_pmodName_1, sess_n_cond_1_pmodName_2, ..., sess_n_cond_1_pmodName_n],
	    [sess_n_cond_2_pmodName_1, sess_n_cond_2_pmodName_2, ..., sess_n_cond_2_pmodName_n],
	    ...,
	    [sess_n_cond_n_pmodName_1, sess_n_cond_n_pmodName_2, ..., sess_n_cond_n_pmodName_n]
          ]
        ] # Use None if there is no parametric modulation for a condition or if there is no condition for a session

        ex. [None, [None, ['pmod1', 'pmod2']]]

- *sess_cond_pmod_values <=> sess.cond.pmod.param* [#label]_
    The values used for the parametric modulation, one for each occurence of the event (a list of items which are a list of items which are a list of items which are a list of items which are a float). In a general way, ``sess_cond_pmod_values`` can be defined as:

    .. code-block:: python

        [
          [
            [
              [sess_1_cond_1_pmod_1_val_1, sess_1_cond_1_pmod_1_val_2, ..., sess_1_cond_1_pmod_1_val_n],
              [sess_1_cond_1_pmod_2_val_1, sess_1_cond_1_pmod_2_val_2, ..., sess_1_cond_1_pmod_2_val_n],
              ...,
              [sess_1_cond_1_pmod_n_val_1, sess_1_cond_1_pmod_n_val_2, ..., sess_1_cond_1_pmod_n_val_n]
            ],
            [
              [sess_1_cond_2_pmod_1_val_1, sess_1_cond_2_pmod_1_val_2, ..., sess_1_cond_2_pmod_1_val_n],
              [sess_1_cond_2_pmod_2_val_1, sess_1_cond_2_pmod_2_val_2, ..., sess_1_cond_2_pmod_2_val_n],
              ...,
              [sess_1_cond_2_pmod_n_val_1, sess_1_cond_2_pmod_n_val_2, ..., sess_1_cond_2_pmod_n_val_n]
            ],
            ...,
            [
              [sess_1_cond_n_pmod_1_val_1, sess_1_cond_n_pmod_1_val_2, ..., sess_1_cond_n_pmod_1_val_n],
              [sess_1_cond_n_pmod_2_val_1, sess_1_cond_n_pmod_2_val_2, ..., sess_1_cond_n_pmod_2_val_n],
              ...,
              [sess_1_cond_n_pmod_n_val_1, sess_1_cond_n_pmod_n_val_2, ..., sess_1_cond_n_pmod_n_val_n]
            ]
          ],
          [
            [
              [sess_2_cond_1_pmod_1_val_1, sess_2_cond_1_pmod_1_val_2, ..., sess_2_cond_1_pmod_1_val_n],
              [sess_2_cond_1_pmod_2_val_1, sess_2_cond_1_pmod_2_val_2, ..., sess_2_cond_1_pmod_2_val_n],
              ...,
              [sess_2_cond_1_pmod_n_val_1, sess_2_cond_1_pmod_n_val_2, ..., sess_2_cond_1_pmod_n_val_n]
            ],
            [
              [sess_2_cond_2_pmod_1_val_1, sess_2_cond_2_pmod_1_val_2, ..., sess_2_cond_2_pmod_1_val_n],
              [sess_2_cond_2_pmod_2_val_1, sess_2_cond_2_pmod_2_val_2, ..., sess_2_cond_2_pmod_2_val_n],
              ...,
              [sess_2_cond_2_pmod_n_val_1, sess_2_cond_2_pmod_n_val_2, ..., sess_2_cond_2_pmod_n_val_n]
            ],
            ...,
            [
              [sess_2_cond_n_pmod_1_val_1, sess_2_cond_n_pmod_1_val_2, ..., sess_2_cond_n_pmod_1_val_n],
              [sess_2_cond_n_pmod_2_val_1, sess_2_cond_n_pmod_2_val_2, ..., sess_2_cond_n_pmod_2_val_n],
              ...,
              [sess_2_cond_n_pmod_n_val_1, sess_2_cond_n_pmod_n_val_2, ..., sess_2_cond_n_pmod_n_val_n]
            ]
          ],
          ...,
          [
            [
              [sess_n_cond_1_pmod_1_val_1, sess_n_cond_1_pmod_1_val_2, ..., sess_n_cond_1_pmod_1_val_n],
              [sess_n_cond_1_pmod_2_val_1, sess_n_cond_1_pmod_2_val_2, ..., sess_n_cond_1_pmod_2_val_n],
              ...,
              [sess_n_cond_1_pmod_n_val_1, sess_n_cond_1_pmod_n_val_2, ..., sess_n_cond_1_pmod_n_val_n]
            ],
            [
              [sess_n_cond_2_pmod_1_val_1, sess_n_cond_2_pmod_1_val_2, ..., sess_n_cond_2_pmod_1_val_n],
              [sess_n_cond_2_pmod_2_val_1, sess_n_cond_2_pmod_2_val_2, ..., sess_n_cond_2_pmod_2_val_n],
              ...,
              [sess_n_cond_2_pmod_n_val_1, sess_n_cond_2_pmod_n_val_2, ..., sess_n_cond_2_pmod_n_val_n]
            ],
            ...,
            [
              [sess_n_cond_n_pmod_1_val_1, sess_n_cond_n_pmod_1_val_2, ..., sess_n_cond_n_pmod_1_val_n],
              [sess_n_cond_n_pmod_2_val_1, sess_n_cond_n_pmod_2_val_2, ..., sess_n_cond_n_pmod_2_val_n],
              ...,
              [sess_n_cond_n_pmod_n_val_1, sess_n_cond_n_pmod_n_val_2, ..., sess_n_cond_n_pmod_n_val_n]
            ]
          ]
        ] # Use None if there is no parametric modulation for a condition or if there is no condition for a session.

        ex. [None, [None, [[36.4, 61.9, 105.1, 178.7], [7.2, 19.6, 65.9, 221.4]]]]

- *sess_cond_pmod_polys <=> sess.cond.pmod.poly* [#label]_
    The polynomial expansion used for the parametric modulation (a list of items which are a list of items which are a list of items which are 1 or 2 or 3 or 4 or 5 or 6). 1st order modulation would model the stick functions and a linear change of the stick function heights over different values of the parameter. Higher order modulation will introduce further columns that contain the stick functions scaled by parameter squared, cubed etc. In a general way, ``sess_cond_pmod_polys`` can be defined as:

    .. code-block:: python

        [
          [
            [sess_1_cond_1_pmod_1_poly, sess_1_cond_1_pmod_2_poly, ..., sess_1_cond_1_pmod_n_poly],
	    [sess_1_cond_2_pmod_1_poly, sess_1_cond_2_pmod_2_poly, ..., sess_1_cond_2_pmod_n_poly],
	    ...,
	    [sess_1_cond_n_pmod_1_poly, sess_1_cond_n_pmod_2_poly, ..., sess_1_cond_n_pmod_n_poly]
          ],
            [sess_2_cond_1_pmod_1_poly, sess_2_cond_1_pmod_2_poly, ..., sess_2_cond_1_pmod_n_poly],
	    [sess_2_cond_2_pmod_1_poly, sess_2_cond_2_pmod_2_poly, ..., sess_2_cond_2_pmod_n_poly],
	    ...,
	    [sess_2_cond_n_pmod_1_poly, sess_2_cond_n_pmod_2_poly, ..., sess_2_cond_n_pmod_n_poly]
          ],
          ...,
          [
            [sess_n_cond_1_pmod_1_poly, sess_n_cond_1_pmod_2_poly, ..., sess_n_cond_1_pmod_n_poly],
	    [sess_n_cond_2_pmod_1_poly, sess_n_cond_2_pmod_2_poly, ..., sess_n_cond_2_pmod_n_poly],
	    ...,
	    [sess_n_cond_n_pmod_1_poly, sess_n_cond_n_pmod_2_poly, ..., sess_n_cond_n_pmod_n_poly]
          ]
        ] # Use None if there is no parametric modulation for a condition or if there is no condition for a session

        ex. [None, [None, [1, 6]]]

- *sess_cond_orth <=> sess.cond.orth* [#label]_
    Orthogonalise regressors within trial types (a list of items which are a list of items which are 0 or 1). One orthogonalise modulation parameter for each condition must be applied.
      | - 0: no othogonalisation
      | - 1: othogonalisation

    In a general way, ``sess_cond_orth`` can be defined as:

    .. code-block:: python

        [
          [sess_1_cond_1_orth, sess_1_cond_2_orth, ..., sess_1_cond_n_orth],
          [sess_2_cond_1_orth, sess_2_cond_2_orth, ..., sess2_cond_n_orth],
          ...,
          [sess_n_cond_1_orth, sess_n_cond_2_orth ..., sess_n_cond_n_orth]
        ] # Use None if there is no condition for a session``

        ex. [None, [1, 0]]

- *sess_multi <=> sess.multi* [#label]_
    A \*.mat file containing details of the multiple experimental conditions for each session (a list of items which are a filename). This \*.mat file must include the following cell arrays (each 1 x n): names, onsets and durations. This option can be used to load in one go all the information that can also be given with the ``sess_cond_names``, ``sess_cond_onsets`` and ``sess_cond_durations``. The time, parametric and orthogonalise effects can also be included (please see spm documentation for more informations). In a general way ``sess_multi`` can be defined as:

    .. code-block:: python

        [sess_1_multi, sess_2_multi, ..., sess_n_multi] # Use None if there is no value for a session or [] if no value for all sessions

        ex. []

- *sess_regress <=> sess.regress* [#label]_
    Additional columns included in the design matrix, which may model effects that would not be convolved with the haemodynamic response (a list of items which are a list of items which are a dictionary with keys which are 'name' or 'val' and with values which are a string or a list of float). In a general way ``sess_regress`` can be defined as:

    .. code-block:: python

        [
          [sess_1_dict_1('name' = string, 'val' = list of float),
           sess_1_dict_2('name' = string, 'val' = list of float),
           ...,
           sess_1_dict_n('name' = string, 'val' = list of float)
          ],
          [sess_2_dict_1('name' = string, 'val' = list of float),
           sess_2_dict_2('name' = string, 'val' = list of float),
           ...,
           sess_2_dict_n('name' = string, 'val' = list of float)
          ],
          ...,
          [sess_n_dict_1('name' = string, 'val' = list of float),
           sess_n_dict_2('name' = string, 'val' = list of float),
           ...,
           sess_n_dict_n('name' = string, 'val' = list of float)
          ]
        ] # Use None if no value for a session or [[]] if no value for all sessions

        ex. [
              [
                {
                 'name': 'reg1',
                 'val': [0.79, 0.98, 0.585, 0.805, 0.53, 1.155, 0.66, 1.14, 0.195, 1.045, 0.82, 0.49, 0.765, 0.67, 0.0, 0.12, 0.955, 0.935, 0.26, 0.865]
	        },
                {
	         'name': 'reg2',
	         'val': [1.58, 1.96, 1.17, 1.61, 1.06, 2.31, 1.32, 2.28, 0.39, 2.09, 1.64, 0.98, 1.53, 1.34, 0.0, 0.24, 1.91, 1.87, 0.52, 1.73]
	        }
	      ],
	      None
	    ]

- *sess_multi_reg <=> sess.multi_reg* [#label]_
    The \*.mat/\*.txt file(s) containing details of multiple regressors (a list of items which are a filename). This option can be used to load in one go all the information that can also be given with the ``sess_regress`` parameter (please see spm documentation for more informations). In a general way ``sess_multi_reg`` can be defined as:

    .. code-block:: python

        [
          [sess_1_file_1, sess_1_file_2, ..., sess_1_file_n],
          [sess_2_file_1, sess_2_file_2, ..., sess_2_file_n],
          ...,
          [sess_n_file_1, sess_n_file_2, ..., sess_n_file_n]
        ] # Use None if no value for a session

        ex. [['/home/username/data/raw_data/file1.mat', '/home/username/data/raw_data/file2.txt'], None]

- *sess_hpf <=> sess.hpf* [#label]_
    High-pass filter (a list of items which are a float). One value per session must be applied. Slow signal drifts with a period longer than the applied value will be removed.

    ::

      ex. [128.0, 427.2]

- *factor_info <=> fact*
    If a factorial design exists, then SPM can automatically generate the contrasts necessary to test for the main effects and interactions. In factorial designs, a factor is a major independent variable and a level is a subdivision of a factor. Define for each factor a dictionary with its name and its number of levels (a list of items which are a dictionary with keys which are ‘name’ and ‘levels’ and with values which are a string or an integer).

    ::

      ex. []

- *bases <=> bases* [#label]_
    To define basic functions for modeling hemodynamic response (a None or a dictionary with keys which are ‘hrf’ or ‘fourier’ or ‘fourier_han’ or ‘gamma’ or ‘fir’ and with values which are a dictionary with keys which are ‘derivs’ or ‘length' or ‘order’ and with values which are a list of 2 elements or a float or an integer).

      | - hrf: Canonical Hemodynamic Response Function (HRF):

          | - derivs (model derivatives): The time derivative allows the peak response to vary by plus or minus a second and the dispersion derivative allows the width of the response to vary by a similar amount (2-element list)

              | - No derivatives: [0, 0]
              | - Time derivatives : [1, 0]
              | - Time and Dispersion derivatives: [1, 1]

      | - fourier (Fourier Set), fourier_han (Fourier Set Hanning), gamma (Gamma Functions), or fir (Finite Impulse Response):

          | - length: The length in seconds of the post-stimulus time window that the basis functions span (a float)
          | - order: Number of basis functions to use (a int)

    ::

      ex. {"hrf":{"derivs":[0,0]}}

- *volterra_expansion_order <=> volt* [#label]_
    Modeling using a Volterra series formulation (one of 1 or 2).
        | - 1: do not model interactions
        | - 2: model interactions

    ::

      ex. 1

- *global_intensity_normalization <=> global* [#label]_
    Global intensity normalisation with scaling or not (one of 'none' or 'scaling').

    ::

      ex. none

- *mask_threshold <=> mthresh* [#label]_
    Masking threshold, defined as proportion of globals (a float).

    ::

      ex. 0.8

- *mask_image <=> mask* [#label]_
    An image for explicitly masking the analysis (a pathlike object or string representing a file). Only those voxels will be analysed.

    ::

      ex. /home/username/data/raw_data/mask_swc1Anat.nii

- *model_serial_correlations <=> cvi* [#label]_
    Serial correlations in fMRI time series due to aliased biorhythms and unmodelled neuronal activity can be accounted for using an autoregressive AR(1) model during Classical (ReML) parameter estimation. Serial correlation can be ignored if the 'none' option is selected. FAST option is only available in SPM12 (one of AR(1), or FAST or none).

    ::

      ex. AR(1)

**Outputs parameters:**

- *spm_mat_file*
    The SPM.mat file (a pathlike object or string representing a file').

    ::

      ex. SPM.mat

-------------

.. [#label] Syntax: mia_processes/nipype Level1Design brick <=> SPM12 fMRI model specification.

	    Usefull links:
	    `SPM12 fMRI model specification <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=63>`_,
	    `nipype Level1Design <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.model.html#level1design>`_
..
  `nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.model.html#level1design>`_











                ```
