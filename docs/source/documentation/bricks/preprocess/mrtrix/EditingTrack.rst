:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==================
EditingTrack brick
==================

Perform various editing operations on track files
-------------------------------------------------

(mrtrix tckedit command)

--------------------------------------

**Mandatory inputs parameters:**

- *in_tracks* (a list of items which are a pathlike object a string representing an existing file)
    Inputs track file(s)(valid extensions: [.tck]).

    ::

      ex. ['/home/username/data/raw_data/tracks.tck']

**Optional inputs with default value parameters:**

- *suffix* (a string, default value is edited, optional)
    Output file suffix

    ::

      ex. 'edited'


- *inverse* (a boolean, default value is False, optional)
    Output the inverse selection of streamlines based on the criteria provided;
    i.e. only those streamlines that fail at least one selection criterion, and/or vertices
    that are outside masks if provided, will be written to file

    ::

      ex. False

- *ends_only* (a boolean, default value is False, optional)
    Only test the ends of each streamline against the provided include/exclude ROIs.

    ::

      ex. False

- *get_tck_weights_out_desc* (a boolean, default value is False, optional)
    Get an output text scalar file containing streamline weights.

    ::

      ex. False

**Optional inputs parameters:**

- *roi_excl* (string representing an existing file or a tuple of the form: (a float, a float, a float, a float), optional)
    Specify an exclusion region of interest, streamlines that enter ANY exclude region will be discarded.

    ::

      ex. '/home/username/data/derived_data/roi.mif'

- *roi_incl* (string representing an existing file or a tuple of the form: (a float, a float, a float, a float), optional)
    Specify an inclusion region of interest, streamlines must traverse ALL inclusion regions to be accepted.

    ::

      ex. '/home/username/data/derived_data/roi.mif'

- *roi_incl_ordered* (string representing an existing file or a tuple of the form: (a float, a float, a float, a float), optional)
    Specify an inclusion region of interest, streamlines must traverse ALL inclusion_ordered regions in the order
    they are specified in order to be accepted.

    ::

      ex. '/home/username/data/derived_data/roi.mif'


- *roi_mask* (string representing an existing file or a tuple of the form: (a float, a float, a float, a float), optional)
    Specify a masking region of interest. If defined,streamlines exiting the mask will be truncated.

    ::

      ex. '/home/username/data/derived_data/DWI_brainmask.mif'

- *maxlength* (a float, optional)
    The maximum length of any streamline in mm.

    ::

      ex. 26.0

- *minlength* (a float, optional)
    The minimum length of any streamline in mm.

    ::

      ex. 2.0

- *number* (an integer, optional)
    The desired number of selected streamlines to be propagated to the output file.

    ::

      ex. 20000

- *skip* (an integer, optional)
    Omit this number of selected streamlines before commencing writing to the output file.

    ::

      ex. 500

- *maxweight* (an integer, optional)
    The maximum weight of any streamline

    ::

      ex. 10

- *minweight* (an integer, optional)
    The minimum weight of any streamline

    ::

      ex. 2

- *tck_weights_in* (string representing an existing file, optional)
    Specify a text scalar file containing the streamline weights

    ::

      ex. '/home/username/data/derived_data/tck_weight.txt'


**Outputs parameters:**

- *tracks_out* (a pathlike object or string representing a file)
    The output track file

    ::

      ex. '/home/username/data/derived_data/tracks_edited.tck'

- *tck_weights_out* (a pathlike object or string representing a file, optional)
    Output text scalar file containing streamline weights

    ::

      ex. '/home/username/data/derived_data/tracks_tck_weight.txt'


-------------

Usefull links:

`mrtrix tckedit <https://mrtrix.readthedocs.io/en/latest/reference/commands/tckedit.html>`_
