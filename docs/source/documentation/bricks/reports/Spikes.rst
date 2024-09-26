:orphan:

.. toctree::

+--------------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+-------------------------------------------+----------------------------------------------------+

============
Spikes brick
============

Computes the number of spikes
-----------------------------

Adapted from mriqc functional workflow (`spikes_mask function <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/functional.py#L939>`_
and `Spikes class <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/functional.py#L223>`_).

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input bold image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

**Optional inputs with default value parameters:**

- *detrend* (a boolean, optional, default value is False)
    Detrend data.

    ::

      ex. False

- *no_zscore* (a boolean, optional, default value is True)
    Do not zscore

    ::

      ex. True

- *out_prefix* (a string, optional, default value is 'spikes')
    Prefix of the output image.

    ::

      ex. 'spikes_'

- *skip_frames* (an integer, optional, default value is 0)
    Number of frames to skip in the beginning of the time series.

    ::

      ex. 0

- *spike_thresh* (a float, optional, default value is 6.0)
    z-score to call one timepoint of one axial slice a spike.

    ::

      ex. 6.0



**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out file with all the spkies detected.

    ::

      ex. '/home/username/data/derived_data/spikes_func.out'

-------------

Useful links:

`mriqc spikes <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/functional.py#L939>`_

`mriqc spikes mask <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/functional.py#L939>`_
