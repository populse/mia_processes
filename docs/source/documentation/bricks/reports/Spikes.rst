:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Spikes brick
============

Computes the number of spikes. 

Adapted from mriqc functional workflow (`spikes <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/functional.py#L939>`_
and `spikes mask <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/functional.py#L939>`_)


--------------------------------------

>>> from mia_processes.reports.preprocess import Spikes

>>> Spikes.help()

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input bold image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

**Optional inputs with default value parameters:**

- *detrend* (a boolean, optional)
    Detrend data.  
    
    ::

      default value. False

- *no_zscore* (a boolean, optional)
    Do not zscore 
    
    ::

      default value. True

- *out_prefix* (a string, optional)
    Prefix of the output image.

    ::

      default value. 'spikes_'

- *skip_frames* (an integer, optional)
    Number of frames to skip in the beginning of the time series.
    
    ::

      default value. 0

- *spike_thresh* (a float, optional)
    z-score to call one timepoint of one axial slice a spike.
    
    ::

      default value. 6.0



**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out file with all the spkies detected.
    
    ::

      ex. '/home/username/data/derived_data/spikes_func.out'

-------------

Usefull links:

`mriqc spikes <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/functional.py#L939>`_
`mriqc spikes mask <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/functional.py#L939>`_
