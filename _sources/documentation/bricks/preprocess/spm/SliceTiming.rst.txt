:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

=================
SliceTiming brick
=================

Bring back every slice of a BOLD image to the reference slice time.
-------------------------------------------------------------------

>>> from mia_processes.bricks.preprocess.spm import SliceTiming
>>> SliceTiming.help()

**Inputs parameters:**

- *in_files <=> scans* [#label]_
    Path of the scans for processing (a list of items which are uncompressed files; valid extensions : [.nii]). The same parameters specified below will be applied to all data in the `in_files` parameter. The user will therefore make sure that all data in `in_files` are compatible (there is currently no verification).

    ::

      ex. ['/home/ArthurBlair/data/raw_data/Func.nii']

- *acquisition*
    Type of the acquisition, either sequential ascending, sequential descending, interleaved (middle-top), interleaved (bottom-up) or interleaved (top-down). Slice ordering is assumed to be from foot to head and bottom slice = 1.

    ::

      ex. sequential ascending

- *num_slices <=> nslices* [#label]_
    Number of slices per volume (an integer). If this parameter is undefined at the initialisation time, its value is looked in the database.

    ::

      ex. 53

- *TR <=> tr* [#label]_
    Repetition time (a float, in seconds). This is the time between volume acquisitions. If this parameter is undefined at the initialisation time, its value is retrieved from the database.

    ::

      ex. 3.000002

- *TA <=> ta* [#label]_
    Time of volume acquisition (a float, in seconds). If this parameter is undefined at the initialisation time, its value is automatically calculated as TR-(TR/num_slices). Must be set to 0 for the slice_order and ref_slice parameters to be automatically calculated in units of time (ms).

    ::

      ex. 2.9433981886792453

- *slice_order <=> so* [#label]_
    Time of volume acquisition (a float, in seconds). If this parameter is undefined at the initialisation time, its value is automatically calculated as TR-(TR/num_slices). Must be set to 0 for the slice_order and ref_slice parameters to be automatically calculated in units of time (ms).

    ::

      ex. [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53]

- *ref_slice <=> refslice* [#label]_
    The chosen reference slice (an integer or a float or an _Undefined or None). It will be used as a time reference to correct the differences in the acquisition time of the slices. According to the acquisition parameter value, ref_slice is initialised to num_slices/2 (acquisition == sequential ascending or acquisition == sequential descending) or 1 (other value of acquisition). If TA is set to 0 then ref_slice will be automatically set in ms, otherwise ref_slice will be automatically calculated in slice index.

    ::

      ex. 27

- *out_prefix <=> prefix* [#label]_
    Prefix of the output image (a string).

    ::

      ex. a

**Outputs parameters:**

- *timed_files*
    The image after the SliceTiming correction (a pathlike object or string representing a file).

    ::

      ex. /home/ArthurBlair/data/derived_data/aFunc.nii

-------------

.. [#label] Syntax: mia_processes/nipype SliceTiming <=> SPM12 Slice Timing.

	    Usefull links:
	    `SPM12 Slice Timing <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=19>`_, 
	    `nipype SliceTiming <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.spm.preprocess.html#slicetiming>`_
