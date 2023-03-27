:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
TShift brick
============

Slice-time correction of bold images (using mean option of the AFNI 3dTshift command).
Shifts voxel time series from input so that separate slices are aligned to the same temporal origin.

You can either used "slice_timing" parameter or "tpattern" parameter to specify slice time. 
If neither of this two parameters are defined, the process return the input file. 


--------------------------------------

>>> from mia_processes.bricks.preprocess.others import TShift

>>> TShift.help()

**Inputs parameters:**

- *in_file* (a string representing an existing file)
    Input bold file to be time-shifted (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

- *slice_encoding_dir* (k or k-, optional)
    Direction in which slice_timing is specified. 
    If negative, slice_timing is defined in reverse order, that is, the first entry corresponds to the slice with the largest index and the final entry corresponds to slice index zero.
    Default is K. 

    ::

      ex. k


- *slice_timing* (a string representing an existing file or a list of floats, optional)
    Time offsets from the volume acquisition onset for each slice. 
    Mutually exclusive with "tpattern" parameter.

    ::

      ex. slice_timing.1D

- *tpattern* (alt+z or altplus or alt+z2 or alt-z or altminus or alt-z2 or seq+z or seqplus or seq-z or seqminus, optional)
    Use specified slice time pattern rather than one in header. 
    One of ‘alt+z’ or ‘altplus’ or ‘alt+z2’ or ‘alt-z’ or ‘altminus’ or ‘alt-z2’ or ‘seq+z’ or ‘seqplus’ or‘seq-z’ or ‘seqminus’.
    Mutually exclusive with "slice_timming" parameter.
    Default is Undefined (ie parameter not used).

    ::

      ex. False

- *ignore* (an integer, optional)
    Ignore the first set of points specified.
    The first ii values will be unchanged in the output.They also will not be used in the detrending or time shifting.
    Default is Undefined (ie parameter not used).

    ::

      ex. 2

- *interpolation* (‘Fourier’ or ‘linear’ or ‘cubic’ or ‘quintic’ or ‘heptic’, optional)
    Different interpolation methods:
      -Fourier = Use a Fourier method (the default: most accurate; slowest).
      -linear  = Use linear (1st order polynomial) interpolation (least accurate).
      -cubic   = Use the cubic (3rd order) Lagrange polynomial interpolation.
      -quintic = Use the quintic (5th order) Lagrange polynomial interpolation.
      -heptic  = Use the heptic (7th order) Lagrange polynomial interpolation.
    Default is 'Fourier'.

    ::

      ex. Fourier

- *rlt* (a boolean, optional)
    Before shifting, remove the mean and linear trend.
    Default is False.

    ::

      ex. False

- *rlt* (a boolean, optional)
    Before shifting, remove the mean and linear trend and later put back the mean.
    Default is False.

    ::

      ex. False

- *tr* (a string, optional)
    Manually set the TR. Add suffix “s” for seconds or “ms” for milliseconds. 
    Default is Undefined (ie parameter not used).

    ::

      ex. 2.5s

- *tslice* (an integer, optional)
    Align each slice to time offset of given slice.
    Mutually exclusive with tzero parameter.
    Default is Undefined (ie parameter not used).

    For example if tslice = 5, each slice will be align to time offset of slice 5.

    ::

      ex. 5

- *tzero* (an float, optional)
    Align each slice to the given time offset.
    The value  must be between the minimum and maximum slice temporal offsets.
    Mutually exclusive with tslice parameter.
    Default is Undefined (ie parameter not used).

    Note that the default alignment time (when tzero is not used) is the average of the 'tpattern' values 
    (either from tpattern option or slice_timing option).

    For example if tzero = 0.0, each slice will be align to time offset 0.0.

    ::

      ex. 0.0

- *output_type* (NIFTI or NIFTI_GZ, optional)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *out_prefix* (a string, optional)
    Prefix of the output image. Default is 'st_corr_'.
    
    ::

        ex. 'st_corr_'

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out image (extensions: [.nii, .nii.gz]).
    
    ::

      ex. '/home/username/data/derived_data/st_corr_func.nii'

-------------

Usefull links:

`AFNI 3dTshift <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html>`_
`AFNI TShift - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.afni.preprocess.html#tshift>`_
