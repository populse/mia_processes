:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

.. |br| raw:: html

   <br />

.. |ws| raw:: html

   &thinsp;

=============
Realign brick
=============

Realigns a time-series of images acquired from the same subject using a least squares approach and a 6 parameters (rigid body) spatial transformation
-----------------------------------------------------------------------------------------------------------------------------------------------------

>>> from mia_processes.preprocess.spm import Realign
>>> Realign.help()

**Inputs parameters:** [#label]_

- *in_files <=> data*
    A list of items which are an existing, uncompressed file (valid extensions: [.img, .nii, .hdr]).

    ::

      ex. ['/home/ArthurBlair/data/raw_data/Func.nii']

- *jobtype*
    One of 'estwrite', 'estimate' or 'write'.

      | estimate: generates realignment_parameters and modified_in_files
      | write: |ws| |ws| |ws| |ws| with write_which == [2, 0] or [1, 0] generates realigned_files
      | |ws| |ws| |ws| |ws| |ws| |ws| |ws| |ws| |ws| |ws| with write_which == [2, 1] generates mean_image and realigned_files
      | |ws| |ws| |ws| |ws| |ws| |ws| |ws| |ws| |ws| |ws| with write_which == [0, 1] generates mean_image
      | estwrite: |ws| with write_which == [2, 0] or [1, 0] generates realignment_parameters, modified_in_files and realigned_files
      | |ws| |ws| |ws| |ws| |ws| |ws| |ws| |ws| |ws| |ws| with write_which == [2, 1] generates realignment_parameters, modified_in_file, mean_image and realigned_files
      | |ws| |ws| |ws| |ws| |ws| |ws| |ws| |ws| |ws| |ws| with write_which == [0, 1] generates realignment_parameters, modified_in_file and mean_image

    ::

      ex. estwrite

- *quality <=> eoptions.quality*
    Quality versus speed trade-off (0.0 <= a floating point number <= 1.0). Highest quality (1)
    gives most precise results, whereas lower qualities gives faster realignment.

    ::

     ex. 0.9

- *separation <=> eoptions.sep*
    Sampling separation in mm in the reference image (a floating point number >= 0.0).
    Smaller sampling distances gives more accurate results, but will be slower.

    ::

      ex. 4.0

- *fwhm <=> eoptions.fwhm*
    The gaussian smoothing kernel width (mm, a floating point number >= 0.0) applied to the images
    before estimating the realignment parameters.

    ::

      ex. 5.0

 - *register_to_mean <=> eoptions.rtm*
     Indicate whether realignment is done to the mean image (True) or to the first
     image (False).

     ::

       ex. True

 - *interp <=> eoptions.interp*
     Degree of b-spline (1 <= a long integer <= 7) used for interpolation. Higher degree
     interpolation methods provide the better interpolation, but they are slower because they use
     more neighbouring voxels.

     ::

       ex. 2

 - **

 - **

 - **

**Outputs parameters:** [#label]_

- *realigned_files*
    If the write_which parameter is equal to [2, 0], [1, 0] or [2, 1] and jobtype parameter is equal to write or estwrite,
    these will be the resliced files (a list of items which are a list of items which are a pathlike object or string
    representing an existing file or a pathlike object or string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/raw_data/rFunc.nii

-------------

.. [#label] Syntax: mia_processes/nipype Realign <=> SPM12 Realign.
	    
	    Usefull links:
	    `SPM12 Realign <https://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=25>`_, 
	    `nipype <https://nipype.readthedocs.io/en/latest/interfaces/generated/interfaces.spm/preprocess.html#realign>`_
