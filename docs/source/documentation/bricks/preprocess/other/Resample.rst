:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

.. line break
.. |br| raw:: html

   <br />

.. thin space
.. |ws1| raw:: html

   &thinsp;

.. em space

.. |ws2| raw:: html

   &emsp;

.. en space

.. |ws3| raw:: html

   &ensp;

.. non-breakable space

.. |ws4| raw:: html

   &nbsp;

==============
Resample brick
==============

Resamples an image to the resolution of a reference image, using N-D spline interpolation
-----------------------------------------------------------------------------------------

>>> from mia_processes.bricks.preprocess.other import Resample
>>> Resample.help()

**Inputs parameters:**

- *reference_image*
    File used as reference to resample the files_to_resample image (a pathlike object or a string representing a file with extension in [.img,
    .nii, .hdr]).

    ::

      ex. /home/ArthurBlair/data/raw_data/swrFunc.nii

- *files_to_resample*
    Files to apply transformation to (a list of pathlike object or string representing a file, valid extensions in [.img, .nii, .hdr]).
    ::

      ex. ['/home/ArthurBlair/data/raw_data/mask_swc1Anat_002.nii']

- *suffix_to_delete*
    The suffix to delete from the files_to_resample, when creating the out_files (a string).
    
    ::

      ex. _002

- *suffix*
    The suffix for the out_files image(s) (a string).

    ::

       ex. _003

- *prefix*
    The prefix for the out_files image(s) (a string).

    ::

      ex. <undefined>

- *interp*
    The order of the spline interpolation (an integer between 0 and 5).

    ::

      ex. 3

**Outputs parameters:**

- *out_files*
    The resulting image after resampling.

    ::

      ex. /home/ArthurBlair/data/raw_data/mask_swc1Anat_003.nii

-------------

NOTE:

- In addition to the resolution settings, it seems the size, the orientation and the position of the reference image is also transferred to the resampled image.
- Currently, this brick use `nibabel.processing.resample_from_to() <https://nipy.org/nibabel/reference/nibabel.processing.html#nibabel.processing.resample_from_to>`_. It seems this method produces `a little noise (non-zero intensity for the background) <https://mail.python.org/pipermail/neuroimaging/2019-January/001902.html>`_.
- Currently, resampling is only done for 3D images and using 3D or 4D images for reference.
	
