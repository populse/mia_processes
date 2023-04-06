:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==============
Resample brick
==============

Resamples an image to the resolution of a reference image, using N-D spline interpolation
-----------------------------------------------------------------------------------------

>>> from mia_processes.bricks.preprocess.other import Resample
>>> Resample.help()

**Mandatory inputs parameters:**

- *reference_image* (a string representing a file with extension)
    File used as reference to resample the files_to_resample image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. /home/username/data/raw_data/swrFunc.nii

- *files_to_resample*  (a list of pathlike object or string representing a file)
    Files to apply transformation to (valid extensions: [.nii, .nii.gz]).
    ::

      ex. ['/home/username/data/raw_data/mask_swc1Anat_002.nii']

**Optional inputs with default value parameters:**

- *interp* (an integer between 0 and 5, optional)
    The order of the spline interpolation.

    ::

      default value. 3

- *prefix* (a string, optional)
    The prefix for the out_files image(s) .

    ::

      default value. <undefined>

- *suffix_to_delete* (a string, optional)
    The suffix to delete from the files_to_resample, when creating the out_files.
    
    ::

      default value. _002

- *suffix* (a string, optional)
    The suffix for the out_files image(s).

    ::

       default value. _003

**Outputs parameters:**

- *out_files*
    The resulting image after resampling.

    ::

      ex. /home/username/data/raw_data/mask_swc1Anat_003.nii

-------------

NOTE:

- In addition to the resolution settings, it seems the size, the orientation and the position of the reference image is also transferred to the resampled image.
- Currently, this brick use `nibabel.processing.resample_from_to() <https://nipy.org/nibabel/reference/nibabel.processing.html#nibabel.processing.resample_from_to>`_. It seems this method produces `a little noise (non-zero intensity for the background) <https://mail.python.org/pipermail/neuroimaging/2019-January/001902.html>`_.
- Currently, resampling is only done for 3D images and using 3D or 4D images for reference.
	
