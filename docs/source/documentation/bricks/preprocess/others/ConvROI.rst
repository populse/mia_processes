:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=============
ConvROI brick
=============

Image convolution with one image
--------------------------------

- Resample the `images_to_convolve` to the size of `convolve_with`.
- Then, convolve each resized `images_to_convolve` with `convolve_with`.
- The “PatientName_data/ROI_data/convROI_BOLD” directory is created to receive the convolution results from the runtime.
  If this directory exists at runtime, it is overwritten.
- To work correctly, the database entry for the `convolve_with` parameter must have the "PatientName" tag filled in.
- To eliminate the near-zero noise generated in the resized image, thresholding (10-5) is performed.
  All values below 10-5 are therefore truncated at 0.

--------------------------------------

**Mandatory inputs parameters:**

- *images_to_convolve*
    A list of images to convolve with `convolve_with` image (a list of image paths).

    ::

      ex. ['/home/username/data/raw_data/ACA_L.nii',
           '/home/username/data/raw_data/ACA_R.nii',
	   '/home/username/data/raw_data/ACM_L.nii',
	   '/home/username/data/raw_data/ACM_R.nii']

- *convolve_with* (a string representing an existing file)
    An image used to convolve with `images_to_convolve` (a string representing a path to a file).

    ::

      ex. '/home/username/data/raw_data/mask.nii'


- *prefix*
    Specify the string to be prepended to `out_images` (a string).

    ::

      ex. conv

**Outputs parameters:**

- *out_images*
    The convolved images (a list of pathlike objects or strings representing a file).

    ::

      ex. [`/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convACA_L.nii',
           `/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convACA_R.nii',
           `/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convACM_L.nii',
           '/home/username/data/derived_data/patient-name_data/ROI_data/convROI_BOLD/convACM_R.nii']
