:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Denoise brick
============

Denoise image using Non-Local Means algorithm (NLMEANS).
The value of a pixel is replaced by an average of a set of other pixel values: the specific patches centered on the other pixels are contrasted to the patch centered on the pixel of interest, a
nd the average only applies to pixels with patches close to the current patch.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input file to denoise (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'

**Optional inputs with default value parameters:**

- *block_radius* (an integer, optional, default value is 5)
    The block size will be 2 x block_radius + 1.

    ::

      ex. 5

- *noise_model* (rician or gaussian, optional, default value is rician)
    Noise distribution model.

    ::

      ex. rician

- *out_prefix* (a string, optional, default value is 'denoise')
    Prefix of the output image.

    ::

      ex. 'denoise_'

- *patch_radius* (an integer, optional, default value is 1)
    The patch size will be 2 x patch_radius + 1.

    ::

      ex. 1

  **Optional inputs parameters:**

- *in_mask* (a string representing an existing file, optional)
    Brain mask.

    ::

      ex. '/home/username/data/derived_data/func_brain_mask.nii'

- *noise_mask* (a string representing an existing file, optional)
    Mask in which the mean signal will be computed.

    ::

      ex. '/home/username/data/derived_data/func_brain_mask.nii'

- *signal_mask* (a string representing an existing file, optional)
    Mask in which the standard deviation of noise will be computed

    ::

      ex. '/home/username/data/derived_data/func_brain_mask.nii'

- *snr* (a float, optional)
    Set manually Signal to Noise Ratio.
    Default is Undefined (ie parameter not used).

    ::

      ex. 260.0

**Outputs parameters:**

- *out_file* (a strings representing a file)
    Out image (extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/denoise_func.nii'

-------------

Usefull links:

`Dipy Denoise - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.dipy.preprocess.html#denoise>`_
`Dipy Denoise <https://dipy.org/documentation/1.6.0./examples_built/denoise_nlmeans/#example-denoise-nlmeans>`_
`NLMEANS article 1 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2881565/>`_
`NLMEANS article 2 <https://hal.science/hal-00645538/document>`_
