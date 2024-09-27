:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===============
FitTensor brick
===============

Diffusion tensor estimation
---------------------------

Convert diffusion-weighted images to tensor images.

| By default, the diffusion tensor (and optionally the kurtosis) is fitted to the log-signal in two steps:
|   - first fit is done using weighted least-squares (WLS) with weights based on the empirical signal intensities (or using ordinary least-squares (OLS) is "ols_option" is used)
|   - second fit is done using iterated weighted least-squares (IWLS) with weights determined by the signal predictions from the previous iteration (number of iteration could be choose with the "number_of_iteration" option).


(mrtrix dwi2tensor command)

*Please note that, in mia_processes, MRtrix bricks required diffusion data in MRtrix .mif format.
The MRConvert brick could be used to convert diffusion data in NIfTI format into MRtrix .mif format.
In this case, bvec and bval file should be specified.*

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input DWI image (valid extensions: [.mif]).

    ::

      ex. '/home/username/data/raw_data/DWI.mif'

**Optional inputs with default value parameters:**

- *estimate_dkt* (a boolean, default value is False, optional)
    Estimate diffusion kurtosis

    ::

      ex. False


- *get_predicted_signal* (a boolean, default value is False, optional)
    Get a file with the predicted signal from the tensor fits

    ::

      ex. False

- *get_output_b0* (a boolean, default value is False, optional)
    Get the put b0 image

    ::

      ex. False

- *ols_option* (a boolean, default value is False, optional)
    Perform initial fit using an ordinary least-squares (OLS) fit, that is, all measurements contribute equally to the fit (instead of using WLS)

    ::

      ex. False

- *number_of_iter* (an integer, default value is 2, optional)
    Number of iterative reweightings for IWLS algorithm.
    If 0 is set, only the first fitting will be done (WLS or OLS is ols_option is used).

    ::

      ex. 2


**Optional inputs parameters:**

- *in_mask* (a string representing an existing file, optional)
    Input mask image, only perform computation within the specified binary brain mas image (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/DWI_brainmask.mif'

**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    | The output diffusion tensor image (DTI).
    | The tensor coefficients are stored in the output image as follows:
    |   - volumes 0-5: D11, D22, D33, D12, D13, D23

    ::

      ex. '/home/username/data/derived_data/DWI_dti.mif'

- *out_dkt* (a pathlike object or string representing a file, optional)
    | The output diffusion kurtosis image (DKI).
    | The coefficients  are stored as follows:
    |  - volumes 0-2: W1111, W2222, W3333
    |  - volumes 3-8: W1112, W1113, W1222, W1333, W2223, W2333
    |  - volumes 9-11: W1122, W1133, W2233
    |  - volumes 12-14: W1123, W1223, W1233

    ::

      ex. '/home/username/data/derived_data/DWI_dki.mif'

- *out_b0* (a pathlike object or string representing a file, optional)
    The output b0 image

    ::

      ex. '/home/username/data/derived_data/DWI_b0.mif'

- *predicted_signal_file* (a pathlike object or string representing a file, optional)
    The output predicted dwi image

    ::

      ex. '/home/username/data/derived_data/DWI_dti_predicted_signal.mif'

-------------

Useful links:

`mrtrix dwi2tensor <https://mrtrix.readthedocs.io/en/latest/reference/commands/dwi2tensor.html>`_

`mrtrix dwi2tensor - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.reconst.html#fittensor>`_
