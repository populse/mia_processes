:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

================
ComputeDKI brick
================

Reconstruction of the diffusion signal with the kurtosis tensor model
---------------------------------------------------------------------

The diffusion kurtosis imaging (DKI) model is an expansion of the diffusion tensor imaging (DTI) model
that allows quantification of the degree to which water diffusion in biological tissues is non-Gaussian.

This brick used functions proposed by Dipy to reconstruct the diffusion signal with the kurtosis tensor model.

Since the diffusion kurtosis model estimates the diffusion tensor, all standard diffusion tensor statistics can be computed:
the fractional anisotropy (FA), the mean diffusivity (MD), the axial diffusivity (AD) and the radial diffusivity (RD).

In addition to the standard diffusion statistics, this brick can be used to estimate the non-Gaussian measures of mean kurtosis (MK),
the axial kurtosis (AK) and the radial kurtosis (RK).

The mean of the kurtosis tensor (mKT) and the kurtosis fractional anisotropy (kFA) are also computed. These measures only depend on the kurtosis tensor.

--------------------------------------

**Mandatory inputs parameters:**

- *in_dwi* (a string representing an existing file)
    Diffusion file (valid extensions: [.nii, .nii.gz]).
    It should be multi-shell data, i.e. data acquired from more than one non-zero b-value.

    ::

      ex. '/home/username/data/raw_data/dwi.nii'

- *dwi_bvec* (a string representing an existing file)
    Bvec file (valid extensions: [.bvec]).

    ::

      ex. '/home/username/data/raw_data/dwi.bvec'

- *dwi_bval* (a string representing an existing file)
    Bval file (valid extensions: [.bval]).

    ::

      ex. '/home/username/data/raw_data/dwi.bval'

- *in_mask* (a string representing an existing file)
    Brain mask file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/dwi_brainmask.nii'


**Outputs parameters:**

- *out_FA* (a strings representing a file)
    The fractional anisotropy (FA) image

    ::

      ex. '/home/username/data/derived_data/dwi_dki_FA.nii'

- *out_MD* (a strings representing a file)
    The mean diffusivity (MD) image

    ::

      ex. '/home/username/data/derived_data/dwi_dki_MD.nii'

- *out_RD* (a strings representing a file)
    The radial diffusivity (RD) image

    ::

      ex. '/home/username/data/derived_data/dwi_dki_RD.nii'

- *out_AD* (a strings representing a file)
    The axial diffusivity (AD) image

    ::

      ex. '/home/username/data/derived_data/dwi_dki_AD.nii'

- *out_MK* (a strings representing a file)
    The mean kurtosis (MK) image

    ::

      ex. '/home/username/data/derived_data/dwi_dki_MK.nii'

- *out_RK* (a strings representing a file)
    The radial kurtosis (RK) image

    ::

      ex. '/home/username/data/derived_data/dwi_dki_RK.nii'

- *out_AK* (a strings representing a file)
    The axial kurtosis (AK) image

    ::

      ex. '/home/username/data/derived_data/dwi_dki_AK.nii'

- *out_mKT* (a strings representing a file)
    The mean of the kurtosis tensor (mKT) image

    ::

      ex. '/home/username/data/derived_data/dwi_dki_mKT.nii'

- *out_kFA* (a strings representing a file)
    The kurtosis fractional anisotropy (kFA) image

    ::

      ex. '/home/username/data/derived_data/dwi_dki_kFA.nii'

-------------

Usefull links:

`Dipy Reconstruction of the diffusion signal with the kurtosis tensor model <https://dipy.org/documentation/1.7.0/examples_built/07_reconstruction/reconst_dki/#sphx-glr-examples-built-07-reconstruction-reconst-dki-py>`_

`Jensen JH 2005 <https://pubmed.ncbi.nlm.nih.gov/15906300/>`_
