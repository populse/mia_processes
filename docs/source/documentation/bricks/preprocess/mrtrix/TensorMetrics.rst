:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===================
TensorMetrics brick
===================

Compute metrics from tensors
----------------------------

Generate maps of tensor-derived parameters as the mean apparent diffusion coefficient (ADC), the fractional anisotropy (FA), the axial diffusivity (AD),
the radial diffusivity (RD), the linearity metric, the planarity metric, the sphericity metric, the selected eigenvalue(s) of the diffusion tensor or
the selected eigenvector(s) of the diffusion tensor
,
(mrtrix tensor2metric command)

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a pathlike object or a string representing an existing file)
    Input tensor image (valid extensions: [.mif]).

    ::

      ex. '/home/username/data/derived_data/DWI_dti.mif'


**Optional inputs with default value parameters:**

- *modulate* (FA or none or eigval, default value is FA, optional)
    Specify how to modulate the magnitude of the eigenvectors.

    ::

      ex. FA

- *get_ad* (a boolean, default value is False, optional)
    Compute the axial diffusivity (AD) of the diffusion tensor.

    ::

      ex. False

- *get_adc* (a boolean, default value is True, optional)
    Compute the mean apparent diffusion coefficient (ADC) of the diffusion tensor.

    ::

      ex. True

- *get_cl* (a boolean, default value is False, optional)
    Compute the linearity metric of the diffusion tensor (one of the three Westin shape metrics).

    ::

      ex. False

- *get_cp* (a boolean, default value is False, optional)
    Compute the planarity metric of the diffusion tensor (one of the three Westin shape metrics).

    ::

      ex. False

- *get_cs* (a boolean, default value is False, optional)
    Compute the sphericity metric of the diffusion tensor (one of the three Westin shape metrics).

    ::

      ex. False

- *get_fa* (a boolean, default value is True, optional)
    Compute the fractional anisotropy (FA)  of the diffusion tensor.

    ::

      ex. True

- *get_rd* (a boolean, default value is False, optional)
    Compute the radial diffusivity (RD) of the diffusion tensor.

    ::

      ex. False

- *get_value* (a boolean, default value is False, optional)
    Compute the selected eigenvalue(s) of the diffusion tensor.

    ::

      ex. False

- *get_vector* (a boolean, default value is True, optional)
    Computethe selected eigenvector(s) of the diffusion tensor.

    ::

      ex. True

**Outputs parameters:**

- *ad_file* (a pathlike object or string representing a file, optional)
    Output AD file

    ::

      ex. '/home/username/data/derived_data/DWI_dti_ad.mif'

- *adc_file* (a pathlike object or string representing a file, optional)
    Output ADC file

    ::

      ex. '/home/username/data/derived_data/DWI_dti_adc.mif'

- *cl_file* (a pathlike object or string representing a file, optional)
    Output linearity metric file

    ::

      ex. '/home/username/data/derived_data/DWI_dti_cl.mif'

- *cp_file* (a pathlike object or string representing a file, optional)
    Output planarity metric file

    ::

      ex. '/home/username/data/derived_data/DWI_dti_cp.mif'

- *cs_file* (a pathlike object or string representing a file, optional)
    Output sphericity metric file

    ::

      ex. '/home/username/data/derived_data/DWI_dti_cs.mif'

- *fa_file* (a pathlike object or string representing a file, optional)
    Output FA file

    ::

      ex. '/home/username/data/derived_data/DWI_dti_fa.mif'

- *rd_file* (a pathlike object or string representing a file, optional)
    Output RD file

    ::

      ex. '/home/username/data/derived_data/DWI_dti_rd.mif'

- *value_file* (a pathlike object or string representing a file, optional)
    Output selected eigenvalue(s) file

    ::

      ex. '/home/username/data/derived_data/DWI_dti_value.mif'

- *vector_file* (a pathlike object or string representing a file, optional)
    Output selected eigenvector(s) file

    ::

      ex. '/home/username/data/derived_data/DWI_dti_vector.mif'

-------------

Usefull links:

`mrtrix tensor2metric <https://mrtrix.readthedocs.io/en/latest/reference/commands/tensor2metric.html>`_

`mrtrix tensor2metric - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.utils.html#tensormetrics>`_
