:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
EpiReg brick
============

Register EPI images to structural images
----------------------------------------

Run epi_reg script (FSL) to register EPI images to structural images.
To use this method, it should be possible to segment the structural image
in order to get a good white matter boundary and the EPI should contain some
intensity contrast between white matter and grey matter.

It is also possible to add fiedmaps to perform simultaneous registration and
EPI distortion-correction.
To use this option you should provided fmap, fmapmag and fmapmagbrain
parameters and information about EPI (see echospacing and pedir parameters).

--------------------------------------------------------------

**Mandatory inputs parameters:**

- *in_epi* (a string representing an existing file)
    Input EPI image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/bold.nii'

- *in_t1* (a string representing an existing file)
    Input structural image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

- *in_t1* (a string representing an existing file)
    Input brain extracted structural image (valid extensions:
    [.nii, .nii.gz]).

    For example, brain extracted image obtained with
    `BetSurfacesExtraction brick <./BetSurfacesExtraction.html>`_

    ::

      ex. '/home/username/data/raw_data/T1w_brain.nii'

**Optional inputs with default value parameters:**

- *suffix* (a string, optional, default value is epi2struct)
    Output base name.

    ::

      ex. 'epi2struct'

- *output_type* ('NIFTI' or 'NIFTI_GZ', optional, default value is NIFTI)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   - NIFTI: \*.nii
    |   - NIFTI_GZ: \*.nii.gz

    ::

      ex. NIFTI

- *no_clean* (a boolean, optional, default value is False)
    Do not clean up intermediate files.

    ::

      ex. False

- *no_fmapreg* (a boolean, optional, default value is False)
    If fmap available, do not perform registration of fmap to T1.

    ::

      ex. False

**Optional inputs :**

- *wmseg* (a string representing an existing file, optional)
    White matter segmentation of T1 image (valid extensions: [.nii, .nii.gz]).

    If no white-matter segmentation is given, the script will run FAST (fsl) to get one.

    ::

      ex. '/home/username/data/raw_data/T1w_wmseg.nii'


- *weight_image* (a string representing an existing file, optional)
    Weighting image (in T1 space) (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/weight_space-T1w.nii'

- *fmap* (a string representing an existing file, optional)
    Fieldmap image (in rad/s) (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/fmap.nii'

- *fmapmag* (a string representing an existing file, optional)
    Fieldmap magnitude image - wholehead(valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/fmap_mag.nii'

- *fmapmagbrain* (a string representing an existing file, optional)
    Fieldmap magnitude image - brain extracted (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/fmap_mag.nii'

- *pedir* (either x, y, z, -x, -y, -z, optional, default is Undefined)
    Phase encoding direction of the EPI

    ::

      ex. x

- *echospacing* (a float, optional, default is Undefined)
    Effective EPI echo spacing.

    ::

      ex. 0.3

**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    Unwarped and coregistered epi input.

    ::

      ex. '/home/username/data/derived_data/epi2struct.nii'

- *epi2str_ma* (a pathlike object or string representing a file)
    Rigid epi-to-structural transform.

    ::

      ex. '/home/username/data/derived_data/epi2struct.mat'

- *epi2str_inv* (a pathlike object or string representing a file, optional)
    Rigid structural-to-epi transform.

    ::

      ex. '/home/username/data/derived_data/epi2struct_inv.mat'

- *out_1vol* (a pathlike object or string representing a file, optional)
    Unwarped and coregistered epi input.

    ::

      ex. '/home/username/data/derived_data/epi2struct_1vol.nii'

- *fmap2epi_mat* (a pathlike object or string representing a file, optional)
    Rigid fieldmap-to-epi transform.

    ::

      ex. '/home/username/data/derived_data/epi2struct_fieldmaprads2epi.mat'

- *fmap2str_mat* (a pathlike object or string representing a file, optional)
    Rigid fieldmap-to-structural transform.

    ::

      ex. '/home/username/data/derived_data/epi2struct_fieldmap2str.mat'


- *fmap_epi* (a pathlike object or string representing a file, optional)
    Fieldmap in epi space.

    ::

      ex. '/home/username/data/derived_data/epi2struct_fieldmaprads2epi.nii'

- *fmap_str* (a pathlike object or string representing a file, optional)
    Fieldmap in structural space.

    ::

      ex. '/home/username/data/derived_data/epi2struct_fieldmaprads2str.nii'

- *fmapmag_str* (a pathlike object or string representing a file, optional)
    Fieldmap magnitude image in structural space.

    ::

      ex. '/home/username/data/derived_data/epi2struct_fieldmap2str.nii'

- *fullwarp* (a pathlike object or string representing a file, optional)
    Warpfield to unwarp epi and transform into structural space.
    ::

      ex. '/home/username/data/derived_data/epi2struct_warp.nii'


- *shiftmap* (a pathlike object or string representing a file, optional)
    Shiftmap in epi space.
    ::

      ex. '/home/username/data/derived_data/epi2struct_fieldmaprads2epi_shift.nii'

- *seg* (a pathlike object or string representing a file, optional)
    White matter, gray matter, csf segmentation.
    ::

      ex. '/home/username/data/derived_data/epi2struct_fast_seg.nii'

- *wmedge* (a pathlike object or string representing a file, optional)
    White matter edges for visualization.
    ::

      ex. '/home/username/data/derived_data/epi2struct_fast_wmedge.nii'

- *wmseg_out* (a pathlike object or string representing a file, optional)
    White matter segmentation used in flirt bbr.
    ::

      ex. '/home/username/data/derived_data/epi2struct_fast_wmseg.nii'

-------------

Useful links:

`FSL epi_reg <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT/UserGuide#epi_reg>`_

`FSL EpiReg - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.fsl.epi.html#epireg>`_
