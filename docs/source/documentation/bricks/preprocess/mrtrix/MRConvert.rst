:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===============
MRConvert brick
===============

Perform conversion between different file types and optionally extract a subset of the input image
--------------------------------------------------------------------------------------------------

*Please note that, in mia_processes, MRtrix bricks required diffusion data in MRtrix .mif format.
This brick could be used to convert diffusion data in NIfTI format into MRtrix .mif format.
In this case, bvec and bval file should be specified.*

(mrtrix mrconvert command)

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a pathlike object or a string representing an existing file)
    Input image (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/DWI.nii.gz'

**Optional inputs with default value parameters:**

- *out_file_format* (NIFTI, NIFTI_GZ or MIF, optional, default value is MIF, optional)
    Format of the output image (NIFTI, NIFTI_GZ or MIF)

    ::

      ex. MIF

- *export_bvec_bval* (a boolean, default is False, optional)
    Export bvec / bval files in FSL format (if bvec/bval in the MIF image)

    ::

      ex. False


**Optional inputs:**

- *axes* (a list of items which are an integer, optional)
    Specify the axes from the input image that will be used to form the output image
    ::

      ex. [0, 1, 2] --> retain the first 3 axes


- *coord* (a list of items which are a tuple of the form (axis, selection), optional)
    Extract data from the input image only at the coordinates specified in the selection along the specified axis.

    ::

      ex. (3, 0) --> extract from axis number 3 (which is the fourth axis since counting begins from 0; this is the axis that steps across image volumes), only coordinate number 0 (i.e. the first volume)

    ::

      ex. (1, 24) --> extract slice number 24 from the axis 1

- *scaling* (a list of items which are a float, optional)
    Specify the data scaling parameters used to rescale the intensity values

    ::

      ex. [0.0, 1.0]

- *vox* (a list of items which are a float or an integer, optional)
    Change the voxel dimensions reported in the output image header

    ::

      ex. [1.25] --> set the voxel size to 1.25 mm for all three spatial axes in the output image.

    ::

      ex. [1, , 3.5] --> change the voxel size along the first and third axes to 1.0mm and 3.5mm respectively, and leave the voxel size along the second axis unchanged

- *suffix* (a string, optional)
    Output file suffix.
    If not provided, for some options a suffix will be added automatically ('scaled' for scaling option, 'vox' for vox option, 'coord' for coord option and 'axes' for axes option)

    ::

      ex. 'extract'

- *bval_scale* (no, yes, optional)
    For diffusion data, enable or disable scaling of diffusion b-values by the square of the corresponding DW gradient norm.
    By default this option is not set and it is choose automatically by mrtrix.

    ::

      ex. 'no'

- *grad_file* (a pathlike object or a string representing an existing file, optional)
    Provide the diffusion-weighted gradient scheme used in the acquisition in a text file (MRTrix format)

    ::

      ex. '/home/username/data/derived_data/diff_mrtrix_format.txt'

- *in_bvec* (a pathlike object or a string representing an existing file, optional)
    Provide Bvecs file in FSL format. If data converted into MIF format, bvec will be added in the output

    ::

      ex. '/home/username/data/raw_data/DWI.bvec'

- *in_bval* (a pathlike object or a string representing an existing file, optional)
    Provide Bvals file in FSL format. If data converted into MIF format, bavl will be added in the output

    ::

      ex. '/home/username/data/raw_data/DWI.bval'


**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    The output image.

    ::

      ex. '/home/username/data/derived_data/DWI.mif'

- *out_bvec* (a pathlike object or string representing a file, optional)
    The out bvec in FSL format

    ::

      ex. '/home/username/data/derived_data/DWI.bvec'

- *out_bval* (a pathlike object or string representing a file, optional)
    The out bvec in FSL format

    ::

      ex. '/home/username/data/derived_data/DWI.bval'

-------------

Usefull links:

`mrtrix mrconvert <https://mrtrix.readthedocs.io/en/latest/reference/commands/mrconvert.html>`_
