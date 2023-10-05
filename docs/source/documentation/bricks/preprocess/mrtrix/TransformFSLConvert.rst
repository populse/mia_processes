:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=========================
TransformFSLConvert brick
=========================

Perform conversion between FSL’s transformation matrix format to mrtrix3’s
---------------------------------------------------------------------------

(mrtrix transformconvert command)

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a pathlike object or a string representing an existing file)
    FLIRT input image (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/DWI.mif'

- *reference* (a pathlike object or a string representing an existing file)
    FLIRT reference image (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/T1w.mif'

- *in_transform* (a pathlike object or a string representing an existing file)
    FLIRT output transformation matrix (valid extensions: [.mif, .nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/diff2struct_fsl.mat'


**Outputs parameters:**

- *out_transform* (a pathlike object or string representing a file)
    Output CSF response text file

    ::

      ex. '/home/username/data/derived_data/diff2struct_fsl.txt'

-------------

Usefull links:

`mrtrix transformconvert <https://mrtrix.readthedocs.io/en/latest/reference/commands/transformconvert.html>`_

`FSL FLIRT <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT>`_
