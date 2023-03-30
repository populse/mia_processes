:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
BetSurfacesExtraction brick
============

 Surfaces (skull, inskull, outskull, outskin) extraction using BET (FSL). 

 Both bet2 and betsurf programs are used in order to get skull and scalp
 surfaces created by betsurf (fsl BET -A option)
 This involves registering to standard space in order to allow betsurf
 to find the standard space masks it needs.

 The mask and mesh files (.vtk) are generated. 

>>> from mia_processes.bricks.preprocess.fsl import BetSurfacesExtraction
>>> BetSurfacesExtraction.help()

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input file to skull strip. An existing file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'

**Optional inputs with default value parameters:**

- *output_type* ('NIFTI' or 'NIFTI_GZ', optional)
    | Format of the output image (one of NIFTI, NIFTI_GZ).
    |   NIFTI: \*.nii
    |   NIFTI_GZ: \*.nii.gz

    ::

      default value. NIFTI


**Outputs parameters:**

- *inskull_mask_file* (a pathlike object or string representing a file)
    Inskull mask (extensions: [.nii, .nii.gz])

    ::

      ex. '/home/username/data/raw_data/T1w_brain_inskull_mask.nii'

- *inskull_mesh_file* (a pathlike object or string representing a file)
    Inskull mesh  (extensions [.vtk])

    ::

      ex. '/home/username/data/raw_data/T1w_brain_inskull_mesh.vtk'

- *out_file* (a pathlike object or string representing a file)
    Skullstripped file (extensions: [.nii, .nii.gz])

    ::

      ex. /home/username/data/raw_data/T1w_brain.nii

- *outskin_mask_file* (a pathlike object or string representing a file)
    Outskin mask (extensions: [.nii, .nii.gz])

    ::

      ex. '/home/username/data/raw_data/T1w_brain_outskin_mask.nii'

- *outskin_mesh_file* (a pathlike object or string representing a file)
    Outskin mesh  (extensions [.vtk])

    ::

      ex. '/home/username/data/raw_data/T1w_brain_outskin_mesh.vtk'

- *outskull_mask_file* (a pathlike object or string representing a file)
    Outskull mask (extensions: [.nii, .nii.gz])

    ::

      ex. '/home/username/data/raw_data/T1w_brain_outskull_mask.nii'

- *outskull_mesh_file* (a pathlike object or string representing a file)
    Outskull mesh (extensions [.vtk])

    ::

      ex. /home/username/data/raw_data/T1w_brain_outskull_mesh.vtk

- *skull_mask_file* (a pathlike object or string representing a file)
    Skull mask (extensions: [.nii, .nii.gz])

    ::

      ex. '/home/username/data/raw_data/T1w_brain_skull_mask.nii'
      

-------------

Usefull links:
`FSL BET <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BET/UserGuide>`_
`FSL FAST - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.fsl.preprocess.html#bet>`_