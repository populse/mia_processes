:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=======================
AssemblyNetDocker brick
=======================

3D Whole Brain MRI Segmentation using AssemblyNet (volBrain)
-------------------------------------------------------------

This brick allows to run the AssemblyNet Docker image.

Based on a large ensemble of convolution neural networks, AssemblyNet segment a T1w image in 133 labels (according to BrainColor protocol).

| If you use this brick:
|   - cite: Pierrick Coupé, Boris Mansencal, Michaël Clément, Rémi Giraud, Baudouin Denis de Senneville, Vinh-Thong Ta, Vincent Lepetit, José V. Manjon `AssemblyNet <https://www.sciencedirect.com/science/article/pii/S1053811920305127?via%3Dihub>`_: A large ensemble of CNNs for 3D whole brain MRI segmentation. NeuroImage, Elsevier, 2020, 219, pp.117026.
|   - check the AssemblyNet's `license <https://github.com/volBrain/AssemblyNet/blob/main/README.md#license>`_


*This brick requires Docker*

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input T1w image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'


**Outputs parameters:**

- *native_t1* (a pathlike object or string representing a file)
    Filtered and normalized T1 image in native space.

    ::

      ex. '/home/username/data/derived_data/native_t1_T1w.nii.gz'

- *mni_t1* (a pathlike object or string representing a file)
    Filtered and normalized T1 image in MNI space.

    ::

      ex. '/home/username/data/derived_data/mni_t1_T1w.nii.gz'

- *native_mask* (a pathlike object or string representing a file)
    Intracranial Cavity mask image in native space.

    ::

      ex. '/home/username/data/derived_data/native_mask_T1w.nii.gz'

- *mni_mask* (a pathlike object or string representing a file)
    Intracranial Cavity mask image in MNI space.

    ::

      ex. '/home/username/data/derived_data/mni_masks_T1w.nii.gz'

- *native_structures* (a pathlike object or string representing a file)
    Structures segmentation in native space.

    ::

      ex. '/home/username/data/derived_data/native_structures_T1w.nii.gz'

- *mni_structures* (a pathlike object or string representing a file)
    Structures segmentation in MNI space.

    ::

      ex. '/home/username/data/derived_data/mni_structures_T1w.nii.gz'

- *native_macrostructures* (a pathlike object or string representing a file)
    Macrostructures segmentation in native space.

    ::

      ex. '/home/username/data/derived_data/native_macrostructures_T1w.nii.gz'

- *mni_macrostructures* (a pathlike object or string representing a file)
    Macrostructures segmentation in MNI space.

    ::

      ex. '/home/username/data/derived_data/mni_macrostructures_T1w.nii.gz'

- *native_lobes* (a pathlike object or string representing a file)
    Lobes segmentation in native space.

    ::

      ex. '/home/username/data/derived_data/native_lobes_T1w.nii.gz'

- *mni_lobes* (a pathlike object or string representing a file)
    Lobes segmentation in MNI space.

    ::

      ex. '/home/username/data/derived_data/mni_lobes_T1w.nii.gz'

- *native_tissues* (a pathlike object or string representing a file)
    Tissues segmentation in native space.

    ::

      ex. '/home/username/data/derived_data/native_tissues_T1w.nii.gz'

- *mni_tissues* (a pathlike object or string representing a file)
    Tissues segmentation in MNI space.

    ::

      ex. '/home/username/data/derived_data/mni_tissues_T1w.nii.gz'

- *matrix_affine* (a pathlike object or string representing a file)
    ITK transformation matrix from native to MNI space

    ::

      ex. '/home/username/data/derived_data/matrix_affine_native_to_mni_T1w.txt'

- *report_pdf* (a pathlike object or string representing a file)
    Volumetry report (PDF format)

    ::

      ex. '/home/username/data/derived_data/report_T1w.pdf'

- *report_csv* (a pathlike object or string representing a file)
    Volumetry report (csv format)

    ::

      ex. '/home/username/data/derived_data/report_T1w.csv'


-------------

Usefull links:

`volBrain Assemblynet <https://github.com/volBrain/AssemblyNet>`_

`brainCOLOR protocol <https://mindboggle.info/braincolor/>`_

`Docker <https://docs.docker.com/get-docker/>`_
