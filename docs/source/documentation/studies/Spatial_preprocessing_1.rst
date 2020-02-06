:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

.. line break
.. |br| raw:: html

   <br />

.. thin space
.. |ws1| raw:: html

   &thinsp;

.. em space

.. |ws2| raw:: html

   &emsp;

.. en space

.. |ws3| raw:: html

   &ensp;

.. non-breakable space

.. |ws4| raw:: html

   &nbsp;

================================
Spatial_preprocessing_1 pipeline
================================

Data pre-processing for cerebrovascular reserve analysis (CVRa) at `CLUNI <http://www.neuroradiologie-grenoble.fr/>`_ - `IRMaGe <https://irmage.univ-grenoble-alpes.fr/>`_ (Grenoble - France)
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

**Pipeline insight**

- Anatomical image: NewSegment -> Normalize12
- Functional images: Realign -> Coregister (to anatomical image) -> Normalize12 -> Smooth

**Inputs parameters**

- *anat_file*
    An anatomical image (ex. 3D T1 sequence sush as T1 turbo field echo). An existing, uncompressed file (valid extensions: [.img, .nii,
    .hdr]).

    ::

      ex. /home/ArthurBlair/data/raw_data/Anat.nii

- *func_files*
    Functional images under hypercapnic challenge (ex. 3D T2* sequence sush as echo planar imaging). A list of items which are an
    existing, uncompressed file (valid extensions: [.img, .nii, .hdr]).

    ::

      ex. ['/home/ArthurBlair/data/raw_data/Func.nii']

- *voxel_sizes_func* [#label]_
    The voxel sizes (x, y & z, in mm) of the written normalised functional images (this the input write_voxel_sizes parameter of the
    Normalize12 brick for the functional images). A list of 3 items which are a float.

    ::

      ex. [2.0, 2.0, 2.0]

**Outputs parameters:**

- *forward_deformation_field*
    Forward deformation field. Could be used for spatially normalising images to MNI space (a list of items which are a pathlike object or
    string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/raw_data/y_Anat.nii

- *bias_field_images*
    The estimated bias field (a list of items which are a pathlike object or string representing an existing file).

    ::

      ex. <undefined>


- *bias_corrected_images*
    The bias corrected images (a list of items which are a pathlike object or string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/raw_data/mAnat.nii

- *native_class_images*
    Native space probability maps (a list of items which are a list of items which are a pathlike object or string representing an existing
    file).

    ::

      ex. [['/home/ArthurBlair/data/raw_data/c1Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/c2Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/c3Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/c4Anat.nii'],
           ['/home/ArthurBlair/data/raw_data/c5Anat.nii']]

- *realignment_parameters*
    The estimated translation and rotation parameters during the realign stage (a list of items which are a pathlike object or string
    representing an existing file).

    ::

      ex. /home/ArthurBlair/data/raw_data/rp_Func.txt

- *normalized_anat*
    The final normalised anatomical image (a list of items which are a pathlike object or string representing an existing file).

    ::

      ex. /home/ArthurBlair/data/raw_data/wAnat.nii

- *smoothed_func*
    The final, realigned then coregistered then normalised then smoothed, functional images (a list of items which are an existing file
    name).

    ::

      ex. /home/ArthurBlair/data/raw_data/swrFunc.nii

- *coregistered_source*
    Coregistered source files, corresponding to ‘source’ images (a list of items which are an existing file name).

    ::

      ex. /home/ArthurBlair/data/raw_data/meanFunc.nii

-------------

.. [#label] Depends of the study; the value given as an example is valid for cevastoc, cevastoc32, etc. CVRa studies at
	    `CLUNI <http://www.neuroradiologie-grenoble.fr/>`_ - `IRMaGe <https://irmage.univ-grenoble-alpes.fr/>`_.
