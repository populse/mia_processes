:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

=====================
Deconv_from_aif brick
=====================

Deconvolution of the tissue response curve with AIF
---------------------------------------------------

- MRI perfusion imaging by deconvolution using an Arterial Input Function
  (AIF) is a method used to evaluate blood flow and perfusion in tissues,
  often in the brain, to diagnose and assess conditions like stroke, tumors,
  or other vascular abnormalities. This method is part of a broader category
  of dynamic susceptibility contrast (DSC) MRI techniques.

- DSC-MRI involves the injection of a contrast agent (usually
  gadolinium-based) and the rapid acquisition of MRI images as the contrast
  passes through the blood vessels. This generates time-series data that
  reflect how the contrast agent is distributed in the tissue over time.

- AIF represents the concentration of contrast agent over time in a feeding
  artery. It is used as a reference to determine the blood supply to
  the tissue.

- Deconvolution is used to separate the tissue response from the AIF, thereby
  allowing the calculation of perfusion parameters like cerebral blood flow
  (CBF), cerebral blood volume (CBV), mean transit time (MTT), time to peak
  (TTP), time to maximum (Tmax) and bolus arrival time (T0).

-----------------------------------------------

**Inputs parameters:**

*func_file*
    T2* functional Magnetic Resonance Imaging (fMRI) experiment recorded
    during gadolinium bolus. Ideally, the data will have been pre-processed
    (realignment, segmentation, etc.). An existing, uncompressed file
    (valid extensions: [.img, .nii, .hdr]).

    ::

      ex. '/home/username/data/raw_data/swrfunc.nii'

*aif_file*
    The AIF (a file in .json format) typically obtained previously with the
    Make_AIF brick.

    ::

      ex. '/home/username/data/raw_data/swrfunc_aif.json'

*mask_file*
    The mask (at the resolution of `func_file`) used for the perfusion
    deconvolution (valid extensions: .nii).

    ::

      ex. '/home/username/data/raw_data/mask_swc1_anat_003.nii'

**Outputs parameters:**

- *CBV_image*
    Cerebral Blood Volume (a file with .nii format) is a parameter that
    measures the total volume of blood present within a given volume of brain
    tissue, in milliliters of blood per 100 grams of brain tissue (mL/100g).
    CBV provides important information about the vascularity and perfusion
    characteristics of brain tissue.
    ::

      ex. '/home/username/data/derived_data/swrfunc_CBV_deconv.nii'

- *CBF_image*
    Cerebral Blood Flow (a file with .nii format) is a parameter that measures
    the rate at which blood is delivered to the brain tissue, in milliliters
    of blood per 100 grams of brain tissue per minute (mL/100g/min). CBF
    provides valuable information about the brain's blood supply and is used
    to assess the adequacy of cerebral perfusion in various neurological
    conditions.

    ::

      ex. '/home/username/data/derived_data/swrfunc_CBF_deconv.nii'

- *MTT_image*
    Mean Transit Time (a file with .nii format) represents the average time
    (s) it takes for blood to pass through a given region of tissue. It is an
    important indicator of the efficiency of blood flow and is used alongside
    other perfusion metrics like CBF and CBV to assess the health of brain
    tissue.

    ::

      ex. '/home/username/data/derived_data/swrfunc_MTT_deconv.nii'

- *TTP_image*
    Time to Peak (a file with .nii format) reflects the time (s) it takes for
    the contrast agent to reach its maximum concentration in a given voxel
    after its arrival. This metric provides insight into the dynamics of blood
    flow and is particularly useful in evaluating conditions like stroke,
    tumors, and other cerebrovascular disorders.
    ::

      ex. '/home/username/data/derived_data/swrfunc_TTP_deconv.nii'

- *Tmax_image*
    The Time to Maximum (a file with .nii format) represents the time delay
    (s) between the arrival of contrast agent in the arterial input function
    (AIF) and the peak of the residue function (the tissue response after
    deconvolution).

    ::

      ex. '/home/username/data/derived_data/swrfunc_Tmax_deconv.nii'

- *T0_image*
    The Bolus Arrival Time (a file with .nii format) represents the time (s)
    at which the contrast agent first arrives at a particular voxel in the
    tissue. It essentially marks the onset of contrast passage through each
    voxel.

    ::

      ex. '/home/username/data/derived_data/swrfunc_T0_deconv.nii'
