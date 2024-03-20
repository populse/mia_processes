:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=======================
NonSteadyDetector brick
=======================

Detect non-steady-state at the beginning of a bold 4D image
-----------------------------------------------------------

Adapted from  `Nipype Confounds <https://github.com/nipy/nipype/blob/f662acfce8def4717e0c3414618f3a5de5913b31/nipype/algorithms/confounds.py#L974>`_.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'


**Outputs parameters:**

- *n_volumes_to_discard* (int)
    Number of volumes to discard.

    ::

      ex. 1

-------------

Usefull links:

`Confunds NonSteadyStateDetector nipype <https://github.com/nipy/nipype/blob/f662acfce8def4717e0c3414618f3a5de5913b31/nipype/algorithms/confounds.py#L974>`_
