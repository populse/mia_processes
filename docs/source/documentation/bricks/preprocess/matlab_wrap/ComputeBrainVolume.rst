:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

========================
ComputeBrainVolume brick
========================

Compute brain volume using an input image.

*This brick requires an exclusive Matlab license*

**Disclaimer**: This brick is provided as a proof of concept for developers.
It shows how to wrap a MATLAB script in a mia processes brick using nipype (and so with a Matlab license) in the absence of a better solution (with MCR).

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input file. An existing file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'


**Outputs parameters:**

- *volume* (an integer)
   Brain volume

    ::

      ex. 4929083

-------------

Usefull links:

`Nipype Matlab <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.matlab.html>`_
