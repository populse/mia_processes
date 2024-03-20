:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==========================
LabelsCorrespondence brick
==========================

Get AssemblyNet labels names or get labels from names
-----------------------------------------------------

In each segmentation file obtain (tissues, lobes...) with AssemblyNet, each ROI is defined by a label (an integer).
This brick is used to obtain the correspondence between the label and the label name, depending on the type of the segmentation file.

Add a list of lables or a list if label names.
Choose between "tissues", "lobes", "structures", "macrostructures" and obtain the corresponding labels or label names.

--------------------------------------

**Optional inputs parameters:**

- *labels_names* (a list of integer or a list of string, optional)
    List of labels or names  for which the corresponding name / label is wanted

    ::

      ex. [1, 2, 3, 4, 5, 6, 7]

    ::

      ex2. ['CSF', 'Cortical', 'Cerebrum WM', 'Subcortocal GM', 'Cerebellum GM', 'Cerebellum WM', 'Brainstem']

- *lobes* (a boolean, optional)
    Get labels for lobes segmentation

    ::

      ex. False

- *macrostructures* (a boolean, optional)
    Get labels for macrostructures segmentation

    ::

      ex. False

- *structures* (a boolean, optional)
    Get labels for structures segmentation

    ::

      ex. False

- *tissues* (a boolean, optional)
    Get labels for tissues segmentation

    ::

      ex. True


**Outputs parameters:**

- *correspondence* (a list of integer or a list of string)
    List of the corresponding labels

    ::

      ex.  ['CSF', 'Cortical', 'Cerebrum WM', 'Subcortocal GM', 'Cerebellum GM', 'Cerebellum WM', 'Brainstem']

    ::

      ex2. [1, 2, 3, 4, 5, 6, 7]



-------------

Usefull links:

`volBrain Assemblynet <https://github.com/volBrain/AssemblyNet>`_

`brainCOLOR protocol <https://mindboggle.info/braincolor/>`_
