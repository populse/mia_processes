:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

===============
GetLabels brick
===============

Get Assemblynet segmentation labels
------------------------------------

Choose between "tissues", "lobes", "structures", "macrostructures" and obtain the corresponding labels and the label names

--------------------------------------

**Optional inputs parameters:**

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

- *labels* (a list of integer)
    List of the corresponding labels

    ::

      ex. [1, 2, 3, 4, 5, 6, 7]

- *names* (a list of string)
     List of the corresponding label names.

    ::

      ex. ['CSF', 'Cortical', 'Cerebrum WM', 'Subcortocal GM', 'Cerebellum GM', 'Cerebellum WM', 'Brainstem']


-------------

Usefull links:

`volBrain Assemblynet <https://github.com/volBrain/AssemblyNet>`_

`brainCOLOR protocol <https://mindboggle.info/braincolor/>`_

