:orphan:

.. toctree::

+-----------------------------+----------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+----------------------------------------+----------------------------------------------------+

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

==================
List_To_File brick
==================

From a list of files, generation of a single file
-------------------------------------------------

*['/home/ArthurBlair/data/Anat.nii', '/home/ArthurBlair/data/Func.nii'] -> List_To_File -> '/home/ArthurBlair/data/Func.nii'*

>>> from mia_processes.tools import List_to_File
>>> List_to_File.help()

**Inputs parameters:**

- *file_list*
    The list of elements to be filtered.

    ::

      ex. ['/home/ArthurBlair/data/Anat.nii', '/home/ArthurBlair/data/Func.nii']

- *index_filter*
    A list of 0 to 1 indexes (integer) for filtering.

    ::

      ex.  [1]

**Outputs parameters:**

- *file*
    The corresponding filtering result (a file).

    ::

      ex. '/home/ArthurBlair/data/Func.nii'

-------------

NOTE: If no index is specified, returns the first element of the list.
