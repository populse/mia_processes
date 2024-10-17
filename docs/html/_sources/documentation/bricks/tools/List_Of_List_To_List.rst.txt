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

==========================
List_Of_List_To_List brick
==========================

From a list of list, generation of a single list
-------------------------------------------------

**Inputs parameters:**

- *list_of_list (a list of list)*
    The list of list to be filtered.

    ::

      ex. [['/home/username/data/Anat.nii'], ['/home/username/data/Func.nii']]

- *index_filter (a list og integer, default value is [1])*
    A list of 0 to 1 indexes (integer) for filtering.

    ::

      ex.  [1]

**Outputs parameters:**

- *list (a list)*
    The corresponding filtering result (a list).

    ::

      ex. ['/home/username/data/Anat.nii']

-------------

NOTE: If no index is specified, returns the first element of the list (the index corresponds to the actual position in the list, i.e. there is no “0” index as in python).
