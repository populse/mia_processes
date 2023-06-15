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

=======================
Filter_Files_List brick
=======================

Selects one or more (slicing) elements from a list
--------------------------------------------------

**Inputs parameters:**

- *in_list*
    The list of elements to be filtered.

    ::

      ex. ['/home/username/data/raw_data/Anat.nii',
           '/home/username/data/raw_data/Func.nii',
           '/home/username/data/raw_data/rp_Func.txt',
           '/home/username/data/raw_data/meanFunc.nii']

- *index_filter*
    A list of 0 index to 2 indices (integer) for filtering. If no index is given, takes the first element of in_list. If there is one index, takes the element at that position in in_list. If two indices, selects the elements between the first and the second index (the index corresponds to the actual position in the list, i.e. there is no "0" index as in python).

    ::

      ex. [2, 4]

**Outputs parameters:**

- *filtered_list*
    The corresponding filtering result (a list).

    ::

      ex. ['/home/username/data/raw_data/Func.nii',
           '/home/username/data/raw_data/rp_Func.txt',
           '/home/username/data/raw_data/meanFunc.nii']
