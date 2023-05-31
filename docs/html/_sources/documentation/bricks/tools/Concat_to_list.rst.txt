:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

=====================
Concat_to_list  brick
=====================

Make an output list corresponding to the concatenation of list1 and list2

    Ex. ['a', 'b', 'c'] and ['d', 'e'] gives:
        ['a', 'b', 'c', d', 'e']

---------------------------------------------------------------------

**Inputs parameters:**

- *list1* (a list)

    ::

      ex. ['ACA', 'ACM']

- *list2* (a list)

    ::

      ex.  ['L', 'R']

**Outputs parameters:**

- *out_list"*
    A list corresponding to the concatenation of list1 and list2.

    ::

      ex. ['ACA', 'ACM', 'L', 'R']
