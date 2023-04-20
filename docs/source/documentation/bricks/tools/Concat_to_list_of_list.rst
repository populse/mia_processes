
+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

===================
Concat_to_list_of_list brick
===================

Make an output list of list containing the iteration of the input list1 with each element of the input list2

    Ex. ['a', 'b', 'c'] and ['_1', '_2'] gives
    [['a', '_1'], ['a', '_2'],
     ['b', '_1'], ['b', '_2'],
     ['c', '_1'], ['c', '_2']]

---------------------------------------------------------------------

**Inputs parameters:**

- *list1* (a list of string, optional, default value is ['ACA', 'ACM', 'ACP', 'PICA', 'ROI_CING', 'ROI_FRON', 'ROI_INSULA', 'ROI_OCC', 'ROI_PAR', 'ROI_STR', 'ROI_TEMP', 'ROI_THA', 'SCA'])
    List 1.

    ::

      ex. ['ACA', 'ACM', 'ACP', 'PICA', 'ROI_CING', 'ROI_FRON', 'ROI_INSULA', 'ROI_OCC', 'ROI_PAR', 'ROI_STR', 'ROI_TEMP', 'ROI_THA', 'SCA']

- *list2* (a list of string, optional, default value is ['_L', '_R'])
    List 2.

    ::

      ex.  ['_L', '_R']

**Outputs parameters:**

- *listOflist* (a list of list)
    Output list.

    ::

      ex.
