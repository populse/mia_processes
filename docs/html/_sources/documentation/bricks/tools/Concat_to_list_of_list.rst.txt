:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

============================
Concat_to_list_of_list brick
============================

Make an output list of list containing the iteration of the input list1 with each element of the input list2
------------------------------------------------------------------------------------------------------------

    Ex. ['a', 'b', 'c'] and ['1', '2'] gives:
        [['a', '1'], ['a', '2'], ['b', '1'], ['b', '2'], ['c', '1'], ['c', '2']]

---------------------------------------------------------------------

**Inputs parameters:**

- *list1*
    A list of string, optional. Default value is ['ACA', 'ACM', 'ACP', 'PICA', 'ROI_CING', 'ROI_FRON', 'ROI_INSULA', 'ROI_OCC', 'ROI_PAR', 'ROI_STR', 'ROI_TEMP', 'ROI_THA', 'SCA'].

    ::

      ex. ['ACA', 'ACM', 'ACP', 'PICA', 'ROI_CING', 'ROI_FRON', 'ROI_INSULA', 'ROI_OCC', 'ROI_PAR', 'ROI_STR', 'ROI_TEMP', 'ROI_THA', 'SCA']

- *list2*
    A list of string, optional. Default value is ['L', 'R']).

    ::

      ex. ['L', 'R']

**Outputs parameters:**

- *listOflist*
    Output list of lists.

    ::

      ex. [['ACA', 'L'] , ['ACA', 'R'],
           ['ACM', 'L'], ['ACM', 'R'],
           ['ACP', 'L'] , ['ACP', 'R'],
	   ['PICA', 'L'], ['PICA', 'R'],
	   ['ROI_CING', 'L'], ['ROI_CING', 'R'],
	   ['ROI_FRON', 'L'], ['ROI_FRON', 'R'],
	   ['ROI_INSULA', 'L'], ['ROI_INSULA', 'R'],
	   ['ROI_OCC', 'L'], ['ROI_OCC', 'R'],
	   ['ROI_PAR', 'L'], ['ROI_PAR', 'R'],
	   ['ROI_STR', 'L'], ['ROI_STR', 'R'],
	   ['ROI_TEMP', 'L'], ['ROI_TEMP', 'R'],
	   ['ROI_THA', 'L'], ['ROI_THA', 'R'],
	   ['SCA', 'L'], ['SCA', 'R']
	  ]
