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
Input_Filter brick
==================

To filter the content of the Data Browser tab or the output data of another brick
---------------------------------------------------------------------------------

*[Preselected data (or not) from the Data Browser or the output data from another brick] -> Input_Filter -> ['/home/ArthurBlair/data/raw_data/Anat.nii']*

>>> from mia_processes.tools import Input_Filter
>>> Input_Filter.help()

**Inputs parameters:**

- *input*
    A list corresponding to the data from the Data Browser or the output data from another brick.

    ::

      ex. ['/home/ArthurBlair/data/Anat.nii', '/home/ArthurBlair/data/Func.nii']

**Outputs parameters:**

- *output*
    A list with the result of the filter applied.

    ::

      ex. ['/home/ArthurBlair/data/Func.nii']

-------------

NOTE: To run properly, this process (node, brick) needs that, on the one hand, in the DataBrowser, the "send documents to the Pipeline Manager" button was clicked previously (with a selection of data made, or not), and on the other hand, at the input_filter brick level, that a right click then the selection of the option "Export to database_scans" has been made. Finally, a right click on the input_filter brick will allow to filter the input data by selecting "open filter".
