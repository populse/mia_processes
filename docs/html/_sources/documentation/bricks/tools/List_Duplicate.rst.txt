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

====================
List_Duplicate brick
====================

From a file name, generating a list containing this file name and the file name itself
---------------------------------------------------------------------------------------

**Inputs parameters:**

- *file_name*
    A string corresponding to an existing path file.

    ::

      ex. /home/username/data/Func.nii

**Outputs parameters:**

- *out_file*
    A string corresponding to an existing path file.

    ::

      ex.  /home/username/data/Func.nii

- *out_list*
    A list with one string element corresponding to an existing path file.

    ::

      ex. ['/home/username/data/Func.nii']
