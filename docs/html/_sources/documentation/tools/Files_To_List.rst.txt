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

===================
Files_To_List brick
===================

From 2 file names, generating a list containing all theses file names
---------------------------------------------------------------------

*'/home/ArthurBlair/data/Anat.nii' +  '/home/ArthurBlair/data/Func.nii' -> Files_To_List -> ['/home/ArthurBlair/data/Anat.nii', '/home/ArthurBlair/data/Func.nii']*

>>> from mia_processes.tools import Files_To_List
>>> Files_To_List.help()

**Inputs parameters:**

- *file1*
    A string corresponding to an existing path file.

    ::

      ex. /home/ArthurBlair/data/Anat.nii

- *file2*
    An optional string corresponding to an existing path file.

    ::

      ex.  /home/ArthurBlair/data/Func.nii

**Outputs parameters:**

- *file_list*
    A list.

    ::

      ex. ['/home/ArthurBlair/data/Anat.nii', '/home/ArthurBlair/data/Func.nii']

-------------

NOTE: If only file1 is specified, returns a list containing file1 only.
