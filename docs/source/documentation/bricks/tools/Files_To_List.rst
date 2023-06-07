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

===================
Files_To_List brick
===================

From 3 file names, generating a list containing all theses file names
---------------------------------------------------------------------

- `file1` is mandatory, while `file2` and `file3` are optional. If one (or both) of the latter two entries is (are) undefined, it will not be present in `file_list`.
- Ex. `file2` is undefined and `file3` is defined, the result for file_list will be [file1, file3].

**Inputs parameters:**

- *file1*
    A string corresponding to a path file.

    ::

      ex. /home/ArthurBlair/data/Anat.nii

- *file2*
    An optional string corresponding to a path file.

    ::

      ex.  /home/ArthurBlair/data/Func.nii

**Outputs parameters:**

- *file_list*
    A list.

    ::

      ex. ['/home/ArthurBlair/data/Anat.nii', '/home/ArthurBlair/data/Func.nii']
