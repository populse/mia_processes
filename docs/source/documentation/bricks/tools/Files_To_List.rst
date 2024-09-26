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

From 3 file names, generate a list containing all these file names
------------------------------------------------------------------

- `file1` is mandatory, while `file2` and `file3` are optional. If one
  (or both) of the latter two entries is (are) undefined, it will not be
  present in `file_list`.
- Ex. `file2` is undefined and `file3` is defined, the result for file_list
  will be [file1, file3].

**Inputs parameters:**

- *file1*
    A mandatory string corresponding to a path file.

    ::

      ex. /home/username/data/raw_data/Anat.nii

- *file2*
    An optional string corresponding to a path file.

    ::

      ex. /home/username/data/raw_data/Func.nii

- *file3*
    An optional string corresponding to a path file.

    ::

      ex. <undefined>

**Outputs parameters:**

- *file_list*
    A list of files.

    ::

      ex. ['/home/username/data/Anat.nii', '/home/username/data/Func.nii']
