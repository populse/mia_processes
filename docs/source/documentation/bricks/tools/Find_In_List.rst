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
Find_In_List brick
===================

From a list of files, select the 1rst element that contains a pattern
---------------------------------------------------------------------


**Inputs parameters:**

- *in_list*
    List of files to parse. A list of items which are an existing, uncompressed file (valid extensions: [.img, .nii, .hdr]).

    ::

      ex. ['/home/username/data/Anat1_001.nii', '/home/username/data/Anat1_002.nii', '/home/username/data/Anat1_003.nii', '/home/username/data/Anat2_002.nii']

- *pattern*
    Pattern to look for (a string).

    ::

      ex. 002

**Outputs parameters:**

- *out_file*
    The first file found in the `in_list` parameter with the `pattern` in its name (a file).

    ::

      ex. /home/username/data/Anat1_002.nii
