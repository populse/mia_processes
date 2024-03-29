:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

=================
Delete_data brick
=================

Delete data from database
-------------------------

- The input of this brick should be an output file from a brick or a pipeline.
  All the outputs from this file history will be removed.

- If "to_keep_filters" is used, the files matching the regex of the filter will be kept.

- If "to_remove_filter" is used, the files matching the regex of the filter will be deleted.

- You can check your regex `here <https://regex101.com/>`_.

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Output file from a brick or a pipeline. The outputs from this file history will be removed.
    Be carrefull, this file will be also deleted, except if your add a regex matching this file in "to_keep_filters" or if you are using "to_remove_file".

    ::

      ex. '/home/username/data/derived_data/sub-01_anatomical_mriqcReport_2023_05_12_16_26_12_58.pdf'

**Optional inputs with default value parameters:**

- *to_keep_filters* (a list of regex, optional, default value is ["(.)*pdf", "(.)*_qc.json", "(.)*desc-carpet(.)*"])
    A list of regex.
    Files that match those regex will be kept and the others files will be deleted. Mutually exclusif with to_remove_filters.

    ::

      ex. ["(.)*pdf", "(.)*_qc.json", "(.)*desc-carpet(.)*"]

- *to_remove_filters* (a list of regex, optional, default value is [])
    A list of regex.
    Files that match those regex will be deleted and the others files will be kept. Mutually exclusif with to_remove_filters.

    ::

      ex. ["(.)*n4c(.)*"]


**Outputs parameters:**


- *files_removed* (a list of string that represents a file)
    List of the file removed by the brick.
    ::

      ex. ['/home/username/data/derived_data/automask_func.nii']
