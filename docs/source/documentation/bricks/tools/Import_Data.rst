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

==================
Import_Data  brick
==================

Import reference data into the current pipeline
-----------------------------------------------

- This brick was originally written to automatically select and load reference data during pipeline operation.

- The `rois_list` parameter is used to filter the data to be imported from a library defined by `lib_dir`.
  If `roi_list` is a list, each element of it will be a filename filter applied to select reference data.
  If `roi_list` is a list of lists, the filters will result from concatenating the elements of each internal
  list with the underscore caracter (e.g. [["foo", "1"], ["faa", "2"]] gives two filters, "foo_1" and "faa_2".

- If `lib_dir` is not set, the default library used will be the miaresources/ROIs/ directory.

- The `file_in_db` file is used only to retrieve the value of the associated `PatientName` tag.

- The reference data is imported into the `output_directory/PatientName_data/ROI_data/raw_data` directory.


**Inputs parameters:**

- *rois_list*
    A list or list of lists of strings, defining the data to be imported (used as a filter).

    ::

      ex. [['ACA', 'L'], ['ACA', 'R'], ['ACM', 'L'], ['ACM', 'R']]

- *lib_dir*
    The path to a data library (if not defined, the default resources path is used, i.e. miaresources/ROIs/)."

    ::

      ex. <undefined>

- *file_in_db*
    A file in database, only used to catch the `PatientName` tag value.

    ::

      ex. /home/username/data/Anat.nii

- *starts_with*
    If True, applies the file filter only to the beginning of the names, otherwise to the whole names (a boolean).

   ::

      ex. True


**Outputs parameters:**

- *rois_files*
    The list of resulting available files (a list of files).

    ::

      ex. ['/home/username/data/derived_data/patient-name_data/ROI_data/raw_data/ACA_L.nii',
           '/home/username/data/derived_data/patient-name_data/ROI_data/raw_data/ACA_R.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/raw_data/ACM_L.nii',
	   '/home/username/data/derived_data/patient-name_data/ROI_data/raw_data/ACM_R.nii']
