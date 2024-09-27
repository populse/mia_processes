:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

=======================
Get_Patient_Name brick
=======================

Get patient name from a file. "PatientName" tag should be filled for this file in the database.

-----------------------------------------------


**Inputs parameters:**

- *in_file* (an existing file):
    An existing file.

    ::

      ex. '/home/username/MIA_projects/data/derived_data/T1w.nii'


**Outputs parameters:**

- *patient_name* (a string):
    Patient name tag as written in the database

    ::

      ex. '001'
