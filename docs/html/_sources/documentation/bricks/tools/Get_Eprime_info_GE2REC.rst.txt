:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

============================
Get_Eprime_info_GE2REC brick
============================


Obtain usuful information for GE2REC pipeline from E-Prime file
----------------------------------------------------------------

The E-Prime experimentation is available on request.


**Inputs parameters:**

- *eprime_file* (an existing xlsx file):
    The E-Prime file obtain during the experiment.

    ::

        ex. '/home/username/MIA_projects/data/downloaded_data/001_eprime.xlsx'

**Outputs parameters:**

- *csv_correct_response* (a file):
    The percentage of good responses given by the subject during the RECO run.

    ::

       ex. '/home/username/MIA_projects/data/downloaded_data/001_eprime_correct_response.xlsx'

- *csv_encodage_reco* (a file):
    | Encoding performance during GE run indirectly determine by the responses given during the RECO run:
    | - if the response is correct: encoding set to 2
    | - if the response is bad: encoding set to 1
    | - if it is a control question: encoding set to 1
    | - during rest: encoding set to 0

    ::

       ex. '/home/username/MIA_projects/data/downloaded_data/001_eprime_encodage_reco.xlsx'

- *sess_regress_level1design* (a list of items which are a list of items which are a dictionary with keys which are 'name' or 'val' and with values which are a string or a list of float):
    Encoding performance during GE run (indirectly determine by the responses given during the RECO run)
    to use in the Level1Design brick for sess_regress parameter.

    ::

       ex. [
             [
                {
                    'name': 'ENCODAGE',
                    'val': [0, 1, 2, 2, 2]
                }
            ]
        ]
