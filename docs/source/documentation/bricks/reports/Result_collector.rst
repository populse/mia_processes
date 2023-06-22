:orphan:

.. toctree::

+--------------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+-------------------------------------------+----------------------------------------------------+

======================
Result_collector brick
======================


Generates files containing summary data for each region of interest
-------------------------------------------------------------------

    - To work correctly, the database entry for the first element of `parameter_files` must
      have the `PatientName` tag filled in.
    - The “PatientName_data/results_aggregation” directory is created to receive the results.
      If this directory exists at runtime, new results can overwrite old results with the same name.
    - Currently, to work correctly, this brick requires the doublet made up of the two hemispheres
      to be present in the parameter_files list and each hemisphere to be represented by the letters
      L (left) and R (right).

      | For example:
      | [/aPath/ACM_R_mean_spmT_BOLD.txt, /aPat/ACM_L_mean_spmT_BOLD.txt, etc.].

      It would be desirable to develop this brick so that it could also be used to collect a single
      territory without any notion of hemisphere (in this case, of course, the brick would not generate
      any laterality indices) => TODO ASAP

--------------------------------------

**Inputs parameters:**

- *parameter_files*
    A list of files, each containing a parameter value. To work correctly, the name of each file must be exactly like this:
        - ``roi``\_ ``hemi``\_ ``calcul``\_ ``param``\_ ``contrast``.txt, where
            - ``roi``: region of interest (ex. ACA)
            - ``hemi``: hemisphere (ex. L)
            - ``calcul``: type of calcul (ex. mean)
            - ``param``: the parameter recorded in the file (ex. spmT)
            - ``contrast``: the type of contrast/effect used (ex. BOLD)

    ::

     ex. ['/home/username/data/raw_data/ACA_L_mean_spmT_BOLD.txt',
          '/home/username/data/raw_data/ACA_R_mean_spmT_BOLD.txt',
	  '/home/username/data/raw_data/ACM_L_mean_spmT_BOLD.txt',
	  '/home/username/data/raw_data/ACM_R_mean_spmT_BOLD.txt',




- *laterality_index*
    | A Boolean to calculate (True) or not (False) the laterality index:
    | (left hemisphere parameter - right hemisphere parameter) / (left hemisphere parameter + right hemisphere parameter)

    ::

      ex. True

- *patient_info*
    A dictionary whose keys/values correspond to information about the patient.

    e.g. {
          | 'PatientName': 'ablair',
          | 'Pathology': 'ACMD',
          | 'Age': 64,
	  | 'Sex': 'M',
	  | 'MR': '3T',
	  | 'Gas': 'BACTAL',
	  | 'GasAdmin': 'MASK'
	 }

    ::

      ex. {'PatientName': <undefined>, 'Pathology': <undefined>, 'Age': <undefined>, 'Sex': <undefined>, 'MR': <undefined>, 'Gas': <undefined>, 'GasAdmin': <undefined>}


**Outputs parameters:**

- *out_files*
    | A list of .xml files containing a summary of the input parameters. The file names generated are constructed as follows:
    | ``contrast``\_ ``calcul``\_ ``param``.txt (e.g. BOLD_std_beta.xls).

    ::

      ex. ['/home/username/data/derived_data/patient-name_data/results_aggregation/BOLD_IL_mean_spmT.xls',
           '/home/username/data/derived_data/patient-name_data/results_aggregation/BOLD_mean_spmT.xls']
