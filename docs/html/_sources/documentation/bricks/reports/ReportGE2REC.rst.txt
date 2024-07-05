:orphan:

.. toctree::

+--------------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+-------------------------------------------+----------------------------------------------------+

==================
ReportGE2REC brick
==================

Generates the report for GE2REC pipeline (langage and memory)
--------------------------------------------------------------


**Inputs parameters:**

- *norm_anat* (a string representing an existing file)
    Normalised anatomical image (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/wAnat.nii'

- *norm_func_gene* (a string representing an existing file)
    Normalised functional image for the generation task (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/wrFunc_gene.nii'

- *norm_func_reco* (a string representing an existing file)
    Normalised functional image for the recognition task (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/wrFunc_reco.nii'

- *norm_func_recall* (a string representing an existing file)
    Normalised functional image for the recall task (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/derived_data/wrFunc_recall.nii'


- *realignment_parameters_gene* (a string representing an existing file)
    Estimated translation and rotation parameters during fMRI recording
    for the generation task
    (valid extensions: .txt).

    ::

      ex. /home/username/data/derived_data/rp_Func_gene.txt

- *realignment_parameters_reco* (a string representing an existing file)
    Estimated translation and rotation parameters during fMRI recording
    for the recognition task
    (valid extensions: .txt).

    ::

      ex. /home/username/data/derived_data/rp_Func_reco.txt

- *realignment_parameters_recall* (a string representing an existing file)
    Estimated translation and rotation parameters during fMRI recording
    for the generation task
    (valid extensions: .txt).

    ::

      ex. /home/username/data/derived_data/rp_Func_recall.txt

- *spmT_gene* (a string representing an existing file)
    A file containing t-statistic values per voxel, which indicate the
    strength of the effect of interest at each brain voxel, derived from the
    general linear model (GLM) analysis performed in SPM for generation task
    (valid extensions:.nii)

    ::

      ex. /home/username/data/derived_dat/patient/stats_gene/spmT_0001.nii

- *spmT_reco* (a string representing an existing file)
    A file containing t-statistic values per voxel, which indicate the
    strength of the effect of interest at each brain voxel, derived from the
    general linear model (GLM) analysis performed in SPM for recognition task
    (valid extensions:.nii)

    ::

      ex. /home/username/data/derived_dat/patient/stats_reco/spmT_0001.nii

- *spmT_recall* (a string representing an existing file)
    A file containing t-statistic values per voxel, which indicate the
    strength of the effect of interest at each brain voxel, derived from the
    general linear model (GLM) analysis performed in SPM for recall task
    (valid extensions:.nii)

    ::

      ex. /home/username/data/derived_dat/patient/stats_recall/spmT_0001.nii

- *li_curves* (a list of string representing an existing file)
    Images (.png) of the lateralization index curve.

    ::

      ex. [/home/username/data/derived_dat/patient/stats_gene/spmT_0002_LI_frontal.png,
      /home/username/data/derived_dat/patient/stats_gene/spmT_0002_LI_temporal.png]

- *correct_answer* (a string representing an existing file)
    A file containing the behavioural performances of the subject

    ::

      ex. /home/username/data/derived_dat/correct_response.csv

- *patient_info*
    A dictionary for entering patient data.
       - PatientName: the patient's code name
       - Pathology: the patient's pathology
       - Age: patient's age
       - Sex: patient's sex
       - LateralizationPathology: patient's pathology lateralization
       - DominantHand: patient's dominant hand

    ::

      ex. {'PatientName': '002', 'Pathology': 'epilepsy', 'Age': 57,
           'Sex': 'F', 'LateralizationPathology': 'L', 'DominantHand': 'R'}

**Outputs parameters:**

- *report* (a strings representing a file)
    The generated report (pdf).

    ::

      ex. /home/username/data/derived_data/sub-1_GE2REC_Report_2024_01_24_09_34_58_08.pdf
