:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

========================================
GE2REC_without_behavioural_data pipeline
========================================

Interactive mapping of language and memory with the GE2REC protocol without behavioural data
--------------------------------------------------------------------------------------------

The purpose of this pipeline is to process functional MRI data from the GE2REC protocol.
This pipeline do not require the E-Prime (E-Prime 3.0 Psychology Software Tools) file of the experiment
(Behavourial data are not analysed).

See `GE2REC pileline <./GE2REC.html>`_ for explanation about GE2REC protocol.

In this pipeline, the fMRI data are analyzed using a GLM approach after pre-processing.
GE and RA runs are analysed as a block design while RECO run is analysed as an event-related design.

Lateralization indexes are computed using an iterative approach on the GE run.

A report is generated at the end of the analysis with the main statistical results obtained.

**Pipeline insight**

| The GE2REC pipeline combines the following pipelines and processes:
|   - `Preprocessing <../preprocess/Bold_spatial_preprocessing3.html>`_
|   - `Level1Design <../../bricks/stats/spm/Level1Design.html>`_
|   - `EstimateModel <../../bricks/stats/spm/EstimateModel.html>`_
|   - `EstimateContrast <../../bricks/stats/spm/EstimateContrast.html>`_
|   - `Compute lateralization index <../../bricks/reports/LateralizationIndexCurve.html>`_
|   - `Report <../../bricks/reports/ReportGE2REC.html>`_


**Inputs parameters**

- *anat_file* (a string representing an existing file)
    An anatomical image
    An existing, uncompressed file (valid extensions: [.nii]).

    ::

      ex. '/home/username/data/raw_data/Anat.nii'

- *func_gene_file* (a string representing an existing file)
    Functional image of the generation run (valid extensions: [.nii]).

    ::

      ex. '/home/username/data/raw_data/Func_gene.nii'

- *func_reco_file* (a string representing an existing file)
    Functional image of the recognition run (valid extensions: [.nii]).

    ::

      ex. '/home/username/data/raw_data/Func_reco.nii'

- *func_recall_file* (a string representing an existing file)
    Functional image of the recall run (valid extensions: [.nii]).

    ::

      ex. '/home/username/data/raw_data/Func_recall.nii'

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

- *report*
    The output generated report (.pdf).

    ::

      ex. '/home/username/data/derived_data/PatientName_GE2REC_Report_2024_03_19_11_01_15_04.pdf'
