:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

===============
GE2REC pipeline
===============

Interactive mapping of language and memory with the GE2REC protocol
--------------------------------------------------------------------

The purpose of this pipeline is to process functional MRI data from the GE2REC protocol.

This pipeline require the E-Prime (E-Prime 3.0 Psychology Software Tools) file of the experiment
(behavioural data are used).

The GE2REC protocol is used to map the interaction between the functions of language and
memory functions (LMN language memory network) at an individual level.

This protocol consists of three interdependent tasks:

- The first task is a block sequence of sentence generation with implicit encoding (GE). Subjects heard words through headphones and have to generate sentences in internal language. During the control periods, a pseudoword is broadcast to the subjects  and they have to listen the pseudoword and not to talk covertly. The run also include rest period (with a fixation cross displayed).

- The second task is a recognition event paradigm (RECO). Images are presented to the subjects in a pseudo-random modeand they have to indicate whether they recognise the images of the objects whose names had been broadcast during the first GE task. The choice was binary. Either the subject indicates that they recognised the image (OLD), or that it was a new item (NEW). The run also included 40 control images showing the button that needed to be pressed in order to control for the motor activations during button pressing.

- The third task is a block paradigm (RA) in auditory modality. Participants hear the words from the GE task and must remember and repeat the sentences previously generated in the first task.

In this pipeline, the fMRI data are analyzed using a GLM approach after pre-processing.
GE and RA runs are analysed as a block design while RECO run is analysed as an event-related design.

Lateralization indexes are computed using an iterative approach on the GE run.

Based on the responses during the RECO run, we calculated behavioral performances during this task.
The encoding performance during GE run was indirectly determined via those responses and
the GE run is analysed using the encoding performance as a regressor.

A report is generated at the end of the analysis with the main statistical results obtained.

**References**

- Banjac S., Interactive cartography of language and memory in patients with focal and pharmaco- resistant pilepsy. Multimodal assessment. Psychology. Université Grenoble Alpes [2020-..], 2021. English. NNT : 2021GRALS037 . tel- 03641888.

- Banjac S and al. Mapping of Language-and-Memory Networks in Patients With Temporal Lobe Epilepsy by Using the GE2REC Protocol. Front Hum Neurosci. 2022 Jan 6;15:752138. doi: 10.3389/fnhum.2021.752138. PMID: 35069148; PMCID:PMC8772037.

- Banjac and al. (2020). Interactive mapping of language and memory with the GE2REC protocol. Brain Imaging Behav. 15, 1562‒1579. doi: 10.1007/s11682-020-00355-x.


**Pipeline insight**

| The GE2REC pipeline combines the following pipelines and processes:
|   - `Preprocessing <../preprocess/Bold_spatial_preprocessing3.html>`_
|   - `Get E-Prime information <../../bricks/tools/Get_Eprime_info_GE2REC.html>`_
|   - `Level1Design <../../bricks/stats/spm/Level1Design.html>`_
|   - `EstimateModel <../../bricks/stats/spm/EstimateModel.html>`_
|   - `EstimateContrast <../../bricks/stats/spm/EstimateContrast.html>`_
|   - `Compute lateralization index <../../bricks/reports/LateralizationIndexCurve.html>`_
|   - `Report <../../bricks/reports/ReportGE2REC.html>`_

.. image:: ../../images/GE2REC_pipeline.png
    :width: 800
    :alt: GE2REC pipeline

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

- *eprime_file* (a string representing an existing file)
    Eprime file (valid extensions: [.xlsx]).

    ::

      ex. '/home/username/data/raw_data/eprime.xlsx'


**Outputs parameters:**

- *report*
    The output generated report (.pdf).

    ::

      ex. '/home/username/data/derived_data/PatientName_GE2REC_Report_2024_03_19_11_01_15_04.pdf'
