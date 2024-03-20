:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

=========================
Bold_signals_roi pipeline
=========================

This pipeline allows to get plots of the BOLD average signal in severals ROI (defined by AssemblyNet).


*This pipeline requires Docker to run AssemblyNetDocker brick*


--------------------------------------

**Pipeline insight**

| Bold_signals_roi pipeline combines the following pipelines and processes:
|   - `AssemblyNetDocker <../../bricks/preprocess/volbrain/AssemblyNetDocker.html>`_
|   - `EpiReg <../../bricks/preprocess/fsl/EpiReg.html>`_
|   - `Automask <../../bricks/preprocess/afni/Automask.html>`_
|   - `BetSurfacesExtraction <../../bricks/preprocess/fsl/BetSurfacesExtraction.html>`_
|   - `ConvertXFM  <../../bricks/preprocess/fsl/ConvertXFM.html>`_
|   - `Flirt <../../bricks/preprocess/fsl/Flirt.html>`_ (default values : apply_xfm = True, interpolation = 'nearestneighbour')
|   - `LabelsCorrespondence <../../bricks/preprocess/volbrain/LabelsCorrespondence.html>`_
|   - `GetLabels <../../bricks/preprocess/volbrain/GetLabels.html>`_
|   - `ExtractSignalROI <../../bricks/preprocess/others/ExtractSignalROI.html>`_
|   - `ExtractROIbyLabel <../../bricks/preprocess/others/ExtractROIbyLabel.html>`_
|   - `PlotSignalROI <../../bricks/reports/PlotSignalROI.html>`_

**Mandatory inputs parameters**

- *bold* (a string representing an existing file)
    A functional image (BOLD). An existing file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/func.nii'


- *anat* (a string representing an existing file)
    A T1w image. An existing file (valid extensions: [.nii, .nii.gz]).

    ::

      ex. '/home/username/data/raw_data/T1w.nii'


**Optional inputs parameters with default values**

- *labels_structures* (a list of integer, default value is [47, 48])
    List of structures labels for which the average signals will be extracted.

    The defaul value is [47, 48], this correponds to le left and the right hippocampus.

    ::

      ex. [47, 48]

**Outputs parameters:**

- *out_png_tissues* (a string representing a file)
    Out png file with a plot of the average signal for each tissues.

    ::

      ex. '/home/username/data/derived_data/func_extracted_signals_tissues.png'


- *out_png_lobes* (a string representing a file)
    Out png file with a plot of the average signal for each lobes.

    ::

      ex. '/home/username/data/derived_data/func_extracted_signals_lobes.png'


- *out_png_macrostructures* (a string representing a file)
    Out png file with a plot of the average signal for each macrostructures.

    ::

      ex. '/home/username/data/derived_data/func_extracted_signals_macrostructures.png'


- *out_png_structures* (a string representing a file)
    Out png file with a plot of the average signal for each structures given in the labels_structures parameter.

    ::

      ex. '/home/username/data/derived_data/func_extracted_signals_47_48.png'

-------------

Usefull links:

`volBrain Assemblynet <https://github.com/volBrain/AssemblyNet>`_

`Docker <https://docs.docker.com/get-docker/>`_
