:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

======================
Mean_stdDev_calc brick
======================

Makes the mean and standard deviation of parametric maps.
---------------------------------------------------------

    - The `parametric_maps` are first convolved with the ROIs corresponding
      to `doublet_list`.
    - ROIs are defined from `doublet_list` parameter as
      doublet_list[0][0] + doublet_list[0][1] + '.nii',
      doublet_list[1][0] + doublet_list[1][1] + '.nii',
      etc.
    - To work correctly, the "/roi\_/`PatientName`/convROI_BOLD2" directory
      must exist and contain the results of a previous convolution (using for
      example the `preprocess.others.ConvROI` brick).
    - To work correctly, the database entry for the `parametric_maps` items must
      have the `PatientName` tag filled in.


--------------------------------------

**Mandatory inputs parameters:**

- *doublet_list* (a list of lists)
    A list of lists containing doublets of strings.

    ::

      ex. [["ROI_OCC", "_L"], ["ROI_OCC", "_R"], ["ROI_PAR", "_l"]]


- *parametric_maps* (a list of existing files)
   A list of files .

    ::

      ex. []


**Outputs parameters:**

- * mean_out_files* (a list of files)
    A list of .txt files with the calculated average for each ROI determined after convolution.

    ::

      ex. []

- * std_out_files* (a list of files)
    A list of .txt files with the standard deviation for each ROI determined after convolution.

    ::

      ex. []
