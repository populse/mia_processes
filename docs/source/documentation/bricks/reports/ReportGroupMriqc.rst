:orphan:

.. toctree::

+--------------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+-------------------------------------------+----------------------------------------------------+

======================
ReportGroupMriqc brick
======================

Generate a group report for mriqc pipeline
------------------------------------------

- Generate a pdf report for all the subjects for which mriqc run for the chosen modality.
- Boxplots are generated for each IQMs.
    - The box extends from the first quartile (Q1) to the third quartile (Q3) of the data, with a line at the median.
    - The whiskers extend from the box by 1.5x the inter-quartile range (IQR).
- A '.tsv' file containing all the IQMs is also generated.

-----------------------------------

**Mandatory inputs parameters:**

- *modality* (anat or func):
    Choose modality between 'anat' and 'func'.

    ::

        ex. 'anat'

**Outputs parameters:**

- *report_pdf* (a string that representing a file):
    The generated report (a pdf file).

    ::

        ex. '/home/username/data/derived_data/.pdf'

- *out_tsv* (a string that representing a file):
    A '.tsv' file  containing all the IQMs.

    ::

        ex. '/home/username/data/derived_data/.tsv'


-------------

Useful links:

`matplotlib boxplot <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.boxplot.html>`_
