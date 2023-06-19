:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

==============================
TemplateFromTemplateFlow brick
==============================

Get template image from TemplateFlow
------------------------------------

TemplateFlow is a repository of neuroimaging templates including spatial mapping across standard space.
Each tempalte is described by a template name (for example: 'MNI152NLin2009cAsym', 'MNIColin27', 'MNIPediatricAsym') and by several optional entities (resolution, suffix, label..).

--------------------------------------

**Optional inputs with default value parameters:**

- *atlas* (a string, optional, default value is '')
    Name of a particular atlas (entity 'atlas' in template path name). Default is ''.
    Example: 'DiFuMo', 'Scahefer2018'...

    ::

      ex. 'DiFuMo'

- *desc* (a string, optional, default value is '')
    Description field (entity 'desc' in template path name).
    Example: 'brain', 'eye', '256dimensions' ...

    ::

      ex. 'brain'

- *in_template* (a string, default value is 'MNI152NLin2009cAsym')
    Template name.

    ::

      ex. 'MNI152NLin2009cAsym'

- *label* (a string, optional, default value is '')
    Label fields.
    Example: 'CSF', 'GM'...

    ::

      ex. 'WM'


- *resolution* (an int, default value is 2)
    Resolution of the template (entity 'res' in template path name).

    ::

      ex. 2

- *suffix* (a string, optional, default value is '')
    Suffix of the template image.
    Example: 'T1w', 'T2w', 'probseg', 'T2map', 'mask'...


    ::

      ex. 'T1W'


**Outputs parameters:**

- *in_template* (a strings representing a file)
    Path of the template.

    ::

      ex. '/home/username/.cache/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz'

-------------

Usefull links:

`TemplateFlow <https://www.templateflow.org/>`_
