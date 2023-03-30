:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
TemplateFromTemplateFlow brick
============

TemplateFlow is a repository of neuroimaging templates including spatial mapping across standard space. 
Each tempalte is described by a template name (for example: 'MNI152NLin2009cAsym', 'MNIColin27', 'MNIPediatricAsym') and by several optional entities (resolustion, suffix, label..). 

This bricks allows to get template image from TemplateFlow.  


--------------------------------------

>>> from mia_processes.bricks.preprocess.others import TemplateFromTemplateFlow

>>> TemplateFromTemplateFlow.help()

**Optional inputs with default value parameters:**

- *atlas* (a string, optional)
    Name of a particular atlas (entity 'atlas' in template path name). Default is ''.
    Example: 'DiFuMo', 'Scahefer2018'...

    ::

      default value. ''
  
  - *desc* (a string, optional)
    Description field (entity 'desc' in template path name). 
    Example: 'brain', 'eye', '256dimensions' ...

    ::

      default value. ''

  - *in_template* (a string)
    Template name. 

    ::

      default value. 'MNI152NLin2009cAsym'

- *label* (a string, optional)
    Label fields.
    Example: 'CSF', 'GM'...

    ::

      default value. ''


- *resolution* (an int)
    Resolution of the template (entity 'res' in template path name). 
    
    ::

      default value. 2

- *suffix* (a string, optional)
    Suffix of the template image. 
    Example: 'T1w', 'T2w', 'probseg', 'T2map', 'mask'...

    
    ::

      default value. ''


**Outputs parameters:**

- *in_template* (a strings representing a file)
    Path of the template. 
    
    ::

      ex. '/home/username/.cache/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz'

-------------

Usefull links:
`TemplateFlow <https://www.templateflow.org/>`_
