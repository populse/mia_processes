:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Template brick
============

Get template image from templateflow

--------------------------------------

>>> from mia_processes.bricks.preprocess.others import Template

>>> Template.help()

**Inputs parameters:**

- *in_template* (a string)
    Template name. Default is 'MNI152NLin2009cAsym'.

    ::

      ex. 'MNI152NLin2009cAsym'

- *resolution* (an int)
    Resolution of the template. Default is 2.
    
    ::

      ex. 2

- *suffix* (a string, optional=True)
   Suffix of output image. Default is ''.
    
    ::

      ex. ''

- *atlas* (a string, optional=True)
    Name of a particular atlas. Default is ''.

    ::

      ex. ''

- *desc* (a string, optional=True)
    Description field. Default is ''.

    ::

      ex. ''

- *label* (a string, optional=True)
    Label fields. Default is ''.

    ::

      ex. ''



**Outputs parameters:**

- *in_template* (a strings representing a file)
    Path of the template 
    
    ::

      ex. ''

-------------

Usefull links:
`SanitizeImage niworflow <https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/header.py#L394>`_
