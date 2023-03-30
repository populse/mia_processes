:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

============
Enhance brick
============

Image enhancing. The aims is to resample signal excess voxels. 

For each images, voxels with signal excess (superior to 99.8%) are set to a random value among the voxel's values superior to the median.

Be carrefull, if the suffix and prefix input parameters are not defined or consist only of one or more white spaces, the input parameter will be overwritten.

Adapted from https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/anatomical.py#L974
--------------------------------------

>>> from mia_processes.bricks.preprocess.others import Enhance
>>> Enhance.help()

**Mandatory inputs parameters:**

- *in_files* (a list of pathlike objects or strings representing an existing file)
    Input images (valid extensions: [.nii, .nii.gz]).
    
    ::

      ex. ['/home/username/data/raw_data/T1w.nii']

**Optional inputs with default value parameters:**

- *prefix* (a string, optional)
    Prefix of the output image. Default is ''.
    
    ::

      default value. ''


- *suffix* (a string, optional)
   Suffix of output image. Default is '_enh'.
    
    ::

      default value. '_enh'


**Outputs parameters:**

- *out_files* (a list of pathlike objects or strings representing a file)
    Output Images (extensions: [.nii, .nii.gz]).
    
    ::

      ex. ['/home/username/data/derived_data/T1w_enh.nii']

-------------


Usefull links:
`Enhance image mriqc <https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/anatomical.py#L974>`_