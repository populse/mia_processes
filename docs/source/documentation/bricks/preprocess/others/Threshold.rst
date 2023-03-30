:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

.. line break
.. |br| raw:: html

   <br />

.. thin space
.. |ws1| raw:: html

   &thinsp;

.. em space

.. |ws2| raw:: html

   &emsp;

.. en space

.. |ws3| raw:: html

   &ensp;

.. non-breakable space

.. |ws4| raw:: html

   &nbsp;

===============
Threshold brick
===============

Makes a binary mask image at a given threshold
----------------------------------------------

>>> from mia_processes.bricks.preprocess.other import Threshold
>>> Threshold.help()

**Mandatory inputs parameters:**

- *in_files*
    Path of the scans to be processed. A list of items which are uncompressed files (valid extensions: [.img, .nii, .hdr]).

    ::

      ex. ['/home/username/data/raw_data/c1Anat.nii',
           '/home/username/data/raw_data/c2Anat.nii',
           '/home/username/data/raw_data/c3Anat.nii',
           '/home/username/data/raw_data/c4Anat.nii',
           '/home/username/data/raw_data/c5Anat.nii']

**Optional inputs with default value parameters:**

- *GM_filtering* (a boolean, optional)
    A boolean to indicate if the input images list (in_files input parameter) is filtered to retain only grey matter (based on the name of the
    image which must contain c1 in the first part, according to the SPM syntax).
    
    ::

      ex. True

- *prefix* (a string, optional)
    Prefix of the output image.

    ::

      ex. ""

- *suffix* (a string, optional)
    Suffix of the output image.

    ::

       ex. _002

- *threshold* (a float between 0 and 1, optional)
    Value of the rate  defining the applied threshold (everything above becomes 1 and everything below becomes
    0, as usual for a mask). The value of the actual threshold is thus the maximum value of the image multiplied by this rate.
    ::

      defalut value. 0.3

**Outputs parameters:**

- *out_files*
    Scan name after application of the threshold. The input files extension (*in_files* parameter) is kept.

    ::

      ex. /home/username/data/derived_data/c1Anat_002.nii
