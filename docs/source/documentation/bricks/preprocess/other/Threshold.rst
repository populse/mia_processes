:orphan:

.. toctree::

+-----------------------------+-------------------------------------------+----------------------------------------------------+
|`Home <../../../index.html>`_|`Documentation <../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------------+-------------------------------------------+----------------------------------------------------+

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

**Inputs parameters:**

- *in_files*
    Path of the scans to be processed. A list of items which are uncompressed files (valid extensions: [.img, .nii, .hdr]).

    ::

      ex. ['/home/ArthurBlair/data/raw_data/c1Anat.nii',
           '/home/ArthurBlair/data/raw_data/c2Anat.nii',
           '/home/ArthurBlair/data/raw_data/c3Anat.nii',
           '/home/ArthurBlair/data/raw_data/c4Anat.nii',
           '/home/ArthurBlair/data/raw_data/c5Anat.nii']

- *threshold*
    Value of the rate (a float between 0 and 1) defining the applied threshold (everything above becomes 1 and everything below becomes
    0, as usual for a mask). The value of the actual threshold is thus the maximum value of the image multiplied by this rate.
    ::

      ex. 0.3

- *GM_filtering*
    A boolean to indicate if the input images list (in_files input parameter) is filtered to retain only grey matter (based on the name of the
    image which must contain c1 in the first part, according to the SPM syntax).
    
    ::

      ex. True

- *suffix*
    Suffix of the output image (a string).

    ::

       ex. _002

- *prefix*
    Prefix of the output image (a string).

    ::

      ex. ""

**Outputs parameters:**

- *out_files*
    Scan name after application of the threshold. The input files extension (*in_files* parameter) is kept.

    ::

      ex. /home/ArthurBlair/data/derived_data/c1Anat_002.nii
