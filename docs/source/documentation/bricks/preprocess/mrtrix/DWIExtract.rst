:orphan:

.. toctree::

+--------------------------------+----------------------------------------------+----------------------------------------------------+
|`Home <../../../../index.html>`_|`Documentation <../../../documentation.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+--------------------------------+----------------------------------------------+----------------------------------------------------+

================
DWIExtract brick
================

Extract shell or b=0 volumes from DWI data
------------------------------------------

Extract diffusion-weighted volumes, b=0 volumes or certain shells from a DWI dataset

(mrtrix dwiextract command)

--------------------------------------

**Mandatory inputs parameters:**

- *in_file* (a string representing an existing file)
    Input DWI image (valid extensions: [.mif]).

    ::

      ex. '/home/username/data/raw_data/DWI.mif'

**Optional inputs with default value parameters:**

- *bzero* (a boolean, default value is True, optional)
    Extract all b=0 volumes

    ::

      ex. True

- *nobzero* (a boolean, default value is False, optional)
    Extract all non b=0 volumes

    ::

      ex. False

- *singleshell* (a boolean, default value is False, optional)
    Force a single-shell (single non b=0 shell) output. This will include b=0 volumes, if present.

    Use with bzero option to enforce presence of b=0 volumes (error if not present) or with no_bzero option to exclude them.

    ::

      ex. False

**Optional inputs parameters:**

- *shell* (a list of items which are a float, optional)
    Specify one or more gradient shells to use during processing.
    Default is Undefined (ie parameter not used)

    ::

      ex. [3000.0]

**Outputs parameters:**

- *out_file* (a pathlike object or string representing a file)
    The output bias corrected DWI image

    ::

      ex. '/home/username/data/derived_data/DWI_bzero.mif'


-------------

Usefull links:

`mrtrix dwiextract <https://mrtrix.readthedocs.io/en/latest/reference/commands/dwiextract.html>`_

`mrtrix dwiextract - nipype <https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.mrtrix3.utils.html#dwiextract>`_
