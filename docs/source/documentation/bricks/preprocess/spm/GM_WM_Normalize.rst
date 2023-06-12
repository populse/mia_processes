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

=====================
GM_WM_Normalize brick
=====================

Normalises only the grey or/and white matter from a set of images
-----------------------------------------------------------------

The GM_WM_Normalize brick is equivalent to the `Normalize12 <Normalize12.html>`_ brick from the mia_processes library, with an additional filter (`in_filter` parameter) for grey and/or white matter, applied to the `apply_to_files` input parameter [#label]_. A list of native space probability maps (`apply_to_files` parameter) and an inverse deformation field (`deformation_file` parameter) may have been, for example, obtained after a previous segmentation. The `jobtype` input parameter is therefore fixed to ``write``.

The `in_filter` parameter can take the following values:
  - ``GM`` (grey matter images will be selected before to be normalised).
  - ``WM`` (white matter images will be selected before to be normalised).
  - ``GM & WM`` (grey and white matter images will be selected before to be normalised).
  - ``GM + WM`` (grey and white matter images will be selected and then added together before to be normalised).

For the `deformation_file`, `apply_to_files`, `write_bounding_box`, `write_voxel_sizes`, `write_interp` input parameters and `normalized_files` output parameter, please see the `Normalize12 <Normalize12.html>`_ brick of the mia_processes library.

-------------

.. [#label]  The selection of the images is made using the syntax of SPM (*c1\** for grey matter and *c2\** for white matter).
