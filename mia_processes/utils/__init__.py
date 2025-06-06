# -*- coding: utf-8 -*-
"""
Utilities and tools used across mia_processes



Contains:
    Modules:
        - report.py
        - tools.py
"""

###############################################################################
# Populse_mia - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
###############################################################################

from .tools import (  # noqa: F401
    PageNumCanvas,
    ReportLine,
    checkFileExt,
    mriqc_get_all_run,
    mriqc_group_iqms_tsv,
    plot_boxplot_points,
    plot_qi2,
    plot_realignment_parameters,
    plot_segmentation,
    plot_slice_planes,
)

# Prevent circular import (PageNumCanvas)
# isort: off
from .report import Report  # noqa: F401

# isort: on
