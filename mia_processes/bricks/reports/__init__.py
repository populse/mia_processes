# -*- coding: utf-8 -*-
"""
The bricks and functions to calculate the values needed in the reports
and to create the report itself.

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from .processes import (  # noqa: F401
    AnatIQMs,
    BoldIQMs,
    BoldIQMsPlot,
    CarpetParcellation,
    ComputeDVARS,
    FramewiseDisplacement,
    Mean_stdDev_calc,
    Result_collector,
    Spikes,
)
from .reporting import (  # noqa: F401
    ReportAnatMriqc,
    ReportFuncMriqc,
    ReportGroupMriqc,
)
