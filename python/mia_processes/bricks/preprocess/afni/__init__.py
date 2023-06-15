# -*- coding: utf-8 -*-
"""
The atomic calculations from afni.

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from .processes import (  # noqa: F401
    GCOR,
    Automask,
    Calc,
    CalcDropTRs,
    Despike,
    FWHMx,
    OutlierCount,
    QualityIndex,
    RefitDeoblique,
    SkullStrip,
    TShift,
    TStatMean,
    Volreg,
)
