# -*- coding: utf-8 -*-
"""
The pipelines used for the CerebroVascular Reactivity evaluation.

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from .ge2rec import Ge2rec  # noqa: F401
from .ge2rec_without_eprime import Ge2rec_without_eprime  # noqa: F401
