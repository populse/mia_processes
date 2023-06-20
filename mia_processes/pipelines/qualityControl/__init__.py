# -*- coding: utf-8 -*-
"""
The pipelines dedicated to the quality control measurement.

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from .anat_mriqc import Anat_mriqc  # noqa: F401
from .bold_mriqc import Bold_mriqc  # noqa: F401
