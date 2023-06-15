# -*- coding: utf-8 -*-
"""
The various pipelines used for the pre-treatment stage.

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from .anat_airmask import Anat_airmask  # noqa: F401
from .anat_headmask import Anat_headmask  # noqa: F401
from .anat_mni_tpms import Anat_mni_tpms  # noqa: F401
from .anat_skullstrip import Anat_skullstrip  # noqa: F401
from .anat_skullstrip_synthstrip import (  # noqa: F401
    Anat_skullstrip_synthstrip,
)
from .anat_spatial_norm import Anat_spatial_norm  # noqa: F401
from .bold_hmc import Bold_hmc  # noqa: F401
from .bold_mni_align import Bold_mni_align  # noqa: F401
from .bold_spatial_preprocessing1 import (  # noqa: F401
    Bold_spatial_preprocessing1,
)
from .bold_spatial_preprocessing2 import (  # noqa: F401
    Bold_spatial_preprocessing2,
)
from .spatial_mask import Spatial_mask  # noqa: F401
