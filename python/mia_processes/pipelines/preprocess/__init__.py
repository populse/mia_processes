# -*- coding: utf-8 -*- #

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from .spatial_preprocessing import Spatial_preprocessing
from .anat_airmask_pipeline import Anat_airmask_pipeline
from .anat_headmask_pipeline import Anat_headmask_pipeline
from .anat_mni_tpms_pipeline import Anat_mni_tpms_pipeline
from .anat_mriqc_pipeline import Anat_mriqc_pipeline
from .anat_skullstrip_pipeline import Anat_skullstrip_pipeline
from .anat_spatial_norm import Anat_spatial_norm
from .bold_hmc_pipeline import Bold_hmc_pipeline
from .bold_mni_align import Bold_mni_align
from .bold_mriqc_pipeline import Bold_mriqc_pipeline
