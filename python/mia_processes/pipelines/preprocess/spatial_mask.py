# -*- coding: utf-8 -*- #

"""The module for the Spatial_mask pipeline.

:Contains:
    :Class:
        - Spatial_mask

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# Capsul import
from capsul.api import Pipeline

# other import
import traits.api as traits

# populse_mia import
#from populse_mia.project.filter import Filter


class Spatial_mask(Pipeline):
    """
    *Makes the grey matter native mask for cerebrovascular reserve analyse*
    User_processes_ECdev.studies.spatial_mask_1.Spatial_mask_1()

    Please, see the complete documentation for the Spatial_mask_1 pipeline in the
    https://gricad-gitlab.univ-grenoble-alpes.fr/condamie/ec_dev web site

    """
 
    def pipeline_definition(self):
        # nodes
        self.add_process("gm_wm_normalize1", "mia_processes.bricks.preprocess.spm.GM_WM_Normalize")

        self.add_process("threshold1", "mia_processes.bricks.preprocess.others.Threshold")
        self.nodes["threshold1"].process.in_files = traits.Undefined
        self.nodes["threshold1"].process.prefix = traits.Undefined

        self.add_process("smooth1", "mia_processes.bricks.preprocess.spm.Smooth")
        self.nodes["smooth1"].process.in_files = traits.Undefined

        self.add_process("threshold2", "mia_processes.bricks.preprocess.others.Threshold")
        self.nodes["threshold2"].process.in_files = traits.Undefined
        self.nodes["threshold2"].process.suffix = traits.Undefined
        self.nodes["threshold2"].process.prefix = 'mask_'

        self.add_process("resample1", "mia_processes.bricks.preprocess.others.Resample_1")
        self.nodes["resample1"].process.files_to_resample = traits.Undefined
        self.nodes["resample1"].process.prefix = traits.Undefined

        # links
        self.export_parameter("gm_wm_normalize1", "apply_to_files", "native_class_images")
        self.export_parameter("gm_wm_normalize1", "deformation_file")
        self.add_link("gm_wm_normalize1.normalized_files->threshold1.in_files")
        self.add_link("threshold1.out_files->smooth1.in_files")
        self.add_link("smooth1.smoothed_files->threshold2.in_files")
        self.add_link("threshold2.out_files->resample1.files_to_resample")
        self.export_parameter("threshold2", "out_files", "mask_002")
        self.export_parameter("resample1", "out_files", "mask_003")
        self.export_parameter("resample1", "reference_image", "smoothed_func")

        # nodes positions
        self.node_position = {
            "resample1": (660.8743625, 314.2747999999999),
            "gm_wm_normalize1": (121.68129999999996, -140.0),
            "threshold2": (402.93549225419565, 203.88756771783517),
            "smooth1": (742.5339029207395, -156.09480055247312),
            "outputs": (827.6561125000001, 151.5811),
            "inputs": (-40.80193413672403, 164.24180575139843),
            "threshold1": (440.7735314209599, -168.49913686833992),
        }


        # nodes dimensions
        self.node_dimension = {
            "inputs": (155.453125, 145.0),
            "normalize_spatial_mask1": (244.0625, 215.0),
            "threshold1": (152.609375, 215.0),
            "smooth1": (218.046875, 215.0),
            "threshold2": (152.609375, 215.0),
            "resample1": (176.84375, 250.0),
            "outputs": (81.11804503946317, 110.0),
        }

        self.do_autoexport_nodes_parameters = False

