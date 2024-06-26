# -*- coding: utf-8 -*-

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

# other import
import traits.api as traits

# Capsul import
from capsul.api import Pipeline


class Spatial_mask(Pipeline):
    """
    *Makes the grey matter native mask for cerebrovascular reserve analyse*

    Please, see the complete documentation for the `Spatial_mask brick in
    the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Spatial_mask.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "gm_wm_normalize",
            "mia_processes.bricks.preprocess.spm.GM_WM_Normalize",
        )

        self.add_process(
            "threshold_1", "mia_processes.bricks.preprocess.others.Threshold"
        )
        self.nodes["threshold_1"].process.in_files = traits.Undefined
        self.nodes["threshold_1"].process.prefix = traits.Undefined

        self.add_process(
            "smooth", "mia_processes.bricks.preprocess.spm.Smooth"
        )
        self.nodes["smooth"].process.in_files = traits.Undefined
        self.nodes["smooth"].process.data_type = 4

        self.add_process(
            "threshold_2", "mia_processes.bricks.preprocess.others.Threshold"
        )
        self.nodes["threshold_2"].process.in_files = traits.Undefined
        self.nodes["threshold_2"].process.suffix = traits.Undefined
        self.nodes["threshold_2"].process.prefix = "mask_"

        self.add_process(
            "resample1", "mia_processes.bricks.preprocess.others.Resample1"
        )
        self.nodes["resample1"].process.files_to_resample = traits.Undefined
        self.nodes["resample1"].process.prefix = traits.Undefined

        # links
        self.export_parameter(
            "gm_wm_normalize", "apply_to_files", "native_class_images"
        )
        self.export_parameter("gm_wm_normalize", "deformation_file")
        self.add_link("gm_wm_normalize.normalized_files->threshold_1.in_files")
        self.add_link("threshold_1.out_files->smooth.in_files")
        self.add_link("smooth.smoothed_files->threshold_2.in_files")
        self.add_link("threshold_2.out_files->resample1.files_to_resample")
        self.export_parameter("resample1", "out_files", "mask_003")
        self.export_parameter("resample1", "reference_image", "smoothed_func")

        # nodes positions
        self.node_position = {
            "resample1": (1067.6501882850332, 222.64611569347358),
            "gm_wm_normalize": (136.9527473844211, 309.8135411411285),
            "threshold_2": (866.6321673811615, 369.0968621492992),
            "smooth": (607.8675032581176, 370.0759775107607),
            "outputs": (1284.411220634047, 244.59809770510992),
            "inputs": (-49.131814528226414, 168.40674594714966),
            "threshold_1": (418.5605170436202, 453.4652656971715),
        }

        # nodes dimensions
        self.node_dimension = {
            "inputs": (157.4068872392954, 111.0),
            "threshold_1": (152.609375, 215.0),
            "smooth": (218.046875, 215.0),
            "threshold_2": (152.609375, 215.0),
            "resample1": (176.84375, 250.0),
            "outputs": (81.11804503946317, 86.0),
            "gm_wm_normalize": (244.0625, 250.0),
        }

        self.do_autoexport_nodes_parameters = False
