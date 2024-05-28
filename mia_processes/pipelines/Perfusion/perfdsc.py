# -*- coding: utf-8 -*-

"""The pipelines library of the mia_processes package.

The Perfdsc pipeline is designed for cerebral perfusion measurement
using Dynamic Susceptibility Contrast MRI.

:Contains:
    :Class:
        - Perfdsc

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################
from capsul.api import Pipeline


class Perfdsc(Pipeline):
    """
    *Cerebral perfusion measurement by DSC MRI*

    Please, see the complete documentation for the
    `Perfdsc pipeline in the mia_processes website:
    <https://populse.github.io/mia_processes/html/documentation/pipelines/Perfusion/Perfdsc.html>`_
    """

    def pipeline_definition(self):
        """Building the pipeline"""

        # nodes
        self.add_process(
            "1_spatial_preprocessing",
            "mia_processes.pipelines.preprocess.bold_spatial_preprocessing1."
            "Bold_spatial_preprocessing1",
        )
        self.nodes["1_spatial_preprocessing"].process.nodes_activation = {
            "newsegment": True,
            "realign": True,
            "list_duplicate": True,
            "normalize12_1": True,
            "normalize12_2": True,
            "smooth": True,
            "coregister": True,
        }
        # We use exactly the same parameters as in Amigo.
        self.nodes["1_spatial_preprocessing"].process.nodes[
            "realign"
        ].process.register_to_mean = False
        self.nodes["1_spatial_preprocessing"].process.nodes[
            "normalize12_1"
        ].process.write_voxel_sizes = [2.0, 2.0, 2.0]
        self.nodes["1_spatial_preprocessing"].process.nodes[
            "normalize12_1"
        ].process.write_bounding_box = [
            [-78.0, -112.0, -50.0],
            [78.0, 76.0, 85.0],
        ]
        self.nodes["1_spatial_preprocessing"].process.nodes[
            "normalize12_1"
        ].process.write_interp = 1
        self.nodes["1_spatial_preprocessing"].process.nodes[
            "normalize12_2"
        ].process.write_bounding_box = [
            [-78.0, -112.0, -50.0],
            [78.0, 76.0, 85.0],
        ]
        self.nodes["1_spatial_preprocessing"].process.nodes[
            "normalize12_2"
        ].process.write_interp = 1
        self.add_process(
            "2_spatial_mask",
            "mia_processes.pipelines.preprocess.spatial_mask.Spatial_mask",
        )
        self.nodes["2_spatial_mask"].process.nodes_activation = {
            "gm_wm_normalize": True,
            "threshold_1": True,
            "smooth": True,
            "threshold_2": True,
            "resample1": True,
        }

        # links
        self.export_parameter(
            "1_spatial_preprocessing", "anat_file", is_optional=False
        )
        self.export_parameter(
            "1_spatial_preprocessing", "func_files", is_optional=False
        )
        self.add_link(
            "1_spatial_preprocessing.native_class_images->"
            "2_spatial_mask.native_class_images"
        )
        self.add_link(
            "1_spatial_preprocessing.forward_deformation_field->"
            "2_spatial_mask.deformation_file"
        )
        self.add_link(
            "1_spatial_preprocessing.smoothed_func->"
            "2_spatial_mask.smoothed_func"
        )
        self.export_parameter(
            "1_spatial_preprocessing",
            "coregistered_source",
            is_optional=False,
        )
        self.export_parameter("2_spatial_mask", "mask_003", is_optional=False)

        # parameters order

        self.reorder_traits(
            ("mask_003", "anat_file", "func_files", "coregistered_source")
        )

        # nodes positions
        self.node_position = {
            "1_spatial_preprocessing": (-369.0, -166.0),
            "2_spatial_mask": (42.0, -76.0),
            "outputs": (337.00279291855725, -76.0),
            "inputs": (-548.5783216389599, -167.0),
        }

        # nodes dimensions
        self.node_dimension = {
            "1_spatial_preprocessing": (239.515625, 355.0),
            "2_spatial_mask": (200.15625, 145.0),
            "outputs": (135.8125, 110.0),
            "inputs": (92.79751223929541, 110.0),
        }

        self.do_autoexport_nodes_parameters = False
