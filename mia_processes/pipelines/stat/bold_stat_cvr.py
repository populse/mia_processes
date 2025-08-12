"""The module for the Bold_stat pipeline.

:Contains:
    :Class:
        - Bold_stat_cvr

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


class Bold_stat_cvr(Pipeline):
    """
    *GLM-based statistical calculation pipeline used in the CVR evaluation*

    Please, see the complete documentation for the `Bold_stat_cvr pipeline
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/stat/Bold_stat_cvr.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "estimatecontrast",
            "mia_processes.bricks.stat.spm.model.EstimateContrast",
        )
        self.add_process(
            "level1design", "mia_processes.bricks.stat.spm.model.Level1Design"
        )
        self.nodes["level1design"].set_plug_value(
            "sess_multi_reg", traits.Undefined
        )
        self.nodes["level1design"].process.sess_multi_reg = traits.Undefined
        self.add_process(
            "estimatemodel",
            "mia_processes.bricks.stat.spm.model.EstimateModel",
        )
        self.add_process(
            "make_a_list", "mia_processes.bricks.tools.Make_A_List"
        )
        self.nodes["make_a_list"].set_plug_value("obj1", traits.Undefined)
        self.nodes["make_a_list"].process.obj1 = traits.Undefined

        # links
        self.export_parameter(
            "make_a_list", "obj1", "regressors", is_optional=False
        )
        self.export_parameter(
            "level1design", "sess_scans", "smoothed_func", is_optional=False
        )
        self.export_parameter(
            "level1design", "mask_image", "mask_003", is_optional=True
        )
        self.export_parameter(
            "estimatecontrast", "spmT_images", is_optional=True
        )
        self.export_parameter(
            "estimatecontrast", "out_spm_mat_file", is_optional=False
        )
        self.add_link("level1design.spm_mat_file->estimatemodel.spm_mat_file")
        self.add_link(
            "estimatemodel.out_spm_mat_file->estimatecontrast.spm_mat_file"
        )
        self.export_parameter("estimatemodel", "beta_images", is_optional=True)
        self.add_link(
            "estimatemodel.beta_images->estimatecontrast.beta_images"
        )
        self.add_link(
            "estimatemodel.residual_image->estimatecontrast.residual_image"
        )
        self.add_link("make_a_list.obj_list->level1design.sess_multi_reg")

        # parameters order

        self.reorder_traits(
            (
                "regressors",
                "smoothed_func",
                "mask_003",
                "spmT_images",
                "out_spm_mat_file",
                "beta_images",
            )
        )

        # nodes positions
        self.node_position = {
            "estimatecontrast": (854.3356190476188, 53.00380952380948),
            "level1design": (171.42247619047617, -102.36076190476183),
            "estimatemodel": (515.0662857142854, 175.10399999999987),
            "make_a_list": (39.18323809523808, 198.2312380952381),
            "inputs": (-158.2314466389599, 157.20000000000002),
            "outputs": (1150.50640218232, 168.01409523809508),
        }

        # nodes dimensions
        self.node_dimension = {
            "estimatecontrast": (247.8125, 425.0),
            "level1design": (290.015625, 880.0),
            "estimatemodel": (246.296875, 390.0),
            "make_a_list": (96.515625, 110.0),
            "inputs": (129.36001223929543, 111.0),
            "outputs": (124.38367003946317, 111.0),
        }

        self.do_autoexport_nodes_parameters = False
