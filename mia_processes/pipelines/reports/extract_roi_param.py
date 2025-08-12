"""blabla.

The purpose of this module is to blabla.

:Contains:
    :Class:
        - Extract_roi_param

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# Other import
import traits.api as traits

# Capsul import
from capsul.api import Pipeline


class Extract_roi_param(Pipeline):
    """
    *Produces gray matter masks for various ROIs and means, standard
    deviations, laterality indexes for beta and spmT values in these ROIs*

    Please, see the complete documentation for the `Extract_roi_param brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/reports/Extract_roi_param.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "concat_to_list_of_list",
            "mia_processes.bricks.tools.tools.Concat_to_list_of_list",
        )
        self.add_process(
            "import_data", "mia_processes.bricks.tools.tools.Import_Data"
        )
        self.add_process(
            "find_in_list_1", "mia_processes.bricks.tools.tools.Find_In_List"
        )
        self.add_process(
            "find_in_list_2", "mia_processes.bricks.tools.tools.Find_In_List"
        )
        self.add_process(
            "files_to_list", "mia_processes.bricks.tools.tools.Files_To_List"
        )
        self.nodes["files_to_list"].set_plug_value("file1", traits.Undefined)
        self.nodes["files_to_list"].set_plug_value("file2", traits.Undefined)
        self.add_process(
            "convroi",
            "mia_processes.bricks.preprocess.others.processing.ConvROI",
        )
        self.add_process(
            "mean_stddev_calc",
            "mia_processes.bricks.reports.processes.Mean_stdDev_calc",
        )
        self.nodes["mean_stddev_calc"].set_plug_value(
            "parametric_maps", traits.Undefined
        )
        self.nodes["mean_stddev_calc"].process.parametric_maps = (
            traits.Undefined
        )
        self.add_process(
            "concat_to_list",
            "mia_processes.bricks.tools.tools.Concat_to_list",
        )
        self.add_process(
            "result_collector",
            "mia_processes.bricks.reports.processes.Result_collector",
        )
        self.nodes["result_collector"].set_plug_value(
            "parameter_files", traits.Undefined
        )
        self.nodes["result_collector"].process.parameter_files = (
            traits.Undefined
        )
        self.nodes["result_collector"].set_plug_value("laterality_index", True)

        # links
        self.export_parameter(
            "find_in_list_1", "in_list", "spmT_images", is_optional=False
        )
        self.add_link("spmT_images->import_data.file_in_db")
        self.export_parameter(
            "find_in_list_2", "in_list", "beta_images", is_optional=False
        )
        self.export_parameter(
            "convroi", "convolve_with", "mask_003", is_optional=False
        )
        self.add_link(
            "concat_to_list_of_list.listOflist->import_data.rois_list"
        )
        self.add_link("import_data.rois_files->convroi.images_to_convolve")
        self.add_link("find_in_list_1.out_file->files_to_list.file1")
        self.add_link("find_in_list_2.out_file->files_to_list.file2")
        self.add_link(
            "files_to_list.file_list->mean_stddev_calc.parametric_maps"
        )
        self.add_link("convroi.out_images->mean_stddev_calc.rois_files")
        self.export_parameter(
            "convroi", "out_images", "conv_roi_masks", is_optional=False
        )
        self.add_link("mean_stddev_calc.mean_out_files->concat_to_list.list1")
        self.add_link("mean_stddev_calc.std_out_files->concat_to_list.list2")
        self.add_link(
            "concat_to_list.out_list->result_collector.parameter_files"
        )
        self.export_parameter(
            "result_collector", "out_files", "xls_files", is_optional=False
        )
        self.export_parameter(
            "result_collector", "patient_info", is_optional=True
        )

        # parameters order

        self.reorder_traits(
            (
                "spmT_images",
                "beta_images",
                "mask_003",
                "patient_info",
                "xls_files",
                "conv_roi_masks",
            )
        )

        # nodes positions
        self.node_position = {
            "concat_to_list_of_list": (600.7683536977197, 728.759454227074),
            "import_data": (838.6779740642814, 754.1903278234439),
            "find_in_list_1": (932.5914095369304, 358.62362954626724),
            "find_in_list_2": (934.9228804667309, 484.75539916550395),
            "files_to_list": (1114.731116701135, 420.0484771498629),
            "convroi": (1065.1429219295185, 724.9587551579114),
            "inputs": (611.1177232874335, 426.308895086078),
            "mean_stddev_calc": (1271.1258849845178, 433.198331290194),
            "concat_to_list": (1548.909413304923, 471.51209287285405),
            "result_collector": (1593.7045438762555, 614.2307688530867),
            "outputs": (1776.5211971944557, 773.6691624619851),
        }

        # nodes dimensions
        self.node_dimension = {
            "concat_to_list_of_list": (143.21875, 110.0),
            "import_data": (147.703125, 180.0),
            "find_in_list_1": (118.84375, 110.0),
            "find_in_list_2": (118.84375, 110.0),
            "files_to_list": (97.640625, 145.0),
            "convroi": (208.0625, 145.0),
            "inputs": (119.17251223929541, 136.0),
            "mean_stddev_calc": (214.171875, 180.0),
            "concat_to_list": (97.96875, 110.0),
            "result_collector": (171.640625, 145.0),
            "outputs": (113.32117003946317, 86.0),
        }

        self.do_autoexport_nodes_parameters = False
