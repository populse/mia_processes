# -*- coding: utf-8 -*-

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
    in the populse.mia_processes website
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
            "resample1",
            "mia_processes.bricks.preprocess.others.processing.Resample1",
        )
        self.add_process(
            "resample2",
            "mia_processes.bricks.preprocess.others.processing.Resample2",
        )
        self.add_process(
            "mean_stddev_calc",
            "mia_processes.bricks.reports.processes.Mean_stdDev_calc",
        )
        self.nodes["mean_stddev_calc"].set_plug_value(
            "parametric_maps", traits.Undefined
        )
        self.nodes[
            "mean_stddev_calc"
        ].process.parametric_maps = traits.Undefined
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
        self.nodes[
            "result_collector"
        ].process.parameter_files = traits.Undefined
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
            "convroi", "convolve_with", "mask_002", is_optional=False
        )
        self.add_link("mask_002->resample1.files_to_resample")
        self.add_link(
            "concat_to_list_of_list.listOflist->import_data.rois_list"
        )
        self.add_link("import_data.rois_files->convroi.images_to_convolve")
        self.add_link("find_in_list_1.out_file->files_to_list.file1")
        self.add_link("find_in_list_1.out_file->resample1.reference_image")
        self.add_link("find_in_list_2.out_file->files_to_list.file2")
        self.add_link(
            "files_to_list.file_list->mean_stddev_calc.parametric_maps"
        )
        self.add_link("convroi.out_images->mean_stddev_calc.rois_files")
        self.add_link("convroi.out_images->resample2.files_to_resample")
        self.export_parameter(
            "convroi", "out_images", "conv_roi_masks", is_optional=False
        )
        self.add_link("resample1.out_files->resample2.reference_image")
        self.export_parameter(
            "resample2", "out_images", "resample2_masks", is_optional=False
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
                "mask_002",
                "patient_info",
                "resample2_masks",
                "xls_files",
                "conv_roi_masks",
            )
        )

        # nodes positions
        self.node_position = {
            "concat_to_list_of_list": (614.2930418166122, 827.1208223644738),
            "import_data": (925.9736882862239, 848.8631446556915),
            "find_in_list_1": (776.4427376188082, 303.29535996898),
            "find_in_list_2": (856.2337859568112, 567.1330449805765),
            "files_to_list": (1026.2058853774747, 483.98336643917287),
            "convroi": (1089.7332639638682, 685.6142079029513),
            "resample1": (1134.4471373729352, 232.24020431310066),
            "resample2": (1419.5562603305748, 274.20917687569397),
            "inputs": (612.3472403891509, 560.3262591732853),
            "mean_stddev_calc": (1432.1926253095098, 449.1820536125214),
            "concat_to_list": (1748.091183783158, 395.2820325663692),
            "result_collector": (1612.1473004020177, 658.4933845149167),
            "outputs": (1922.8337322988386, 777.3577137671373),
        }

        # nodes dimensions
        self.node_dimension = {
            "concat_to_list_of_list": (143.21875, 110.0),
            "import_data": (147.703125, 180.0),
            "find_in_list_1": (118.84375, 110.0),
            "find_in_list_2": (118.84375, 110.0),
            "files_to_list": (97.640625, 145.0),
            "convroi": (208.0625, 145.0),
            "resample1": (176.84375, 250.0),
            "resample2": (193.921875, 145.0),
            "inputs": (119.17251223929541, 136.0),
            "mean_stddev_calc": (214.171875, 145.0),
            "concat_to_list": (97.96875, 110.0),
            "result_collector": (171.640625, 145.0),
            "outputs": (126.07117003946317, 111.0),
        }

        self.do_autoexport_nodes_parameters = False
