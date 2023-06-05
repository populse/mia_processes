# -*- coding: utf-8 -*-

"""blabla.

The purpose of this module is to blabla.

:Contains:
    :Class:
        - Rois_1

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


class Rois_1(Pipeline):
    """
    * Blabla.

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "concat_to_list_of_list_1",
            "mia_processes.bricks.tools.tools.Concat_to_list_of_list",
        )
        self.add_process(
            "import_data_1", "mia_processes.bricks.tools.tools.Import_Data"
        )
        self.add_process(
            "find_in_list_1", "mia_processes.bricks.tools.tools.Find_In_List"
        )
        self.add_process(
            "find_in_list_2", "mia_processes.bricks.tools.tools.Find_In_List"
        )
        self.add_process(
            "files_to_list_1", "mia_processes.bricks.tools.tools.Files_To_List"
        )
        self.nodes["files_to_list_1"].set_plug_value("file1", traits.Undefined)
        self.nodes["files_to_list_1"].set_plug_value("file2", traits.Undefined)
        self.add_process(
            "convroi_1",
            "mia_processes.bricks.preprocess.others.processing.ConvROI",
        )
        self.add_process(
            "resample_1_1",
            "mia_processes.bricks.preprocess.others.processing.Resample_1",
        )
        self.add_process(
            "resample_2_1",
            "mia_processes.bricks.preprocess.others.processing.Resample_2",
        )
        self.add_process(
            "mean_stddev_calc_1",
            "mia_processes.bricks.reports.processes.Mean_stdDev_calc",
        )
        self.nodes["mean_stddev_calc_1"].set_plug_value(
            "parametric_maps", traits.Undefined
        )
        self.nodes[
            "mean_stddev_calc_1"
        ].process.parametric_maps = traits.Undefined
        self.add_process(
            "concat_to_list_1",
            "mia_processes.bricks.tools.tools.Concat_to_list",
        )
        self.add_process(
            "result_collector_1",
            "mia_processes.bricks.reports.processes.Result_collector",
        )
        self.nodes["result_collector_1"].set_plug_value(
            "parameter_files", traits.Undefined
        )
        self.nodes[
            "result_collector_1"
        ].process.parameter_files = traits.Undefined

        # links
        self.export_parameter(
            "find_in_list_1", "in_list", "spmT_images", is_optional=False
        )
        self.add_link("spmT_images->import_data_1.file_in_db")
        self.export_parameter(
            "find_in_list_2", "in_list", "beta_images", is_optional=False
        )
        self.export_parameter(
            "convroi_1", "convolve_with", "mask_002", is_optional=False
        )
        self.add_link("mask_002->resample_1_1.files_to_resample")
        self.add_link(
            "concat_to_list_of_list_1.listOflist->import_data_1.rois_list"
        )
        self.add_link("import_data_1.rois_files->convroi_1.images_to_convolve")
        self.add_link("find_in_list_1.out_file->files_to_list_1.file1")
        self.add_link("find_in_list_1.out_file->resample_1_1.reference_image")
        self.add_link("find_in_list_2.out_file->files_to_list_1.file2")
        self.add_link(
            "files_to_list_1.file_list->mean_stddev_calc_1.parametric_maps"
        )
        self.add_link("convroi_1.out_images->mean_stddev_calc_1.rois_files")
        self.add_link("convroi_1.out_images->resample_2_1.files_to_resample")
        self.export_parameter(
            "convroi_1", "out_images", "conv_roi_masks", is_optional=False
        )
        self.add_link("resample_1_1.out_files->resample_2_1.reference_image")
        self.export_parameter(
            "resample_2_1", "out_images", "resample2_masks", is_optional=False
        )
        self.add_link(
            "mean_stddev_calc_1.mean_out_files->concat_to_list_1.list1"
        )
        self.add_link(
            "mean_stddev_calc_1.std_out_files->concat_to_list_1.list2"
        )
        self.add_link(
            "concat_to_list_1.out_list->result_collector_1.parameter_files"
        )
        self.export_parameter(
            "result_collector_1", "out_files", "xls_files", is_optional=False
        )
        self.export_parameter(
            "result_collector_1", "patient_info", is_optional=True
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
            "concat_to_list_of_list_1": (
                -626.7365627381539,
                655.1783812050373,
            ),
            "import_data_1": (-359.6572976506295, 715.6133821215734),
            "find_in_list_1": (-508.4657057744808, 134.84299088223963),
            "find_in_list_2": (-361.3432851561171, 431.0012168458302),
            "files_to_list_1": (-151.78323445420767, 347.27763061397314),
            "convroi_1": (-180.2234309904108, 523.195153022265),
            "resample_1_1": (-218.19904973701512, 70.29956616326433),
            "resample_2_1": (264.02963177482127, 39.46765963517973),
            "inputs": (-623.9590860471251, 440.277979255905),
            "mean_stddev_calc_1": (39.318714641281616, 237.85728764977785),
            "concat_to_list_1": (312.67056990213575, 254.65280031900323),
            "result_collector_1": (203.0050153882558, 670.0575258538529),
            "outputs": (500.615875306018, 441.4937568209974),
        }

        # nodes dimensions
        self.node_dimension = {
            "concat_to_list_of_list_1": (157.21875, 110.0),
            "import_data_1": (147.703125, 180.0),
            "find_in_list_1": (118.84375, 110.0),
            "find_in_list_2": (118.84375, 110.0),
            "files_to_list_1": (97.640625, 145.0),
            "convroi_1": (208.0625, 145.0),
            "resample_1_1": (176.84375, 250.0),
            "resample_2_1": (193.921875, 145.0),
            "inputs": (119.17251223929541, 136.0),
            "mean_stddev_calc_1": (214.171875, 145.0),
            "concat_to_list_1": (110.71875, 110.0),
            "result_collector_1": (171.640625, 145.0),
            "outputs": (126.07117003946317, 111.0),
        }

        self.do_autoexport_nodes_parameters = False
