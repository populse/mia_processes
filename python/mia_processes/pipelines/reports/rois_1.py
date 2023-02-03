# -*- coding: utf-8 -*- #

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

# Capsul import
from capsul.api import Pipeline

# Other import
import traits.api as traits


class Rois_1(Pipeline):

    def pipeline_definition(self):
        # nodes
        self.add_process("find_in_list1",
                         "mia_processes.bricks.tools.Find_In_List")
        self.add_process("resample1",
                         "mia_processes.bricks.preprocess.others.Resample_1")
        self.nodes["resample1"].process.reference_image = traits.Undefined
        self.nodes["resample1"].process.files_to_resample = traits.Undefined
        self.add_process("roi_list_generator1",
                         "mia_processes.bricks.tools.Concat_to_list_of_list")
        self.add_process("conv_roi1",
                         "mia_processes.bricks.preprocess.others.Conv_ROI")
        self.nodes["conv_roi1"].process.in_image = traits.Undefined
        self.add_process("resample2",
                         "mia_processes.bricks.preprocess.others.Resample_2")
        self.add_process("files_to_list1",
                         "mia_processes.bricks.tools.tools.Files_To_List")
        self.add_process("find_in_list2",
                         "mia_processes.bricks.tools.Find_In_List")
        self.add_process(
                        "result_collector1",
                        "mia_processes.bricks.reports.Result_collector")
        self.add_process("mean_&_stdDev_calc1",
                         "mia_processes.bricks.reports.Mean_stdDev_calc")

        # links
        self.export_parameter("find_in_list1", "in_list", "spmT_images")
        self.export_parameter("resample1", "files_to_resample", "mask_002")
        self.add_link("mask_002->conv_roi1.in_image")
        self.export_parameter("find_in_list2", "in_list", "beta_images")
        self.add_link("find_in_list1.out_file->resample1.reference_image")
        self.add_link("find_in_list1.out_file->files_to_list1.file1")
        self.add_link("resample1.out_files->resample2.reference_image")
        self.add_link("roi_list_generator1.listOflist->resample2.doublet_list")
        self.add_link("roi_list_generator1.listOflist->conv_roi1.doublet_list")
        self.add_link("roi_list_generator1.listOflist->"
                      "result_collector1.doublet_list")
        self.add_link("roi_list_generator1.listOflist->"
                      "mean_&_stdDev_calc1.doublet_list")
        self.export_parameter("conv_roi1", "out_images", "conv_roi_masks")
        self.export_parameter("resample2", "out_images", "resample2_masks")
        self.add_link("files_to_list1.file_list->"
                      "mean_&_stdDev_calc1.parametric_maps")
        self.add_link("files_to_list1.file_list->result_collector1.parametric_maps")
        self.add_link("find_in_list2.out_file->files_to_list1.file2")
        self.export_parameter("result_collector1", "out_files", "xls_files")
        self.add_link("mean_&_stdDev_calc1.mean_out_files->"
                      "result_collector1.mean_in_files")
        self.add_link("mean_&_stdDev_calc1.std_out_files->"
                      "result_collector1.std_in_files")
        self.export_parameter("result_collector1", "patient_info",
                              is_optional=True)

        # nodes positions
        self.node_position = {
            "conv_roi1": (352.35910249999995, 425.5439),
            "find_in_list1": (187.4372275, 102.63739999999996),
            "inputs": (0.0, 164.05629999999996),
            "resample1": (341.58566499999995, 0.0),
            "find_in_list2": (187.4372275, 224.54389999999995),
            "result_collector1": (681.214865, 256.02479999999997),
            "resample2": (518.45744, 177.0811),
            "mean_&_stdDev_calc1": (506.1527525, 389.375),
            "outputs": (852.45024, 318.52479999999997),
            "roi_list_generator1": (162.710665, 346.4437),
            "files_to_list1": (356.53878999999995, 196.625),
        }

        self.do_autoexport_nodes_parameters = False
