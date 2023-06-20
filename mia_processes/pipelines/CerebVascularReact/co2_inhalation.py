# -*- coding: utf-8 -*-

"""The pipelines library of the mia_processes package.

The CO2_inhalation pipeline is designed for CVR measurement
using hypercapnia challenge.

:Contains:
    :Class:
        - CO2_inhalation

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from capsul.api import Pipeline


class CO2_inhalation(Pipeline):
    """
    *CVR measurement under hypercapnia challenge*

    Please, see the complete documentation for the
    `CO2_inhalation pipeline in the populse.mia_processes website:
    <https://populse.github.io/mia_processes/html/documentation/pipelines/CerebVascularReact/CO2_inhalation.html>`_
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
        self.add_process(
            "3_boldStat",
            "mia_processes.pipelines.stat.bold_stat_cvr.Bold_stat_cvr",
        )
        self.nodes["3_boldStat"].process.nodes_activation = {
            "estimatecontrast": True,
            "level1design": True,
            "estimatemodel": True,
            "make_a_list": True,
        }
        self.add_process(
            "4_extract_roi_param",
            "mia_processes.pipelines.reports.extract_roi_param."
            "Extract_roi_param",
        )
        self.nodes["4_extract_roi_param"].process.nodes_activation = {
            "concat_to_list_of_list": True,
            "import_data": True,
            "find_in_list_1": True,
            "find_in_list_2": True,
            "files_to_list": True,
            "convroi": True,
            "resample1": True,
            "resample2": True,
            "mean_stddev_calc": True,
            "concat_to_list": True,
            "result_collector": True,
        }
        self.add_process(
            "files_to_list", "mia_processes.bricks.tools.tools.Files_To_List"
        )

        # links
        self.export_parameter(
            "1_spatial_preprocessing", "anat_file", is_optional=False
        )
        self.export_parameter(
            "1_spatial_preprocessing", "func_files", is_optional=False
        )
        self.export_parameter(
            "files_to_list", "file1", "regressor_physio", is_optional=False
        )
        self.export_parameter(
            "4_extract_roi_param", "patient_info", is_optional=True
        )
        self.export_parameter(
            "1_spatial_preprocessing", "bias_field_images", is_optional=True
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
            "1_spatial_preprocessing.realignment_parameters->"
            "files_to_list.file2"
        )
        self.export_parameter(
            "1_spatial_preprocessing", "normalized_anat", is_optional=True
        )
        self.add_link(
            "1_spatial_preprocessing.smoothed_func->3_boldStat.smoothed_func"
        )
        self.add_link(
            "1_spatial_preprocessing.smoothed_func->"
            "2_spatial_mask.smoothed_func"
        )
        self.export_parameter(
            "1_spatial_preprocessing", "coregistered_source", is_optional=False
        )
        self.add_link("2_spatial_mask.mask_002->4_extract_roi_param.mask_002")
        self.add_link("2_spatial_mask.mask_002->3_boldStat.mask_002")
        self.export_parameter("2_spatial_mask", "mask_003", is_optional=False)
        self.add_link(
            "3_boldStat.spmT_images->4_extract_roi_param.spmT_images"
        )
        self.export_parameter(
            "3_boldStat", "out_spm_mat_file", is_optional=False
        )
        self.add_link(
            "3_boldStat.beta_images->4_extract_roi_param.beta_images"
        )
        self.export_parameter(
            "4_extract_roi_param", "resample2_masks", is_optional=False
        )
        self.export_parameter(
            "4_extract_roi_param", "xls_files", is_optional=False
        )
        self.export_parameter(
            "4_extract_roi_param", "conv_roi_masks", is_optional=False
        )
        self.add_link("files_to_list.file_list->3_boldStat.regressors")

        # parameters order

        self.reorder_traits(
            (
                "anat_file",
                "func_files",
                "regressor_physio",
                "out_spm_mat_file",
                "bias_field_images",
                "normalized_anat",
                "coregistered_source",
                "mask_003",
                "conv_roi_masks",
                "patient_info",
                "resample2_masks",
                "xls_files",
            )
        )

        # nodes positions
        self.node_position = {
            "inputs": (-275.0790595529901, 183.49872617882536),
            "1_spatial_preprocessing": (
                -164.59507406061374,
                -222.31943258070515,
            ),
            "2_spatial_mask": (183.75989123102062, -212.8388047341995),
            "3_boldStat": (120.47624500640808, 271.95641405347675),
            "4_extract_roi_param": (528.7543856377755, -56.221297992058865),
            "files_to_list": (-52.58039247862696, 250.55123446892662),
            "outputs": (805.7712478465897, -179.76368248093885),
        }

        # nodes dimensions
        self.node_dimension = {
            "inputs": (138.46938723929543, 180.0),
            "1_spatial_preprocessing": (236.515625, 320.0),
            "2_spatial_mask": (200.15625, 145.0),
            "3_boldStat": (215.375, 145.0),
            "4_extract_roi_param": (206.875, 180.0),
            "files_to_list": (97.640625, 145.0),
            "outputs": (138.83679503946317, 320.0),
        }

        self.do_autoexport_nodes_parameters = False
