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
    `CO2_inhalation pipeline in the mia_processes website:
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
        self.add_process(
            "reportco2inhalcvr",
            "mia_processes.bricks.reports.reporting.ReportCO2inhalCvr",
        )
        self.nodes["reportco2inhalcvr"].process.beta_vmin = 0.01
        self.nodes["reportco2inhalcvr"].process.beta_vmax = 0.25
        self.nodes["reportco2inhalcvr"].process.spmT_vmin = 3.0
        self.nodes["reportco2inhalcvr"].process.spmT_vmax = 12
        self.add_process(
            "list_to_file_1", "mia_processes.bricks.tools.tools.List_To_File"
        )
        self.add_process(
            "make_cvr_reg_physio_1",
            "mia_processes.bricks.tools.tools." "Make_CVR_reg_physio",
        )

        # links
        self.export_parameter(
            "1_spatial_preprocessing", "anat_file", is_optional=False
        )
        self.export_parameter(
            "1_spatial_preprocessing", "func_files", is_optional=False
        )
        self.add_link("func_files->make_cvr_reg_physio_1.func_file")
        self.export_parameter(
            "reportco2inhalcvr", "patient_info", is_optional=True
        )
        self.add_link("patient_info->4_extract_roi_param.patient_info")
        self.export_parameter(
            "make_cvr_reg_physio_1", "trigger_data", is_optional=True
        )
        self.export_parameter(
            "make_cvr_reg_physio_1", "physio_data", is_optional=True
        )
        self.nodes["1_spatial_preprocessing"].process.trait(
            "bias_field_images"
        ).userlevel = 1
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
        self.add_link(
            "1_spatial_preprocessing.realignment_parameters->"
            "reportco2inhalcvr.realignment_parameters"
        )
        self.add_link(
            "1_spatial_preprocessing.normalized_anat->"
            "reportco2inhalcvr.norm_anat"
        )
        self.add_link(
            "1_spatial_preprocessing.smoothed_func->3_boldStat.smoothed_func"
        )
        self.add_link(
            "1_spatial_preprocessing.smoothed_func->"
            "2_spatial_mask.smoothed_func"
        )
        self.nodes["1_spatial_preprocessing"].process.trait(
            "coregistered_source"
        ).userlevel = 1
        self.export_parameter(
            "1_spatial_preprocessing", "coregistered_source", is_optional=False
        )
        self.add_link(
            "1_spatial_preprocessing.normalized_func->"
            "reportco2inhalcvr.norm_func"
        )
        self.add_link("2_spatial_mask.mask_002->4_extract_roi_param.mask_002")
        self.add_link("2_spatial_mask.mask_002->3_boldStat.mask_002")
        self.nodes["2_spatial_mask"].process.trait("mask_003").userlevel = 1
        self.export_parameter("2_spatial_mask", "mask_003", is_optional=False)
        self.add_link(
            "3_boldStat.spmT_images->4_extract_roi_param.spmT_images"
        )
        self.nodes["3_boldStat"].process.trait(
            "out_spm_mat_file"
        ).userlevel = 1
        self.export_parameter(
            "3_boldStat", "out_spm_mat_file", is_optional=False
        )
        self.add_link("3_boldStat.beta_images->list_to_file_1.file_list")
        self.add_link(
            "3_boldStat.beta_images->4_extract_roi_param.beta_images"
        )
        self.add_link("3_boldStat.spmT_images->reportco2inhalcvr.spmT_image")
        self.nodes["4_extract_roi_param"].process.trait(
            "resample2_masks"
        ).userlevel = 1
        self.export_parameter(
            "4_extract_roi_param", "resample2_masks", is_optional=False
        )
        self.export_parameter(
            "4_extract_roi_param", "xls_files", is_optional=False
        )
        self.nodes["4_extract_roi_param"].process.trait(
            "conv_roi_masks"
        ).userlevel = 1
        self.export_parameter(
            "4_extract_roi_param", "conv_roi_masks", is_optional=False
        )
        self.add_link("files_to_list.file_list->3_boldStat.regressors")
        self.export_parameter("reportco2inhalcvr", "report", is_optional=True)
        self.add_link("list_to_file_1.file->reportco2inhalcvr.beta_image")
        self.add_link("make_cvr_reg_physio_1.cvr_reg->files_to_list.file1")
        self.add_link(
            "make_cvr_reg_physio_1.cvr_reg->"
            "reportco2inhalcvr.regressor_physio"
        )

        # parameters order
        self.reorder_traits(
            (
                "anat_file",
                "func_files",
                "trigger_data",
                "physio_data",
                "out_spm_mat_file",
                "bias_field_images",
                "coregistered_source",
                "mask_003",
                "conv_roi_masks",
                "patient_info",
                "resample2_masks",
                "xls_files",
                "report",
            )
        )

        # nodes positions
        self.node_position = {
            "inputs": (-261.5363166564018, 251.21244066176695),
            "1_spatial_preprocessing": (
                -221.77554406843106,
                -205.76719126265274,
            ),
            "2_spatial_mask": (123.5699228017393, -98.47786471856485),
            "3_boldStat": (382.3026076737819, 162.10972167003825),
            "4_extract_roi_param": (581.4206080133966, -206.69621906526237),
            "files_to_list": (198.71272571362283, 330.3029426377244),
            "reportco2inhalcvr": (748.1071959262655, 258.2245925859173),
            "outputs": (971.2936610271133, -146.65919984483403),
            "list_to_file_1": (455.4447133457113, 554.0646315898613),
            "make_cvr_reg_physio_1": (-61.59780851025405, 457.28542588968276),
        }

        # nodes dimensions
        self.node_dimension = {
            "inputs": (141.46938723929543, 180.0),
            "1_spatial_preprocessing": (236.515625, 355.0),
            "2_spatial_mask": (200.15625, 145.0),
            "3_boldStat": (215.375, 145.0),
            "4_extract_roi_param": (206.875, 180.0),
            "files_to_list": (97.640625, 145.0),
            "outputs": (138.83679503946317, 320.0),
            "reportco2inhalcvr": (212.484375, 495.0),
            "list_to_file_1": (117.75, 110.0),
            "make_cvr_reg_physio_1": (174.125, 145.0),
        }

        self.do_autoexport_nodes_parameters = False
