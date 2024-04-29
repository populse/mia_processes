# -*- coding: utf-8 -*-

"""The pipelines library of the mia_processes package.

The GE2REc pipeline

:Contains:
    :Class:
        - Ge2rec

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from capsul.api import Pipeline


class Ge2rec(Pipeline):
    """
    *GE2REC pipelinee*

    Please, see the complete documentation for the
    `Ge2rec pipeline in the mia_processes website:
    <https://populse.github.io/mia_processes/html/documentation/pipelines/FunctionallangageMemory/Ge2rec.html>`_
    """

    def pipeline_definition(self):
        """Building the pipeline"""

        # nodes
        # PREPROCESSING
        self.add_process(
            "files_to_list", "mia_processes.bricks.tools.tools.Files_To_List"
        )
        self.add_process(
            "1_preprocessing",
            "mia_processes.pipelines.preprocess.bold_spatial_preprocessing1."
            "Bold_spatial_preprocessing1",
        )
        self.nodes["1_preprocessing"].process.nodes_activation = {
            "newsegment": True,
            "realign": True,
            "list_duplicate": True,
            "normalize12_1": True,
            "normalize12_2": True,
            "smooth": True,
            "coregister": True,
        }
        self.nodes["1_preprocessing"].process.nodes[
            "normalize12_1"
        ].process.write_voxel_sizes = [1.0, 1.0, 1.0]
        self.nodes["1_preprocessing"].process.nodes[
            "normalize12_2"
        ].process.write_voxel_sizes = [3.0, 3.0, 3.0]
        self.nodes["1_preprocessing"].process.trait(
            "coregistered_source"
        ).userlevel = 1
        self.export_parameter(
            "1_preprocessing", "coregistered_source", is_optional=False
        )

        self.add_process(
            "filter_files_list_1",
            "mia_processes.bricks.tools.tools.Filter_Files_List",
        )
        self.add_process(
            "filter_files_list_2",
            "mia_processes.bricks.tools.tools.Filter_Files_List",
        )
        self.add_process(
            "filter_files_list_3",
            "mia_processes.bricks.tools.tools.Filter_Files_List",
        )
        self.nodes["filter_files_list_3"].process.index_filter = [2]
        self.add_process(
            "filter_files_list_4",
            "mia_processes.bricks.tools.tools.Filter_Files_List",
        )
        self.nodes["filter_files_list_4"].process.index_filter = [2]
        self.add_process(
            "filter_files_list_5",
            "mia_processes.bricks.tools.tools.Filter_Files_List",
        )
        self.nodes["filter_files_list_5"].process.index_filter = [3]
        self.add_process(
            "filter_files_list_6",
            "mia_processes.bricks.tools.tools.Filter_Files_List",
        )
        self.nodes["filter_files_list_6"].process.index_filter = [3]
        self.add_process(
            "make_a_list_1", "mia_processes.bricks.tools.tools.Make_A_List"
        )
        self.add_process(
            "make_a_list_2", "mia_processes.bricks.tools.tools.Make_A_List"
        )
        self.add_process(
            "make_a_list_3", "mia_processes.bricks.tools.tools.Make_A_List"
        )

        # STATS gene
        self.add_process(
            "2a_gene_level1design",
            "mia_processes.bricks.stat.spm.model.Level1Design",
        )
        self.nodes["2a_gene_level1design"].process.out_dir_name = "stats_gene"
        self.nodes["2a_gene_level1design"].process.timing_units = "secs"
        self.nodes["2a_gene_level1design"].process.sess_cond_names = [
            ["GENE", "CONTROL"]
        ]
        self.nodes["2a_gene_level1design"].process.sess_cond_onsets = [
            [
                [0.0, 90.0, 180.0, 270.0, 360.0],
                [50.0, 140.0, 230.0, 320.0, 410.0],
            ]
        ]
        self.nodes["2a_gene_level1design"].process.sess_cond_durations = [
            [[40.0, 40.0, 40.0, 40.0, 40.0], [40.0, 40.0, 40.0, 40.0, 40.0]]
        ]
        self.nodes["2a_gene_level1design"].process.sess_cond_tmod = [[0, 0]]
        self.nodes["2a_gene_level1design"].process.sess_cond_orth = [[1, 1]]
        self.nodes["2a_gene_level1design"].process.sess_hpf = [128.0]
        self.add_process(
            "2b_gene_estimateModel",
            "mia_processes.bricks.stat.spm.model." "EstimateModel",
        )
        self.add_process(
            "2c_gene_estimateContrast",
            "mia_processes.bricks.stat.spm.model." "EstimateContrast",
        )
        self.nodes["2c_gene_estimateContrast"].process.T_contrast_names = [
            "GENE-CONTROL",
            "GENE",
        ]
        self.nodes["2c_gene_estimateContrast"].process.T_condition_names = [
            ["GENE", "CONTROL"],
            ["GENE"],
        ]
        self.nodes["2c_gene_estimateContrast"].process.T_contrast_weights = [
            [1.0, -1.0],
            [2.0],
        ]

        # STATS reco
        self.add_process(
            "3a_reco_level1design",
            "mia_processes.bricks.stat.spm.model.Level1Design",
        )
        self.nodes["3a_reco_level1design"].process.out_dir_name = "stats_reco"
        self.nodes["3a_reco_level1design"].process.timing_units = "secs"
        self.nodes["3a_reco_level1design"].process.sess_cond_names = [
            ["OLD", "CONTROL", "NEW"]
        ]
        self.nodes["3a_reco_level1design"].process.sess_cond_onsets = [
            [
                [
                    5.0,
                    7.5,
                    12.5,
                    20.0,
                    35.0,
                    40.0,
                    52.5,
                    57.5,
                    62.5,
                    75.0,
                    82.5,
                    92.5,
                    95.0,
                    107.5,
                    115.0,
                    125.0,
                    130.0,
                    132.5,
                    142.5,
                    152.5,
                    155.0,
                    187.5,
                    190.0,
                    210.0,
                    242.5,
                    247.5,
                    255.0,
                    257.5,
                    275.0,
                    277.5,
                    285.0,
                    335.0,
                    337.5,
                    342.5,
                    350.0,
                    352.5,
                    365.0,
                    385.0,
                    387.5,
                    397.0,
                ],
                [
                    2.5,
                    25.0,
                    27.5,
                    47.5,
                    50.0,
                    72.5,
                    77.5,
                    80.0,
                    105.0,
                    110.0,
                    120.0,
                    150.0,
                    160.0,
                    162.5,
                    172.5,
                    180.0,
                    182.5,
                    202.5,
                    212.5,
                    217.5,
                    225.0,
                    235.0,
                    237.5,
                    252.5,
                    267.5,
                    287.5,
                    290.0,
                    302.5,
                    307.5,
                    312.5,
                    315.0,
                    322.5,
                    327.5,
                    347.5,
                    355.0,
                    360.0,
                    372.5,
                    375.0,
                    380.0,
                    382.5,
                ],
                [
                    0.0,
                    10.0,
                    17.5,
                    22.5,
                    37.5,
                    42.5,
                    45.0,
                    65.0,
                    85.0,
                    90.0,
                    102.5,
                    117.5,
                    122.5,
                    127.5,
                    137.5,
                    140.0,
                    147.5,
                    157.5,
                    170.0,
                    175.0,
                    185.0,
                    197.5,
                    222.5,
                    230.0,
                    232.5,
                    245.0,
                    265.0,
                    270.0,
                    272.5,
                    280.0,
                    282.5,
                    295.0,
                    300.0,
                    305.0,
                    317.5,
                    320.0,
                    330.0,
                    345.0,
                    357.5,
                    362.5,
                ],
            ]
        ]
        self.nodes["3a_reco_level1design"].process.sess_cond_durations = [
            [[0.0], [0.0], [0.0]]
        ]
        self.nodes["3a_reco_level1design"].process.sess_cond_tmod = [[0, 0, 0]]
        self.nodes["3a_reco_level1design"].process.sess_cond_orth = [[1, 1, 1]]
        self.nodes["3a_reco_level1design"].process.sess_hpf = [128.0]
        self.add_process(
            "3b_reco_estimateModel",
            "mia_processes.bricks.stat.spm.model." "EstimateModel",
        )
        self.add_process(
            "3c_reco_estimateContrast",
            "mia_processes.bricks.stat.spm.model." "EstimateContrast",
        )
        self.nodes["3c_reco_estimateContrast"].process.T_contrast_names = [
            "OLD-CONTROL",
            "OLD",
        ]
        self.nodes["3c_reco_estimateContrast"].process.T_condition_names = [
            ["OLD", "CONTROL"],
            ["OLD"],
        ]
        self.nodes["3c_reco_estimateContrast"].process.T_contrast_weights = [
            [1.0, -1.0],
            [1.0],
        ]

        # STATS RECALL
        self.add_process(
            "4a_recall_level1design",
            "mia_processes.bricks.stat.spm.model.Level1Design",
        )
        self.nodes["4a_recall_level1design"].process.out_dir_name = (
            "stats_rappel"
        )
        self.nodes["4a_recall_level1design"].process.timing_units = "secs"
        self.nodes["4a_recall_level1design"].process.sess_cond_names = [
            ["RAPPEL"]
        ]
        self.nodes["4a_recall_level1design"].process.sess_cond_onsets = [
            [
                [0.0, 50.0, 100.0, 150.0, 200.0],
                [40.0, 90.0, 140.0, 190.0, 240.0],
            ]
        ]
        self.nodes["4a_recall_level1design"].process.sess_cond_durations = [
            [[40.0, 40.0, 40.0, 40.0, 40.0], [10.0, 10.0, 10.0, 10.0, 10.0]]
        ]
        self.nodes["4a_recall_level1design"].process.sess_cond_orth = [[1]]
        self.nodes["4a_recall_level1design"].process.sess_hpf = [128.0]
        self.add_process(
            "4b_recall_estimateModel",
            "mia_processes.bricks.stat.spm.model.EstimateModel",
        )
        self.add_process(
            "4c_recall_estimateContrast",
            "mia_processes.bricks.stat.spm.model.EstimateContrast",
        )
        self.nodes["4c_recall_estimateContrast"].process.T_contrast_names = [
            "RAPPEL"
        ]
        self.nodes["4c_recall_estimateContrast"].process.T_condition_names = [
            ["RAPPEL"]
        ]

        # links
        self.export_parameter(
            "1_preprocessing", "anat_file", "anat_file", is_optional=False
        )
        self.export_parameter(
            "files_to_list", "file1", "func_gene_file", is_optional=False
        )
        self.export_parameter(
            "files_to_list", "file2", "func_reco_file", is_optional=False
        )
        self.export_parameter(
            "files_to_list", "file3", "func_recall_file", is_optional=False
        )
        self.add_link("files_to_list.file_list->1_preprocessing.func_files")
        self.add_link(
            "1_preprocessing.smoothed_func->filter_files_list_1.in_list"
        )
        self.add_link(
            "1_preprocessing.realignment_parameters->"
            "filter_files_list_2.in_list"
        )
        self.add_link(
            "1_preprocessing.smoothed_func->filter_files_list_3.in_list"
        )
        self.add_link(
            "1_preprocessing.realignment_parameters->"
            "filter_files_list_4.in_list"
        )
        self.add_link(
            "1_preprocessing.smoothed_func->filter_files_list_5.in_list"
        )
        self.add_link(
            "1_preprocessing.realignment_parameters->"
            "filter_files_list_6.in_list"
        )
        self.add_link(
            "filter_files_list_1.filtered_list->"
            "2a_gene_level1design.sess_scans"
        )
        self.add_link("filter_files_list_2.filtered_list->make_a_list_1.obj1")
        self.add_link(
            "make_a_list_1.obj_list->2a_gene_level1design.sess_multi_reg"
        )
        self.add_link(
            "filter_files_list_3.filtered_list->"
            "3a_reco_level1design.sess_scans"
        )
        self.add_link("filter_files_list_4.filtered_list->make_a_list_2.obj1")
        self.add_link(
            "make_a_list_2.obj_list->3a_reco_level1design.sess_multi_reg"
        )
        self.add_link(
            "filter_files_list_5.filtered_list->"
            "4a_recall_level1design.sess_scans"
        )
        self.add_link("filter_files_list_6.filtered_list->make_a_list_3.obj1")
        self.add_link(
            "make_a_list_3.obj_list->4a_recall_level1design.sess_multi_reg"
        )
        self.export_parameter(
            "1_preprocessing", "normalized_anat", is_optional=True
        )

        self.add_link(
            "2a_gene_level1design.spm_mat_file->"
            "2b_gene_estimateModel.spm_mat_file"
        )
        self.add_link(
            "2b_gene_estimateModel.out_spm_mat_file->"
            "2c_gene_estimateContrast.spm_mat_file"
        )
        self.add_link(
            "2b_gene_estimateModel.beta_images->"
            "2c_gene_estimateContrast.beta_images"
        )
        self.add_link(
            "2b_gene_estimateModel.residual_image->"
            "2c_gene_estimateContrast.residual_image"
        )
        self.export_parameter(
            "2c_gene_estimateContrast",
            "out_spm_mat_file",
            "out_spm_mat_file_gene",
            is_optional=False,
        )

        self.add_link(
            "3a_reco_level1design.spm_mat_file->"
            "3b_reco_estimateModel.spm_mat_file"
        )
        self.add_link(
            "3b_reco_estimateModel.out_spm_mat_file->"
            "3c_reco_estimateContrast.spm_mat_file"
        )
        self.add_link(
            "3b_reco_estimateModel.beta_images->"
            "3c_reco_estimateContrast.beta_images"
        )
        self.add_link(
            "3b_reco_estimateModel.residual_image->"
            "3c_reco_estimateContrast.residual_image"
        )
        self.export_parameter(
            "3c_reco_estimateContrast",
            "out_spm_mat_file",
            "out_spm_mat_file_reco",
            is_optional=False,
        )

        self.add_link(
            "4a_recall_level1design.spm_mat_file->"
            "4b_recall_estimateModel.spm_mat_file"
        )
        self.add_link(
            "4b_recall_estimateModel.out_spm_mat_file->"
            "4c_recall_estimateContrast.spm_mat_file"
        )
        self.add_link(
            "4b_recall_estimateModel.beta_images->"
            "4c_recall_estimateContrast.beta_images"
        )
        self.add_link(
            "4b_recall_estimateModel.residual_image->"
            "4c_recall_estimateContrast.residual_image"
        )
        self.export_parameter(
            "4c_recall_estimateContrast",
            "out_spm_mat_file",
            "out_spm_mat_file_rappel",
            is_optional=False,
        )

        # parameters order
        self.reorder_traits(
            (
                "anat_file",
                "func_gene_file",
                "func_reco_file",
                "func_recall_file",
                "out_spm_mat_file_gene",
                "out_spm_mat_file_reco",
                "out_spm_mat_file_rappel",
            )
        )

        # nodes positions
        self.node_position = {
            "inputs": (-900.0, -260),
            "files_to_list": (-686, -158),
            "1_preprocessing": (-700.0, -314),
            "filter_files_list_1": (-65, -671.8044843685886),
            "filter_files_list_2": (-65, -302.38649739872164),
            "filter_files_list_3": (-65, 260.67801499889794),
            "filter_files_list_4": (-65, 647.9710658544036),
            "filter_files_list_5": (-65, 1185.7125710807024),
            "filter_files_list_6": (-465, 1580.4535652218892),
            "make_a_list_1": (312.31121725366285, -280.0518231067058),
            "make_a_list_2": (338.23780514863375, 418.7650015170469),
            "make_a_list_3": (239.47112906454464, 2239.654481700065),
            "2a_gene_level1design": (30, -864),
            "2b_gene_estimateModel": (425, -864),
            "2c_gene_estimateContrast": (800, -864),
            "3a_reco_level1design": (30, 37),
            "3b_reco_estimateModel": (425, 37),
            "3c_reco_estimateContrast": (800, 37),
            "4a_recall_level1design": (30, 937),
            "4b_recall_estimateModel": (425, 937),
            "4c_recall_estimateContrast": (800, 937),
            "outputs": (1300, 37),
        }

        # nodes dimensions
        self.node_dimension = {
            "inputs": (152.5318872392954, 140.0),
            "1_preprocessing": (294.734375, 320.0),
            "2_preprocessing_reco": (384.59375, 145.0),
            "3_preprocessing_recall": (384.59375, 145.0),
            "2a_gene_level1design": (349.03125, 880.0),
            "2b_gene_estimateModel": (297.53125, 390.0),
            "2c_gene_estimateContrast": (304.0625, 425.0),
            "3a_reco_level1design": (349.03125, 880.0),
            "3b_reco_estimateModel": (297.53125, 390.0),
            "3c_reco_estimateContrast": (304.0625, 425.0),
            "4a_recall_level1design": (349.03125, 880.0),
            "4b_recall_estimateModel": (297.53125, 390.0),
            "4c_recall_estimateContrast": (304.0625, 425.0),
            "files_to_list": (118.8125, 145.0),
            "filter_files_list_1": (202.390625, 110.0),
            "filter_files_list_2": (199.390625, 110.0),
            "filter_files_list_3": (199.390625, 110.0),
            "filter_files_list_4": (199.390625, 110.0),
            "filter_files_list_5": (199.390625, 110.0),
            "filter_files_list_6": (199.390625, 110.0),
            "make_a_list_1": (114.890625, 110.0),
            "make_a_list_2": (114.890625, 110.0),
            "make_a_list_3": (114.890625, 110.0),
            "outputs": (205.10242003946317, 115.0),
        }

        self.do_autoexport_nodes_parameters = False
