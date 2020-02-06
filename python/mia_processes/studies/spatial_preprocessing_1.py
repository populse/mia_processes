# -*- coding: utf-8 -*- #

"""The pipelines library of the mia_processes package.

The objective of this module is to propose pipelines built from atomic bricks
or other pipelines provided in the mia_processes library.

:Contains:
    :Class:
        - Spatial_preprocessing_1

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# capsul import
from capsul.api import Pipeline


class Spatial_preprocessing_1(Pipeline):
    """
    *Data pre-processing for cerebrovascular reserve analysis (CVRa) at CLUNI - IRMaGe (Grenoble - France)*
    Please, see the complete documention for the `Spatial_preprocessing_1 brick in the populse.mia_processes web site:
    <https://populse.github.io/mia_processes/html/documentation/studies/Spatial_preprocessing_1.html>`_
    """

    def pipeline_definition(self):
        # nodes
        self.add_process("newsegment1", "mia_processes.preprocess.spm.spatial_preprocessing.NewSegment")
        self.add_process("realign1", "mia_processes.preprocess.spm.spatial_preprocessing.Realign")
        self.add_process("list_duplicate1", "mia_processes.tools.tools.List_Duplicate")
        self.add_process("normalize12_1", "mia_processes.preprocess.spm.spatial_preprocessing.Normalize12")
        self.add_process("normalize12_2", "mia_processes.preprocess.spm.spatial_preprocessing.Normalize12")
        self.add_process("smooth1", "mia_processes.preprocess.spm.spatial_preprocessing.Smooth")
        self.add_process("coregister1", "mia_processes.preprocess.spm.spatial_preprocessing.Coregister")

        # links
        self.export_parameter("list_duplicate1", "file_name", "anat_file")
        self.export_parameter("realign1", "in_files", "func_files")
        self.export_parameter("normalize12_2", "write_voxel_sizes", "voxel_sizes_func")
        self.export_parameter("newsegment1", "forward_deformation_field")
        self.add_link("newsegment1.forward_deformation_field->normalize12_1.deformation_file")
        self.add_link("newsegment1.forward_deformation_field->normalize12_2.deformation_file")
        self.export_parameter("newsegment1", "bias_field_images")
        self.export_parameter("newsegment1", "bias_corrected_images")
        self.export_parameter("newsegment1", "native_class_images")
        self.add_link("realign1.realigned_files->coregister1.apply_to_files")
        self.add_link("realign1.mean_image->coregister1.source")
        self.export_parameter("realign1", "realignment_parameters")
        self.add_link("list_duplicate1.out_file->coregister1.target")
        self.add_link("list_duplicate1.out_list->newsegment1.channel_files")
        self.add_link("list_duplicate1.out_list->normalize12_1.apply_to_files")
        self.export_parameter("normalize12_1", "normalized_files", "normalized_anat")
        self.add_link("normalize12_2.normalized_files->smooth1.in_files")
        self.export_parameter("smooth1", "smoothed_files", "smoothed_func")
        self.add_link("coregister1.coregistered_files->normalize12_2.apply_to_files")
        self.export_parameter("coregister1", "coregistered_source")

        # nodes positions
        self.node_position = {
            "outputs": (286.6713, 405.9444000000001),
            "inputs": (-747.6544000000002, 676.9191999999996),
            "normalize12_1": (-105.31994999999998, 361.1199999999999),
            "list_duplicate1": (-698.5690125000001, 322.7524000000001),
            "normalize12_2": (11.760050000000234, 814.1444),
            "newsegment1": (-529.162, 239.17460000000008),
            "realign1": (-572.1639500000001, 611.8702999999997),
            "coregister1": (-254.54776249999983, 717.4246000000004),
            "smooth1": (293.10992500000015, 837.1643999999998),
        }

        self.do_autoexport_nodes_parameters = False
