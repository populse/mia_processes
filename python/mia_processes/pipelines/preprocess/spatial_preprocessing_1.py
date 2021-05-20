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

# other import
import traits.api as traits


class Spatial_preprocessing_1(Pipeline):
    """
    *Data pre-processing for cerebrovascular reserve analysis (CVRa) at CLUNI - IRMaGe (Grenoble - France)*
    Please, see the complete documention for the `Spatial_preprocessing_1 brick in the populse.mia_processes web site:
    <https://populse.github.io/mia_processes/documentation/pipelines/preprocess/Spatial_preprocessing_1.html>`_
    """

    def pipeline_definition(self):
        # nodes
        self.add_process("newsegment1", "mia_processes.bricks.preprocess.spm.spatial_preprocessing.NewSegment")
        self.nodes["newsegment1"].process.channel_files = traits.Undefined

        self.add_process("realign1", "mia_processes.bricks.preprocess.spm.spatial_preprocessing.Realign")

        self.add_process("list_duplicate1", "mia_processes.bricks.tools.tools.List_Duplicate")

        self.add_process("normalize12_1", "mia_processes.bricks.preprocess.spm.spatial_preprocessing.Normalize12")
        self.nodes["normalize12_1"].process.apply_to_files = traits.Undefined

        self.add_process("normalize12_2", "mia_processes.bricks.preprocess.spm.spatial_preprocessing.Normalize12")
        self.nodes["normalize12_2"].process.apply_to_files = traits.Undefined
        self.nodes["normalize12_2"].process.write_voxel_sizes = [2.0, 2.0, 2.0]

        self.add_process("smooth1", "mia_processes.bricks.preprocess.spm.spatial_preprocessing.Smooth")
        self.nodes["smooth1"].process.in_files = traits.Undefined

        self.add_process("coregister1", "mia_processes.bricks.preprocess.spm.spatial_preprocessing.Coregister")
        self.nodes["coregister1"].process.source = traits.Undefined
        self.nodes["coregister1"].process.apply_to_files = traits.Undefined


        # links
        self.export_parameter("list_duplicate1", "file_name", "anat_file")
        self.export_parameter("realign1", "in_files", "func_files")
        self.export_parameter("newsegment1", "bias_corrected_images")
        self.export_parameter("newsegment1", "bias_field_images")
        self.export_parameter("newsegment1", "native_class_images")
        self.export_parameter("newsegment1", "forward_deformation_field")
        self.add_link("newsegment1.forward_deformation_field->normalize12_2.deformation_file")
        self.add_link("newsegment1.forward_deformation_field->normalize12_1.deformation_file")
        self.add_link("realign1.realigned_files->coregister1.apply_to_files")
        self.add_link("realign1.mean_image->coregister1.source")
        self.export_parameter("realign1", "realignment_parameters")
        self.add_link("list_duplicate1.out_file->coregister1.target")
        self.add_link("list_duplicate1.out_list->normalize12_1.apply_to_files")
        self.add_link("list_duplicate1.out_list->newsegment1.channel_files")
        self.export_parameter("normalize12_1", "normalized_files", "normalized_anat")
        self.add_link("normalize12_2.normalized_files->smooth1.in_files")
        self.export_parameter("smooth1", "smoothed_files", "smoothed_func")
        self.export_parameter("coregister1", "coregistered_source")
        self.add_link("coregister1.coregistered_files->normalize12_2.apply_to_files")


        # nodes positions
        self.node_position = {
            "outputs": (431.5335139107567, 281.77678807649437),
            "inputs": (-703.3088243130337, 598.5753496196925),
            "normalize12_1": (20.32584777973804, 161.56490940865137),
            "list_duplicate1": (-698.5690125000001, 322.7524000000001),
            "normalize12_2": (177.31686589800802, 797.8843555814453),
            "newsegment1": (-453.77452133215724, 199.26358188173032),
            "realign1": (-501.21102890085405, 716.821495792487),
            "coregister1": (-115.59829201417193, 792.8120786678431),
            "smooth1": (508.92506000990295, 860.8153736997155),
        }


        # nodes dimensions
        self.node_dimension = {
            "inputs": (279.359375, 90.0),
            "newsegment1": (332.46875, 355.0),
            "realign1": (275.21875, 530.0),
            "list_duplicate1": (129.84375, 110.0),
            "normalize12_1": (291.421875, 530.0),
            "normalize12_2": (291.421875, 530.0),
            "smooth1": (218.046875, 215.0),
            "coregister1": (224.171875, 460.0),
            "outputs": (239.515625, 320.0),
        }


        self.do_autoexport_nodes_parameters = False
