# -*- coding: utf-8 -*-

"""The pipelines library of the mia_processes package.

The objective of this module is to propose pipelines built from atomic bricks
or other pipelines provided in the mia_processes library.

:Contains:
    :Class:
        - Bold_spatial_preprocessing1

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

# capsul import
from capsul.api import Pipeline


class Bold_spatial_preprocessing1(Pipeline):
    """
    *Data pre-processing for cerebrovascular reserve analysis (CVR)*

    Please, see the complete documentation for the
    `Bold_spatial_preprocessing1 brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Bold_spatial_preprocessing1.html>`_
    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "newsegment",
            "mia_processes.bricks.preprocess.spm."
            "spatial_preprocessing.NewSegment",
        )
        self.nodes["newsegment"].process.channel_files = traits.Undefined

        self.add_process(
            "realign",
            "mia_processes.bricks.preprocess.spm.spatial_preprocessing."
            "Realign",
        )

        self.add_process(
            "list_duplicate",
            "mia_processes.bricks.tools.tools.List_Duplicate",
        )

        self.add_process(
            "normalize12_1",
            "mia_processes.bricks.preprocess.spm.spatial_preprocessing."
            "Normalize12",
        )
        self.nodes["normalize12_1"].process.apply_to_files = traits.Undefined

        self.add_process(
            "normalize12_2",
            "mia_processes.bricks.preprocess.spm.spatial_preprocessing."
            "Normalize12",
        )
        self.nodes["normalize12_2"].process.apply_to_files = traits.Undefined
        self.nodes["normalize12_2"].process.write_voxel_sizes = [2.0, 2.0, 2.0]

        self.add_process(
            "smooth",
            "mia_processes.bricks.preprocess.spm.spatial_preprocessing.Smooth",
        )
        self.nodes["smooth"].process.in_files = traits.Undefined

        self.add_process(
            "coregister",
            "mia_processes.bricks.preprocess.spm.spatial_preprocessing."
            "Coregister",
        )
        self.nodes["coregister"].process.source = traits.Undefined
        self.nodes["coregister"].process.apply_to_files = traits.Undefined

        # links
        self.export_parameter("list_duplicate", "file_name", "anat_file")
        self.export_parameter("realign", "in_files", "func_files")
        self.export_parameter("newsegment", "bias_corrected_images")
        self.export_parameter("newsegment", "bias_field_images")
        self.export_parameter("newsegment", "native_class_images")
        self.export_parameter("newsegment", "forward_deformation_field")
        self.add_link(
            "newsegment.forward_deformation_field->"
            "normalize12_2.deformation_file"
        )
        self.add_link(
            "newsegment.forward_deformation_field->"
            "normalize12_1.deformation_file"
        )
        self.add_link("realign.realigned_files->coregister.apply_to_files")
        self.add_link("realign.mean_image->coregister.source")
        self.export_parameter("realign", "realignment_parameters")
        self.add_link("list_duplicate.out_file->coregister.target")
        self.add_link("list_duplicate.out_list->normalize12_1.apply_to_files")
        self.add_link("list_duplicate.out_list->newsegment.channel_files")
        self.export_parameter(
            "normalize12_1", "normalized_files", "normalized_anat"
        )
        self.add_link("normalize12_2.normalized_files->smooth.in_files")
        self.export_parameter("smooth", "smoothed_files", "smoothed_func")
        self.export_parameter("coregister", "coregistered_source")
        self.add_link(
            "coregister.coregistered_files->normalize12_2.apply_to_files"
        )

        # nodes positions
        self.node_position = {
            "outputs": (431.5335139107567, 281.77678807649437),
            "inputs": (-703.3088243130337, 598.5753496196925),
            "normalize12_1": (20.32584777973804, 161.56490940865137),
            "list_duplicate": (-698.5690125000001, 322.7524000000001),
            "normalize12_2": (177.31686589800802, 797.8843555814453),
            "newsegment": (-453.77452133215724, 199.26358188173032),
            "realign": (-501.21102890085405, 716.821495792487),
            "coregister": (-115.59829201417193, 792.8120786678431),
            "smooth": (508.92506000990295, 860.8153736997155),
        }

        # nodes dimensions
        self.node_dimension = {
            "inputs": (279.359375, 90.0),
            "newsegment": (332.46875, 355.0),
            "realign": (275.21875, 530.0),
            "list_duplicate": (129.84375, 110.0),
            "normalize12_1": (291.421875, 530.0),
            "normalize12_2": (291.421875, 530.0),
            "smooth": (218.046875, 215.0),
            "coregister": (224.171875, 460.0),
            "outputs": (239.515625, 320.0),
        }

        self.do_autoexport_nodes_parameters = False
