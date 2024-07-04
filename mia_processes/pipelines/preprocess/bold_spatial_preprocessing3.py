# -*- coding: utf-8 -*-

"""The pipelines library of the mia_processes package.

The objective of this module is to propose pipelines built from atomic bricks
or other pipelines provided in the mia_processes library.

:Contains:
    :Class:
        - Bold_spatial_preprocessing3

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################
import os

from capsul.api import Pipeline
from populse_mia.software_properties import Config


class Bold_spatial_preprocessing3(Pipeline):
    """
    **

    Please, see the complete documentation for the
    `Bold_spatial_preprocessing3 brick in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Bold_spatial_preprocessing3.html>`_
    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "newsegment",
            "mia_processes.bricks.preprocess.spm."
            "spatial_preprocessing.NewSegment",
        )
        self.nodes["newsegment"].process.channel_info = (
            0.001,
            60.0,
            (False, True),
        )

        config = Config()

        tpm_path = os.path.join(
            config.get_resources_path(), "spm12", "tpm", "TPM.nii"
        )
        tissues_list = [
            ((tpm_path, 1), 1, (True, False), (False, False)),
            ((tpm_path, 2), 1, (True, False), (False, False)),
            ((tpm_path, 3), 2, (True, False), (False, False)),
            ((tpm_path, 4), 3, (True, False), (False, False)),
            ((tpm_path, 5), 4, (True, False), (False, False)),
            ((tpm_path, 6), 2, (True, False), (False, False)),
        ]
        self.nodes["newsegment"].process.tissues = tissues_list

        self.add_process(
            "realign",
            "mia_processes.bricks.preprocess.spm."
            "spatial_preprocessing.Realign",
        )
        self.add_process(
            "normalize12_1",
            "mia_processes.bricks.preprocess.spm."
            "spatial_preprocessing.Normalize12",
        )
        self.nodes["normalize12_1"].process.write_voxel_sizes = [1.0, 1.0, 1.0]
        self.add_process(
            "smooth",
            "mia_processes.bricks.preprocess.spm."
            "spatial_preprocessing.Smooth",
        )
        self.add_process(
            "coregister",
            "mia_processes.bricks.preprocess.spm."
            "spatial_preprocessing.Coregister",
        )
        self.add_process(
            "slicetiming",
            "mia_processes.bricks.preprocess.spm."
            "spatial_preprocessing.SliceTiming",
        )
        self.add_process(
            "normalize12_2",
            "mia_processes.bricks.preprocess.spm."
            "spatial_preprocessing.Normalize12",
        )
        self.nodes["normalize12_2"].process.write_voxel_sizes = [3.0, 3.0, 3.0]

        # links
        self.export_parameter(
            "coregister", "source", "anat_file", is_optional=False
        )
        self.export_parameter(
            "realign", "in_files", "func_files", is_optional=False
        )
        self.export_parameter(
            "slicetiming", "acquisition", "st_acquisition", is_optional=True
        )
        self.export_parameter(
            "normalize12_2",
            "write_voxel_sizes",
            "func_out_voxel_sizes",
            is_optional=True,
        )
        self.export_parameter(
            "normalize12_1",
            "write_voxel_sizes",
            "anat_out_voxel_sizes",
            is_optional=True,
        )
        self.add_link(
            "newsegment.bias_corrected_images->normalize12_1.apply_to_files"
        )
        self.export_parameter(
            "newsegment", "bias_corrected_images", is_optional=True
        )
        self.export_parameter(
            "newsegment", "bias_field_images", is_optional=True
        )
        self.export_parameter(
            "newsegment", "native_class_images", is_optional=True
        )
        self.add_link(
            "newsegment.forward_deformation_field->"
            "normalize12_1.deformation_file"
        )
        self.add_link(
            "newsegment.forward_deformation_field->"
            "normalize12_2.deformation_file"
        )
        self.export_parameter(
            "newsegment", "forward_deformation_field", is_optional=True
        )
        self.add_link("realign.realigned_files->slicetiming.in_files")
        self.add_link("realign.mean_image->coregister.target")
        self.export_parameter(
            "realign", "realignment_parameters", is_optional=True
        )
        self.export_parameter(
            "normalize12_1",
            "normalized_files",
            "normalized_anat",
            is_optional=True,
        )
        self.add_link("normalize12_2.normalized_files->smooth.in_files")
        self.export_parameter(
            "normalize12_2",
            "normalized_files",
            "normalized_func",
            is_optional=True,
        )
        self.export_parameter(
            "smooth", "smoothed_files", "smoothed_func", is_optional=False
        )
        self.add_link(
            "coregister.coregistered_source->newsegment.channel_files"
        )
        self.add_link("slicetiming.timed_files->normalize12_2.apply_to_files")

        # parameters order

        self.reorder_traits(
            (
                "anat_file",
                "func_files",
                "bias_corrected_images",
                "bias_field_images",
                "native_class_images",
                "forward_deformation_field",
                "realignment_parameters",
                "normalized_anat",
                "smoothed_func",
                "normalized_func",
                "st_acquisition",
                "func_out_voxel_sizes",
                "anat_out_voxel_sizes",
            )
        )

        # nodes positions
        self.node_position = {
            "outputs": (856.5917331726254, 147.0485132211699),
            "inputs": (-748.8507763768052, 501.7987014841776),
            "normalize12_1": (268.90900279449147, 180.5407227685563),
            "normalize12_2": (312.0451407533323, 809.2698435973882),
            "newsegment": (-374.0761052205569, 87.30628305829183),
            "realign": (-734.6135332276835, 781.3392612161631),
            "coregister": (-81.44182796634325, 732.0894759161474),
            "smooth": (784.074353728523, 809.5806776279725),
            "slicetiming": (-343.46222181427765, 937.4051799792994),
        }

        # nodes dimensions
        self.node_dimension = {
            "inputs": (187.8125, 215.0),
            "newsegment": (407.296875, 355.0),
            "realign": (330.375, 530.0),
            "normalize12_1": (351.640625, 565.0),
            "normalize12_2": (351.640625, 565.0),
            "smooth": (259.109375, 215.0),
            "coregister": (271.03125, 460.0),
            "outputs": (216.84375, 320.0),
            "slicetiming": (188.25, 320.0),
        }

        self.do_autoexport_nodes_parameters = False
