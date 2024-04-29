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

from capsul.api import Pipeline


class Bold_spatial_preprocessing3(Pipeline):
    """
    *Data pre-processing for cerebrovascular reserve analysis (CVR)*

    Please, see the complete documentation for the
    `Bold_spatial_preprocessing1 brick in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Bold_spatial_preprocessing1.html>`_
    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "realign",
            "mia_processes.bricks.preprocess."
            "spm.spatial_preprocessing.Realign",
        )
        self.add_process(
            "normalize12",
            "mia_processes.bricks.preprocess."
            "spm.spatial_preprocessing.Normalize12",
        )
        self.add_process(
            "smooth",
            "mia_processes.bricks.preprocess."
            "spm.spatial_preprocessing.Smooth",
        )
        self.add_process(
            "coregister",
            "mia_processes.bricks.preprocess."
            "spm.spatial_preprocessing.Coregister",
        )

        # links
        self.export_parameter(
            "coregister", "target", "anat_file", is_optional=False
        )
        self.export_parameter(
            "realign", "in_files", "func_files", is_optional=False
        )
        self.export_parameter(
            "normalize12",
            "deformation_file",
            "deformation_file_T1_MNI",
            is_optional=False,
        )
        self.add_link("realign.realigned_files->coregister.apply_to_files")
        self.add_link("realign.mean_image->coregister.source")
        self.export_parameter(
            "realign", "realignment_parameters", is_optional=True
        )
        self.add_link("normalize12.normalized_files->smooth.in_files")
        self.export_parameter(
            "normalize12",
            "normalized_files",
            "normalized_func",
            is_optional=True,
        )
        self.export_parameter(
            "smooth", "smoothed_files", "smoothed_func", is_optional=False
        )
        self.export_parameter(
            "coregister", "coregistered_source", is_optional=False
        )
        self.add_link(
            "coregister.coregistered_files->normalize12.apply_to_files"
        )

        # parameters order

        self.reorder_traits(
            (
                "anat_file",
                "func_files",
                "realignment_parameters",
                "smoothed_func",
                "coregistered_source",
                "normalized_func",
                "deformation_file_T1_MNI",
            )
        )

        # nodes positions
        self.node_position = {
            "outputs": (565.1211661051032, 302.32873456793226),
            "inputs": (-676.576610318852, 499.14176836383615),
            "normalize12": (209.65051626919063, 683.6824024812697),
            "realign": (-558.3410913273092, 658.4832223737252),
            "coregister": (-150.20722756102208, 272.981908700011),
            "smooth": (604.8341436366131, 704.9631128063114),
        }

        # nodes dimensions
        self.node_dimension = {
            "inputs": (203.3125, 145.0),
            "realign": (330.375, 530.0),
            "normalize12": (354.640625, 565.0),
            "smooth": (259.109375, 215.0),
            "coregister": (271.03125, 460.0),
            "outputs": (196.28125, 180.0),
        }

        self.do_autoexport_nodes_parameters = False
