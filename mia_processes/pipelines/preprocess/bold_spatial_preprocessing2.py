# -*- coding: utf-8 -*-

"""The pipelines library of the mia_processes package.

The objective of this module is to propose pipelines built from atomic bricks
or other pipelines provided in the mia_processes library.

:Contains:
    :Class:
        - Bold_spatial_preprocessing2

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################


from capsul.api import Pipeline


class Bold_spatial_preprocessing2(Pipeline):
    """
    *Bold data pre-processing*

    Please, see the complete documentation for the
    `Bold_spatial_preprocessing2 brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Bold_spatial_preprocessing2.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "normalize12_T1_estwrite",
            "mia_processes.bricks.preprocess.spm.spatial_preprocessing."
            "Normalize12",
        )
        self.nodes["normalize12_T1_estwrite"].process.jobtype = "estwrite"
        self.nodes[
            "normalize12_T1_estwrite"
        ].process.bias_regularization = 0.0001
        self.nodes["normalize12_T1_estwrite"].process.bias_fwhm = 30
        self.nodes["normalize12_T1_estwrite"].process.write_interp = 4

        self.add_process(
            "realign",
            "mia_processes.bricks.preprocess.spm.spatial_preprocessing."
            "Realign",
        )

        self.add_process(
            "coregister",
            "mia_processes.bricks.preprocess.spm.spatial_preprocessing."
            "Coregister",
        )

        self.add_process(
            "normalize12_func",
            "mia_processes.bricks.preprocess.spm.spatial_preprocessing."
            "Normalize12",
        )
        self.nodes["normalize12_func"].process.write_interp = 4

        self.add_process(
            "smooth",
            "mia_processes.bricks.preprocess.spm.spatial_preprocessing.Smooth",
        )

        self.add_process(
            "slicetiming",
            "mia_processes.bricks.preprocess.spm.spatial_preprocessing."
            "SliceTiming",
        )

        # links
        self.export_parameter(
            "slicetiming", "in_files", "func_files", is_optional=False
        )
        self.export_parameter(
            "coregister", "source", "anat_image", is_optional=False
        )
        self.export_parameter(
            "normalize12_func",
            "write_voxel_sizes",
            "write_voxel_sizes_irmf",
            is_optional=True,
        )
        self.export_parameter(
            "slicetiming", "acquisition", "acquisition_type", is_optional=True
        )
        self.add_link(
            "normalize12_T1_estwrite.deformation_field->normalize12_func."
            "deformation_file"
        )
        self.export_parameter(
            "normalize12_T1_estwrite",
            "normalized_files",
            "normalized_anat",
            is_optional=True,
        )
        self.add_link(
            "realign.realigned_files->normalize12_func.apply_to_files"
        )
        self.add_link("realign.mean_image->coregister.target")
        self.export_parameter(
            "realign", "realignment_parameters", is_optional=True
        )
        self.add_link(
            "coregister.coregistered_source->normalize12_T1_estwrite."
            "image_to_align"
        )
        self.add_link("normalize12_func.normalized_files->smooth.in_files")
        self.export_parameter(
            "smooth", "smoothed_files", "smoothed_func", is_optional=False
        )
        self.add_link("slicetiming.timed_files->realign.in_files")

        # parameters order
        self.reorder_traits(
            (
                "func_files",
                "anat_image",
                "write_voxel_sizes_irmf",
                "smoothed_func",
                "realignment_parameters",
                "normalized_anat",
                "acquisition_type",
            )
        )

        # nodes positions
        self.node_position = {
            "inputs": (-1330.6501789831684, -288.15260556753907),
            "normalize12_T1_estwrite": (
                -628.1150643192523,
                -38.853194422614536,
            ),
            "realign": (-719.1690583940288, -645.6232551039974),
            "coregister": (-1051.81769119588, 3.0859261395252133),
            "normalize12_func": (-95.96325602200594, -623.0809886962197),
            "smooth": (399.6456470428327, -243.53983169313443),
            "outputs": (829.5548202409861, 109.14575051340813),
            "slicetiming": (-975.8860749597629, -612.2454479258214),
        }

        # nodes dimensions
        self.node_dimension = {
            "slicetiming": (188.25, 320.0),
            "inputs": (199.609375, 140.0),
            "realign_1": (330.375, 530.0),
            "coregister_1": (271.03125, 460.0),
            "normalize12_func": (351.109375, 530.0),
            "normalize12_T1_estwrite": (351.109375, 530.0),
            "realign": (330.375, 530.0),
            "coregister": (271.03125, 460.0),
            "smooth": (259.109375, 215.0),
            "outputs": (196.28125, 145.0),
        }

        self.do_autoexport_nodes_parameters = False
