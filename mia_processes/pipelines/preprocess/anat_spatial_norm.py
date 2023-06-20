# -*- coding: utf-8 -*-
"""
The Anat_spatial_norm pipeline.

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


class Anat_spatial_norm(Pipeline):
    """
    *Spatial normalization to MNI (using 'MNI152NLin2009cAsym' template)*

    Adapted from `mriqc 22.06 anatomical workflow
    <https://github.com/nipreps/mriqc/blob/5a0f0408bd0c176dbc46088c6ffe279269180f3f/mriqc/workflows/anatomical.py#L282>`_

    Please, see the complete documentation for the
    `Anat_spatial_norm pipeline in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Anat_spatial_norm.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "mask_moving_image",
            "mia_processes.bricks.preprocess." "others.processing.Mask",
        )
        self.add_process(
            "mask_fixed_image",
            "mia_processes.bricks.preprocess." "others.processing.Mask",
        )
        self.add_process(
            "affine_initializer",
            "mia_processes.bricks.preprocess."
            "ants.processes.AffineInitializer",
        )
        self.add_process(
            "registration",
            "mia_processes.bricks.preprocess." "ants.processes.Registration",
        )
        self.nodes["registration"].process.convergence_threshold = [
            1e-07,
            1e-08,
        ]
        self.nodes["registration"].process.convergence_window_size = [15, 5, 3]
        self.nodes[
            "registration"
        ].process.interpolation = "LanczosWindowedSinc"
        self.nodes["registration"].process.metric = ["Mattes", "Mattes"]
        self.nodes["registration"].process.metric_weight = [1.0, 1.0]
        self.nodes["registration"].process.number_of_iterations = [[20], [15]]
        self.nodes["registration"].process.radius_or_number_of_bins = [56, 56]
        self.nodes["registration"].process.sampling_percentage = [0.2, 0.1]
        self.nodes["registration"].process.sampling_strategy = [
            "Random",
            "Random",
        ]
        self.nodes["registration"].process.shrink_factors = [[2], [1]]
        self.nodes["registration"].process.smoothing_sigmas = [[4.0], [2.0]]
        self.nodes["registration"].process.transform_parameters = [
            (1.0,),
            (1.0,),
        ]
        self.nodes["registration"].process.transforms = ["Rigid", "Affine"]
        self.nodes["registration"].process.use_histogram_matching = [
            False,
            True,
        ]
        self.nodes["registration"].process.num_threads = int(
            os.getenv("OMP_NUM_THREADS", os.cpu_count())
        )
        self.nodes["affine_initializer"].process.num_threads = int(
            os.getenv("OMP_NUM_THREADS", os.cpu_count())
        )
        self.add_process(
            "template_mask",
            "mia_processes.bricks.preprocess."
            "others.processing.TemplateFromTemplateFlow",
        )
        self.nodes["template_mask"].process.suffix = "mask"
        self.nodes["template_mask"].process.desc = "brain"
        self.nodes["template_mask"].process.in_template = "MNI152NLin2009cAsym"
        self.nodes["template_mask"].process.resolution = 2
        self.add_process(
            "template",
            "mia_processes.bricks.preprocess."
            "others.processing.TemplateFromTemplateFlow",
        )
        self.nodes["template"].process.suffix = "T1w"
        self.nodes["template"].process.in_template = "MNI152NLin2009cAsym"
        self.nodes["template"].process.resolution = 2

        # links
        self.export_parameter(
            "mask_moving_image", "in_file", "moving_image", is_optional=False
        )
        self.export_parameter(
            "mask_moving_image", "mask_file", "moving_mask", is_optional=False
        )
        self.add_link(
            "mask_moving_image.out_file->" "affine_initializer.moving_image"
        )
        self.add_link("mask_moving_image.out_file->registration.moving_image")
        self.add_link(
            "mask_fixed_image.out_file->" "affine_initializer.fixed_image"
        )
        self.add_link("mask_fixed_image.out_file->registration.fixed_image")
        self.add_link(
            "affine_initializer.out_file->"
            "registration.initial_moving_transform"
        )
        self.export_parameter(
            "registration", "composite_transform", is_optional=True
        )
        self.export_parameter(
            "registration", "inverse_composite_transform", is_optional=False
        )
        self.export_parameter(
            "registration", "warped_image", is_optional=False
        )
        self.add_link("template_mask.template->mask_fixed_image.mask_file")
        self.add_link("template.template->mask_fixed_image.in_file")

        # parameters order
        self.reorder_traits(
            (
                "moving_image",
                "moving_mask",
                "composite_transform",
                "inverse_composite_transform",
                "warped_image",
            )
        )

        # default and initial values
        self.template = "MNI152NLin2009cAsym"
        self.template_res = 2

        # nodes positions
        self.node_position = {
            "mask_moving_image": (-336.81600000000003, -20.0),
            "mask_fixed_image": (-347.05600000000004, 229.51999999999995),
            "affine_initializer": (-63.07200000000006, 208.21600000000007),
            "registration": (231.85599999999988, 130.06400000000008),
            "template_mask": (-352.6961689849344, 451.25877932328956),
            "template": (-609.4425831662135, 454.5343161027413),
            "inputs": (-960.0347138972587, 611.731382193135),
            "outputs": (880.2747499999998, 130.06400000000008),
        }

        # nodes dimensions
        self.node_dimension = {
            "mask_moving_image": (180.3125, 180.0),
            "mask_fixed_image": (158.53125, 180.0),
            "affine_initializer": (184.234375, 145.0),
            "registration": (446.9375, 705.0),
            "template_mask": (174.640625, 215.0),
            "template": (180.640625, 215.0),
            "inputs": (227.765625, 635.0),
            "outputs": (226.78125, 145.0),
        }

        self.do_autoexport_nodes_parameters = False
