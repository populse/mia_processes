# -*- coding: utf-8 -*-
"""
The Bold_mni_align pipeline.

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

import os

import traits.api as traits
from capsul.api import Pipeline


class Bold_mni_align(Pipeline):
    """
    *Registered a mean functional image EPI in MNI space*

    For checking out registration.

    Please, see the complete documentation for the
    `Bold_mni_align pipeline in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Bold_mni_align.html>`_


    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "affineinitializer",
            "mia_processes.bricks.preprocess."
            "ants.processes.AffineInitializer",
        )
        self.nodes["affineinitializer"].process.num_threads = int(
            os.getenv("OMP_NUM_THREADS", os.cpu_count())
        )
        self.add_process(
            "registration",
            "mia_processes.bricks.preprocess." "ants.processes.Registration",
        )
        self.nodes["registration"].process.fixed_image_masks = traits.Undefined
        self.nodes["registration"].process.convergence_threshold = [
            1e-06,
            1e-06,
            1e-06,
        ]
        self.nodes["registration"].process.convergence_window_size = [
            20,
            20,
            10,
        ]
        self.nodes[
            "registration"
        ].process.interpolation = "LanczosWindowedSinc"
        self.nodes["registration"].process.metric = ["Mattes", "Mattes", "CC"]
        self.nodes["registration"].process.metric_weight = [1, 1, 1]
        self.nodes["registration"].process.number_of_iterations = [
            [10000, 1000, 100],
            [500, 250, 100],
            [100, 30, 20],
        ]
        self.nodes["registration"].process.radius_or_number_of_bins = [
            56,
            56,
            4,
        ]
        self.nodes["registration"].process.sampling_percentage = [
            0.25,
            0.25,
            1.0,
        ]
        self.nodes["registration"].process.sampling_strategy = [
            "Regular",
            "Regular",
            "None",
        ]
        self.nodes["registration"].process.shrink_factors = [
            [4, 2, 1],
            [8, 4, 2],
            [8, 4, 2],
        ]
        self.nodes["registration"].process.smoothing_sigmas = [
            [4.0, 2.0, 1.0],
            [4.0, 2.0, 0.0],
            [3.0, 2.0, 1.0],
        ]
        self.nodes["registration"].process.transform_parameters = [
            (0.05,),
            (0.08,),
            (0.1, 3.0, 0.0),
        ]
        self.nodes["registration"].process.transforms = [
            "Rigid",
            "Affine",
            "SyN",
        ]
        self.nodes["registration"].process.num_threads = int(
            os.getenv("OMP_NUM_THREADS", os.cpu_count())
        )
        self.add_process(
            "n4biasfieldcorrection",
            "mia_processes.bricks.preprocess."
            "ants.processes.N4BiasFieldCorrection",
        )
        self.nodes["n4biasfieldcorrection"].process.dimension = 3
        self.add_process(
            "applytransforms",
            "mia_processes.bricks.preprocess."
            "ants.processes.ApplyTransforms",
        )
        self.nodes["applytransforms"].process.interpolation = "MultiLabel"
        self.add_process(
            "template",
            "mia_processes.bricks.preprocess."
            "others.processing.TemplateFromTemplateFlow",
        )
        self.nodes["template"].process.suffix = "boldref"
        self.nodes["template"].process.desc = "fMRIPrep"
        self.add_process(
            "template_mask",
            "mia_processes.bricks.preprocess."
            "others.processing.TemplateFromTemplateFlow",
        )
        self.nodes["template_mask"].process.suffix = "mask"
        self.nodes["template_mask"].process.desc = "brain"
        self.add_process(
            "seg_template",
            "mia_processes.bricks.preprocess."
            "others.processing.TemplateFromTemplateFlow",
        )
        self.nodes["seg_template"].process.suffix = "dseg"
        self.nodes["seg_template"].process.desc = "carpet"

        # links
        self.export_parameter(
            "registration", "moving_image_masks", "epi_mask", is_optional=True
        )
        self.export_parameter(
            "n4biasfieldcorrection", "in_file", "epi_mean", is_optional=False
        )
        self.add_link("epi_mean->applytransforms.reference_image")
        self.export_parameter(
            "template", "in_template", "template", is_optional=True
        )
        self.add_link("template->template_mask.in_template")
        self.add_link("template->seg_template.in_template")
        self.export_parameter(
            "template_mask", "resolution", "template_res", is_optional=True
        )
        self.add_link("template_res->template.resolution")
        self.export_parameter(
            "seg_template", "resolution", "seg_template_res", is_optional=True
        )
        self.export_parameter(
            "registration", "transforms", "reg_transforms", is_optional=True
        )
        self.export_parameter(
            "registration",
            "transform_parameters",
            "reg_transform_parameters",
            is_optional=True,
        )
        self.export_parameter(
            "registration", "metric", "reg_metric", is_optional=True
        )
        self.export_parameter(
            "registration",
            "metric_weight",
            "reg_metric_weight",
            is_optional=True,
        )
        self.export_parameter(
            "registration",
            "shrink_factors",
            "reg_shrink_factors",
            is_optional=True,
        )
        self.export_parameter(
            "registration",
            "smoothing_sigmas",
            "reg_smoothing_sigmas",
            is_optional=True,
        )
        self.export_parameter(
            "registration",
            "number_of_iterations",
            "reg_number_of_iterations",
            is_optional=True,
        )
        self.export_parameter(
            "registration",
            "radius_or_number_of_bins",
            "reg_radius_or_number_of_bins",
            is_optional=True,
        )
        self.export_parameter(
            "registration",
            "convergence_threshold",
            "reg_convergence_threshold",
            is_optional=True,
        )
        self.export_parameter(
            "registration",
            "convergence_window_size",
            "reg_convergence_window_size",
            is_optional=True,
        )
        self.export_parameter(
            "registration",
            "sampling_percentage",
            "reg_sampling_percentage",
            is_optional=True,
        )
        self.export_parameter(
            "registration",
            "sampling_strategy",
            "reg_sampling_strategy",
            is_optional=True,
        )
        self.export_parameter(
            "registration",
            "interpolation",
            "reg_interpolation",
            is_optional=True,
        )
        self.add_link(
            "affineinitializer.out_file->"
            "registration.initial_moving_transform"
        )
        self.export_parameter(
            "registration", "composite_transform", is_optional=True
        )
        self.export_parameter(
            "registration", "inverse_composite_transform", is_optional=True
        )
        self.add_link(
            "registration.inverse_composite_transform->"
            "applytransforms.transforms"
        )
        self.export_parameter(
            "registration", "warped_image", "epi_mni", is_optional=False
        )
        self.add_link(
            "n4biasfieldcorrection.out_file->" "registration.moving_image"
        )
        self.add_link(
            "n4biasfieldcorrection.out_file->" "affineinitializer.moving_image"
        )
        self.export_parameter(
            "n4biasfieldcorrection", "bias_image", is_optional=True
        )
        self.export_parameter(
            "applytransforms", "output_image", "epi_parc", is_optional=False
        )
        self.add_link("template.template->registration.fixed_image")
        self.add_link("template.template->affineinitializer.fixed_image")
        self.add_link("template_mask.template->registration.fixed_image_masks")
        self.add_link("seg_template.template->applytransforms.input_image")

        # parameters order
        self.reorder_traits(
            (
                "epi_mask",
                "epi_mean",
                "epi_parc",
                "template",
                "template_res",
                "epi_mni",
                "composite_transform",
                "inverse_composite_transform",
                "seg_template_res",
                "bias_image",
                "reg_transforms",
                "reg_transform_parameters",
                "reg_metric",
                "reg_metric_weight",
                "reg_shrink_factors",
                "reg_smoothing_sigmas",
                "reg_number_of_iterations",
                "reg_radius_or_number_of_bins",
                "reg_convergence_threshold",
                "reg_convergence_window_size",
                "reg_sampling_percentage",
                "reg_sampling_strategy",
                "reg_interpolation",
            )
        )

        # default and initial values
        self.template = "MNI152NLin2009cAsym"
        self.template_res = 2
        self.seg_template_res = 1

        # nodes positions
        self.node_position = {
            "affineinitializer": (135.11715839999994, 47.53198080000004),
            "registration": (631.8374400000001, 182.46527999999978),
            "n4biasfieldcorrection": (-390.90908160000004, 151.86078719999995),
            "applytransforms": (1139.6390399999998, 1051.0175231999995),
            "inputs": (-638.3967713999999, 726.1401600000002),
            "outputs": (1283.0001249999996, 573.5807999999995),
            "template": (-295.4831104000001, 375.1539712000002),
            "template_mask": (-108.99864319999998, 505.4964480000002),
            "seg_template": (711.2617855999997, 1029.3064064),
        }

        # nodes dimensions
        self.node_dimension = {
            "affineinitializer": (184.234375, 145.0),
            "registration": (440.9375, 705.0),
            "n4biasfieldcorrection": (209.78125, 145.0),
            "applytransforms": (242.328125, 215.0),
            "inputs": (227.765625, 670.0),
            "outputs": (226.78125, 215.0),
            "template": (174.640625, 215.0),
            "template_mask": (174.640625, 215.0),
            "seg_template": (180.640625, 215.0),
        }

        self.do_autoexport_nodes_parameters = False
