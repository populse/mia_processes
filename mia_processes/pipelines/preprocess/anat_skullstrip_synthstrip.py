# -*- coding: utf-8 -*-
"""
The .Anat_skullstrip_synthstrip pipeline.

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


class Anat_skullstrip_synthstrip(Pipeline):
    """
    *Create a brain-extraction workflow using SynthStrip (Freesurfer)*

    Adapted from `mriqc v22.06 synthstrip workflow
    <https://github.com/nipreps/mriqc/blob/5a0f0408bd0c176dbc46088c6ffe279269180f3f/mriqc/workflows/anatomical.py#L849>`_

    Please, see the complete documentation for the
    `Anat_skullstrip_synthstrip pipeline in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Anat_skullstrip_synthstrip.html>`_


    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "pre_n4biasfieldcor",
            "mia_processes.bricks.preprocess."
            "ants.processes.N4BiasFieldCorrection",
        )
        self.add_process(
            "pre_clip",
            "mia_processes.bricks.preprocess."
            "others.processing.IntensityClip",
        )
        self.add_process(
            "synthstrip",
            "mia_processes.bricks.preprocess."
            "freesurfer.processes.SynthStripMriqc",
        )
        self.add_process(
            "post_n4biasfieldcor",
            "mia_processes.bricks.preprocess."
            "ants.processes.N4BiasFieldCorrection",
        )
        self.add_process(
            "mask", "mia_processes.bricks.preprocess." "others.processing.Mask"
        )

        self.nodes["pre_n4biasfieldcor"].process.dimension = 3
        self.nodes["pre_n4biasfieldcor"].process.out_prefix = "pre_n4c_"
        self.nodes["pre_n4biasfieldcor"].process.rescale_intensities = True
        self.nodes["pre_n4biasfieldcor"].process.num_threads = int(
            os.getenv("OMP_NUM_THREADS", os.cpu_count())
        )
        self.nodes["post_n4biasfieldcor"].process.n_iterations = [50] * 4
        self.nodes["post_n4biasfieldcor"].process.out_prefix = "post_n4c_"
        self.nodes["post_n4biasfieldcor"].process.dimension = 3
        self.nodes["post_n4biasfieldcor"].process.num_threads = int(
            os.getenv("OMP_NUM_THREADS", os.cpu_count())
        )
        self.nodes["mask"].process.suffix = ""
        self.nodes["mask"].process.prefix = "ss_"

        # links
        self.export_parameter("pre_clip", "in_file", is_optional=False)
        self.add_link("pre_n4biasfieldcor.out_file->synthstrip.in_file")
        self.add_link("pre_clip.out_file->pre_n4biasfieldcor.in_file")
        self.add_link("pre_clip.out_file->post_n4biasfieldcor.in_file")
        self.add_link(
            "synthstrip.out_mask->" "post_n4biasfieldcor.weight_image"
        )
        self.add_link("synthstrip.out_mask->mask.mask_file")
        self.add_link("post_n4biasfieldcor.out_file->mask.in_file")
        self.export_parameter(
            "synthstrip", "out_mask", "out_mask_synthstrip", is_optional=True
        )
        self.export_parameter(
            "post_n4biasfieldcor", "bias_image", is_optional=True
        )
        self.export_parameter(
            "post_n4biasfieldcor",
            "out_file",
            "out_corrected",
            is_optional=True,
        )
        self.export_parameter(
            "mask", "out_file", "out_brain", is_optional=True
        )

        self.reorder_traits(
            (
                "in_file",
                "out_mask_synthstrip",
                "out_brain",
                "out_corrected",
                "bias_image",
            )
        )

        # nodes positions
        self.node_position = {
            "inputs": (-1151.61875, -639.0),
            "pre_n4biasfieldcor": (-797.0, -502.0),
            "pre_clip": (-1022.0, -417.0),
            "synthstrip": (-561.0, -692.0),
            "post_n4biasfieldcor": (-406.0, -434.0),
            "mask": (-62.0, -528.0),
            "outputs": (-31.97187500000001, -749.0),
        }

        # nodes dimensions
        self.node_dimension = {
            "intensityclip_1": (168.03125, 215.0),
            "n4biasfieldcorrection_1": (236.484375, 180.0),
            "n4biasfieldcorrection_2": (236.484375, 180.0),
            "synthstrip_1": (176.53125, 110.0),
            "mask_1": (149.453125, 180.0),
            "inputs": (99.328125, 110.0),
            "pre_n4biasfieldcor": (236.484375, 180.0),
            "pre_clip": (168.03125, 215.0),
            "synthstrip": (179.53125, 110.0),
            "post_n4biasfieldcor": (239.484375, 180.0),
            "mask": (149.453125, 180.0),
            "outputs": (165.421875, 215.0),
        }

        self.do_autoexport_nodes_parameters = False
