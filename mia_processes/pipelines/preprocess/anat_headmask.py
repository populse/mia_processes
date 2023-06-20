# -*- coding: utf-8 -*-
"""
The Anat_headmask pipeline..

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from capsul.api import Pipeline


class Anat_headmask(Pipeline):
    """
    *Compute head mask from skull stripped structural image using "Denoise"
    from Dipy*

    Please, see the complete documentation for the
    `Anat_headmask pipeline in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Anat_headmask.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "estimateSNR",
            "mia_processes.bricks.preprocess." "others.processing.EstimateSNR",
        )
        self.add_process(
            "denoise",
            "mia_processes.bricks.preprocess." "dipy.processes.Denoise",
        )
        self.add_process(
            "enhance",
            "mia_processes.bricks.preprocess." "others.processing.Enhance",
        )
        self.add_process(
            "gradient_threshold",
            "mia_processes.bricks.preprocess."
            "others.processing.GradientThreshold",
        )

        # links
        self.export_parameter(
            "estimateSNR", "in_file", "in_file", is_optional=False
        )
        self.export_parameter(
            "gradient_threshold", "seg_file", is_optional=False
        )
        self.add_link("in_file->enhance.in_files")
        self.add_link("seg_file->estimateSNR.seg_file")
        self.add_link("estimateSNR.out_snr->denoise.snr")
        self.add_link("enhance.out_files->denoise.in_file")
        self.add_link("denoise.out_file->gradient_threshold.in_file")
        self.export_parameter(
            "gradient_threshold", "out_file", is_optional=False
        )

        # parameters order

        self.reorder_traits(("in_file", "seg_file", "out_file"))

        # nodes positions
        self.node_position = {
            "estimateSNR": (-100.0, -70.0),
            "denoise": (1.0, -70.0),
            "inputs": (-346.415625, 4.0),
            "enhance": (-180.0, 182.0),
            "gradient_threshold": (216.0, 125.0),
            "outputs": (487.028125, 125.0),
        }

        # nodes dimensions
        self.node_dimension = {
            "estimateSNR": (178.265625, 60.0),
            "denoise": (157.84375, 180.0),
            "inputs": (86.503125, 110.0),
            "enhance": (141.546875, 145.0),
            "gradient_threshold": (178.265625, 180.0),
            "outputs": (79.0625, 60.0),
        }

        self.do_autoexport_nodes_parameters = False
