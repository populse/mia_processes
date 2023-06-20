# -*- coding: utf-8 -*-
"""
The Anat_mni_tpms pipeline.

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from capsul.api import Pipeline


class Anat_mni_tpms(Pipeline):
    """
    *Get templates tissues probability map (white matter - WM -,
    grey matter - GM - and cerebrospinal - CSF -)*

    From 'MNI152NLin2009cAsym' template
    (resolution 1 mm) and register them in subject space.

    Adapted from `mriqc <https://mriqc.readthedocs.io/en/22.0.6//>`_
    anatomical workflow.

    Please, see the complete documentation for the
    `Anat_mni_tpms pipeline in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Anat_mni_tpms.html>`_
    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "template_CSF",
            "mia_processes.bricks.preprocess."
            "others.processing.TemplateFromTemplateFlow",
        )
        self.nodes["template_CSF"].process.in_template = "MNI152NLin2009cAsym"
        self.nodes["template_CSF"].process.resolution = 1
        self.nodes["template_CSF"].process.suffix = "probseg"
        self.nodes["template_CSF"].process.label = "CSF"
        self.add_process(
            "applytransforms_CSF",
            "mia_processes.bricks.preprocess."
            "ants.processes.ApplyTransforms",
        )
        self.nodes["applytransforms_CSF"].process.out_prefix = "csf_"
        self.nodes["applytransforms_CSF"].process.float = True
        self.add_process(
            "template_GM",
            "mia_processes.bricks.preprocess"
            ".others.processing.TemplateFromTemplateFlow",
        )
        self.nodes["template_GM"].process.in_template = "MNI152NLin2009cAsym"
        self.nodes["template_GM"].process.resolution = 1
        self.nodes["template_GM"].process.suffix = "probseg"
        self.nodes["template_GM"].process.label = "GM"
        self.add_process(
            "applytransforms_GM",
            "mia_processes.bricks.preprocess."
            "ants.processes.ApplyTransforms",
        )
        self.nodes["applytransforms_GM"].process.out_prefix = "gm_"
        self.nodes["applytransforms_GM"].process.float = True
        self.add_process(
            "template_WM",
            "mia_processes.bricks.preprocess."
            "others.processing.TemplateFromTemplateFlow",
        )
        self.nodes["template_WM"].process.in_template = "MNI152NLin2009cAsym"
        self.nodes["template_WM"].process.resolution = 1
        self.nodes["template_WM"].process.suffix = "probseg"
        self.nodes["template_WM"].process.label = "WM"
        self.add_process(
            "applytransforms_WM",
            "mia_processes.bricks.preprocess."
            "ants.processes.ApplyTransforms",
        )
        self.nodes["applytransforms_WM"].process.out_prefix = "wm_"
        self.nodes["applytransforms_WM"].process.float = True
        self.add_process(
            "files_to_list", "mia_processes.bricks.tools.tools.Files_To_List"
        )

        # links
        self.export_parameter(
            "applytransforms_CSF",
            "reference_image",
            "in_ras",
            is_optional=False,
        )
        self.add_link("in_ras->applytransforms_GM.reference_image")
        self.add_link("in_ras->applytransforms_WM.reference_image")
        self.export_parameter(
            "applytransforms_WM",
            "transforms",
            "inverse_composite_transform",
            is_optional=False,
        )
        self.add_link(
            "inverse_composite_transform->" "applytransforms_CSF.transforms"
        )
        self.add_link(
            "inverse_composite_transform->" "applytransforms_GM.transforms"
        )
        self.add_link("template_CSF.template->applytransforms_CSF.input_image")
        self.add_link("applytransforms_CSF.output_image->files_to_list.file1")
        self.add_link("template_GM.template->applytransforms_GM.input_image")
        self.add_link("template_WM.template->applytransforms_WM.input_image")
        self.add_link("applytransforms_WM.output_image->files_to_list.file3")
        self.add_link("applytransforms_GM.output_image->files_to_list.file2")
        self.export_parameter(
            "files_to_list", "file_list", "mni_tpms", is_optional=False
        )

        # parameters order
        self.reorder_traits(
            ("in_ras", "inverse_composite_transform", "mni_tpms")
        )

        # default and initial values
        self.template = "MNI152NLin2009cAsym"
        self.template_res = 1
        self.template_suffix = "probseg"

        # nodes positions
        self.node_position = {
            "inputs": (-1029.0321499999993, -164.47839999999997),
            "template_CSF": (-658.9776000000002, -507.5264),
            "applytransforms_CSF": (-384.7936, -512.6864),
            "template_GM": (-657.7967999999998, -230.95039999999992),
            "template_WM": (-669.3920000000002, 52.54400000000007),
            "applytransforms_WM": (-382.3792000000001, 47.635200000000054),
            "applytransforms_GM": (-385.4448, -223.1727999999999),
            "files_to_list": (-35.82719999999998, -185.29919999999998),
            "outputs": (169.08355000000003, -150.73919999999998),
        }

        # nodes dimensions
        self.node_dimension = {
            "inputs": (244.08125, 215.0),
            "template_CSF": (177.640625, 250.0),
            "applytransforms_CSF": (242.328125, 215.0),
            "template_GM": (177.640625, 250.0),
            "template_WM": (177.640625, 250.0),
            "applytransforms_WM": (242.328125, 215.0),
            "applytransforms_GM": (242.328125, 215.0),
            "files_to_list": (119.125, 145.0),
            "outputs": (91.09375, 75.0),
        }

        self.do_autoexport_nodes_parameters = False
