# -*- coding: utf-8 -*-
"""
The Bold_hmc pipeline.

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from capsul.api import Pipeline


class Bold_hmc(Pipeline):
    """
    *Head-motion correction based on 3dvolreg from AFNI*

    Please, see the complete documentation for the
    `Bold_hmc pipeline in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Bold_hmc.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "despike",
            "mia_processes.bricks.preprocess." "afni.processes.Despike",
        )
        self.add_process(
            "deoblique",
            "mia_processes.bricks.preprocess." "afni.processes.RefitDeoblique",
        )
        self.add_process(
            "volreg",
            "mia_processes.bricks.preprocess." "afni.processes.Volreg",
        )
        self.nodes["volreg"].process.twopass = True
        self.nodes["despike"].process.despike = False
        self.nodes["deoblique"].process.deoblique = False

        # links
        self.export_parameter("despike", "in_file", is_optional=False)
        self.export_parameter("despike", "despike", is_optional=True)
        self.export_parameter("deoblique", "deoblique", is_optional=True)
        self.add_link("despike.out_file->deoblique.in_file")
        self.add_link("deoblique.out_file->volreg.in_file")
        self.export_parameter("volreg", "out_file", is_optional=False)
        self.export_parameter("volreg", "oned_file", is_optional=False)

        # parameters order
        self.reorder_traits(
            ("in_file", "despike", "deoblique", "out_file", "oned_file")
        )

        # nodes positions
        self.node_position = {
            "despike": (21.799999999999983, 294.59999999999997),
            "deoblique": (186.20000000000005, 439.9999999999999),
            "volreg": (353.6, 514.3999999999999),
            "inputs": (-575.73325, 294.9679999999999),
            "outputs": (633.9562500000001, 514.3999999999999),
        }

        # nodes dimensions
        self.node_dimension = {
            "despike": (165.65625, 180.0),
            "deoblique": (151.03125, 110.0),
            "volreg": (176.625, 250.0),
            "inputs": (110.26875, 152.0),
            "outputs": (90.03125, 83.0),
        }

        self.do_autoexport_nodes_parameters = False
