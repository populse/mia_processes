"""The pipelines library of the mia_processes package.

The Perfdsc pipeline is designed for cerebral perfusion measurement
using Dynamic Susceptibility Contrast MRI.

:Contains:
    :Class:
        - Perfdsc

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################
from capsul.api import Pipeline


class Perfdsc(Pipeline):
    """
    *Cerebral perfusion measurement by DSC MRI*

    Please, see the complete documentation for the
    `Perfdsc pipeline in the mia_processes website:
    <https://populse.github.io/mia_processes/html/documentation/pipelines/Perfusion/Perfdsc.html>`_
    """

    def pipeline_definition(self):
        """Building the pipeline"""

        # nodes
        self.add_process(
            "spatial_preprocessing",
            "mia_processes.pipelines.preprocess.bold_spatial_preprocessing1."
            "Bold_spatial_preprocessing1",
        )
        self.nodes["spatial_preprocessing"].process.nodes_activation = {
            "newsegment": True,
            "realign": True,
            "list_duplicate": True,
            "normalize12_1": True,
            "normalize12_2": True,
            "smooth": True,
            "coregister": True,
        }
        # We use exactly the same parameters as in Amigo.
        self.nodes["spatial_preprocessing"].process.nodes[
            "realign"
        ].process.register_to_mean = False
        # We can do the same with:
        # self.nodes["spatial_preprocessing"].process.nodes[
        # "realign"].set_plug_value("register_to_mean", False)
        self.nodes["spatial_preprocessing"].process.nodes[
            "normalize12_1"
        ].process.write_voxel_sizes = [1.0, 1.0, 1.0]
        self.nodes["spatial_preprocessing"].process.nodes[
            "normalize12_1"
        ].process.write_bounding_box = [
            [-78.0, -112.0, -50.0],
            [78.0, 76.0, 85.0],
        ]
        self.nodes["spatial_preprocessing"].process.nodes[
            "normalize12_1"
        ].process.write_interp = 1
        self.nodes["spatial_preprocessing"].process.nodes[
            "normalize12_2"
        ].process.write_voxel_sizes = [2.0, 2.0, 2.0]
        self.nodes["spatial_preprocessing"].process.nodes[
            "normalize12_2"
        ].process.write_bounding_box = [
            [-78.0, -112.0, -50.0],
            [78.0, 76.0, 85.0],
        ]
        self.nodes["spatial_preprocessing"].process.nodes[
            "normalize12_2"
        ].process.write_interp = 1
        self.add_process(
            "spatial_mask",
            "mia_processes.pipelines.preprocess.spatial_mask.Spatial_mask",
        )
        self.nodes["spatial_mask"].process.nodes_activation = {
            "gm_wm_normalize": True,
            "threshold_1": True,
            "smooth": True,
            "threshold_2": True,
            "resample1": True,
        }
        self.add_process(
            "make_aif", "mia_processes.bricks.tools.tools.Make_AIF"
        )
        self.add_process(
            "deconv_from_aif",
            "mia_processes.bricks.tools.tools.Deconv_from_aif",
        )
        self.add_process(
            "reportperfdsc",
            "mia_processes.bricks.reports.reporting.ReportPerfDsc",
        )
        self.nodes["reportperfdsc"].process.norm_anat_slices_gap = 5
        self.nodes["reportperfdsc"].process.norm_func_inf_slice_start = 13
        self.nodes["reportperfdsc"].process.norm_func_slices_gap = 2
        self.nodes["reportperfdsc"].process.CBV_vmin = 0.0
        self.nodes["reportperfdsc"].process.CBV_vmax = 8.0
        self.nodes["reportperfdsc"].process.CBF_vmin = 0.0
        self.nodes["reportperfdsc"].process.CBF_vmax = 120.0
        self.nodes["reportperfdsc"].process.Tmax_vmin = 0.0
        self.nodes["reportperfdsc"].process.Tmax_vmax = 6.0
        self.nodes["reportperfdsc"].process.MTT_vmin = 0.0
        self.nodes["reportperfdsc"].process.MTT_vmax = 8.0
        # links
        self.export_parameter(
            "spatial_preprocessing", "anat_file", is_optional=False
        )
        self.export_parameter(
            "spatial_preprocessing", "func_files", is_optional=False
        )
        self.add_link(
            "spatial_preprocessing.native_class_images->"
            "spatial_mask.native_class_images"
        )
        self.add_link(
            "spatial_preprocessing.forward_deformation_field->"
            "spatial_mask.deformation_file"
        )
        self.add_link(
            "spatial_preprocessing.realignment_parameters->"
            "reportperfdsc.realignment_parameters"
        )
        self.add_link(
            "spatial_preprocessing.normalized_anat->" "reportperfdsc.norm_anat"
        )
        self.nodes["spatial_preprocessing"].process.trait(
            "coregistered_source"
        ).userlevel = 1
        self.export_parameter(
            "spatial_preprocessing", "coregistered_source", is_optional=False
        )
        self.add_link(
            "spatial_preprocessing.smoothed_func->"
            "spatial_mask.smoothed_func"
        )
        self.add_link(
            "spatial_preprocessing.smoothed_func->make_aif.func_file"
        )
        self.add_link(
            "spatial_preprocessing.normalized_func->" "reportperfdsc.norm_func"
        )
        self.add_link(
            "spatial_preprocessing.smoothed_func->deconv_from_aif.func_file"
        )
        self.add_link("spatial_mask.mask_003->deconv_from_aif.mask_file")
        self.add_link("make_aif.aif_file->reportperfdsc.aif_file")
        self.add_link("make_aif.aif_file->deconv_from_aif.aif_file")
        self.add_link("deconv_from_aif.CBV_image->reportperfdsc.CBV_image")
        self.add_link("deconv_from_aif.CBF_image->reportperfdsc.CBF_image")
        self.add_link("deconv_from_aif.MTT_image->reportperfdsc.MTT_image")
        self.add_link("deconv_from_aif.Tmax_image->reportperfdsc.Tmax_image")
        self.nodes["deconv_from_aif"].process.trait("TTP_image").userlevel = 1
        self.export_parameter(
            "deconv_from_aif", "TTP_image", is_optional=False
        )
        self.nodes["deconv_from_aif"].process.trait("T0_image").userlevel = 1
        self.export_parameter("deconv_from_aif", "T0_image", is_optional=False)
        self.export_parameter("reportperfdsc", "report", is_optional=True)
        self.export_parameter(
            "reportperfdsc", "patient_info", is_optional=True
        )

        # parameters order
        self.reorder_traits(
            (
                "anat_file",
                "func_files",
                "coregistered_source",
                "TTP_image",
                "T0_image",
                "patient_info",
                "report",
            )
        )

        # nodes positions
        self.node_position = {
            "spatial_preprocessing": (-247.52, 68.59),
            "spatial_mask": (220.61, -129.10),
            "make_aif": (73.01, 379.44),
            "deconv_from_aif": (552.04, 318.076),
            "reportperfdsc": (925.85, -308.25),
            "inputs": (-473.60, 181.38),
            "outputs": (1245.23, 248.28),
        }

        # nodes dimensions
        self.node_dimension = {
            "spatial_preprocessing": (236.5, 355.0),
            "spatial_mask": (200.2, 145.0),
            "make_aif": (172.2, 250.0),
            "deconv_from_aif": (215.3, 320.0),
            "reportperfdsc": (212.5, 1265.0),
            "inputs": (96.8, 110.0),
            "outputs": (59.6, 75.0),
        }

        self.do_autoexport_nodes_parameters = False
