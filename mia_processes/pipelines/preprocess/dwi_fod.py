# -*- coding: utf-8 -*-

"""The pipelines library of the mia_processes package.

The objective of this module is to propose pipelines built from atomic bricks
or other pipelines provided in the mia_processes library.

:Contains:
    :Class:
        - Dwi_fod

"""

from capsul.api import Pipeline


class Dwi_fod(Pipeline):

    """
    *DWI fiber orientation distribution*

    Please, see the complete documentation for the
    `Dwi_fod pipeline in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Dwi_fod.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "ResponseFunctionEstimation",
            "mia_processes.bricks.preprocess."
            "mrtrix.processes.ResponseSDDhollander",
        )
        self.add_process(
            "FODEstimation",
            "mia_processes.bricks.preprocess."
            "mrtrix.processes.ConstrainedSphericalDeconvolution",
        )
        self.nodes["FODEstimation"].process.algorithm = "msmt_csd"
        self.nodes["FODEstimation"].process.get_predicted_signal = True
        self.add_process(
            "IntensityNormalization",
            "mia_processes.bricks.preprocess.mrtrix.processes.MTnormalise",
        )
        self.add_process(
            "FilesToList", "mia_processes.bricks.tools.tools.Files_To_List"
        )
        self.add_process(
            "ListToFile_1", "mia_processes.bricks.tools.tools.List_To_File"
        )
        self.nodes["ListToFile_1"].process.index_filter = [1]
        self.add_process(
            "ListToFile_2", "mia_processes.bricks.tools.tools.List_To_File"
        )
        self.nodes["ListToFile_2"].process.index_filter = [2]
        self.add_process(
            "ListToFile_3", "mia_processes.bricks.tools.tools.List_To_File"
        )
        self.nodes["ListToFile_3"].process.index_filter = [3]

        # links
        self.export_parameter(
            "ResponseFunctionEstimation",
            "in_file",
            "in_dwi",
            is_optional=False,
        )
        self.add_link("in_dwi->FODEstimation.in_file")
        self.export_parameter(
            "FODEstimation", "in_mask", "brain_mask", is_optional=False
        )
        self.add_link("brain_mask->IntensityNormalization.mask")
        self.add_link(
            "ResponseFunctionEstimation.csf_file->FODEstimation.csf_txt"
        )
        self.add_link(
            "ResponseFunctionEstimation.gm_file->FODEstimation.gm_txt"
        )
        self.add_link(
            "ResponseFunctionEstimation.wm_file->FODEstimation.wm_txt"
        )
        self.export_parameter(
            "ResponseFunctionEstimation",
            "voxels_image",
            "voxels_fianl_selection",
            is_optional=True,
        )
        self.add_link("FODEstimation.csf_odf->FilesToList.file1")
        self.add_link("FODEstimation.gm_odf->FilesToList.file2")
        self.add_link("FODEstimation.wm_odf->FilesToList.file3")
        self.export_parameter(
            "FODEstimation", "predicted_signal_file", is_optional=True
        )
        self.add_link(
            "IntensityNormalization.out_files->ListToFile_1.file_list"
        )
        self.add_link(
            "IntensityNormalization.out_files->ListToFile_2.file_list"
        )
        self.add_link(
            "IntensityNormalization.out_files->ListToFile_3.file_list"
        )
        self.add_link("FilesToList.file_list->IntensityNormalization.in_files")
        self.export_parameter(
            "ListToFile_1", "file", "csf_fod_norm", is_optional=False
        )
        self.export_parameter(
            "ListToFile_2", "file", "gm_fod_norm", is_optional=False
        )
        self.export_parameter(
            "ListToFile_3", "file", "wm_fod_norm", is_optional=False
        )

        # parameters order

        self.reorder_traits(
            (
                "in_dwi",
                "voxels_fianl_selection",
                "brain_mask",
                "predicted_signal_file",
                "csf_fod_norm",
                "gm_fod_norm",
                "wm_fod_norm",
            )
        )

        # nodes positions
        self.node_position = {
            "ResponseFunctionEstimation": (
                -1262.912368640373,
                152.4026270677166,
            ),
            "FODEstimation": (-933.6347733561715, -153.54350458013056),
            "IntensityNormalization": (
                -385.88257304070976,
                -94.08005523992625,
            ),
            "FilesToList": (-561.7630640303614, -71.30995697733732),
            "inputs": (-1358.4880421537448, -98.81792065248086),
            "outputs": (198.02645870077032, 57.60027545906087),
            "ListToFile_1": (-79.00017223905586, -104.28022735555368),
            "ListToFile_2": (-82.16017912861807, 47.40010334343346),
            "ListToFile_3": (-82.16017912861801, 210.14045815588844),
        }

        # nodes dimensions
        self.node_dimension = {
            "dwidenoise_1": (156.234375, 180.0),
            "ResponseFunctionEstimation": (252.125, 285.0),
            "FODEstimation": (334.625, 390.0),
            "IntensityNormalization": (229.1875, 250.0),
            "FilesToList": (118.8125, 145.0),
            "inputs": (106.984375, 110.0),
            "outputs": (178.296875, 215.0),
            "ListToFile_1": (139.703125, 110.0),
            "ListToFile_2": (139.703125, 110.0),
            "ListToFile_3": (139.703125, 110.0),
        }

        self.do_autoexport_nodes_parameters = False
