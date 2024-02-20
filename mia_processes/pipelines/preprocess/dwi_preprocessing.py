# -*- coding: utf-8 -*-

"""The pipelines library of the mia_processes package.

The objective of this module is to propose pipelines built from atomic bricks
or other pipelines provided in the mia_processes library.

:Contains:
    :Class:
        - Dwi_preprocessing

"""
from capsul.api import Pipeline


class Dwi_preprocessing(Pipeline):
    """
    *DWI preprocessing with reverse b0*

    Please, see the complete documentation for the
    `Dwi_preprocessing pipeline in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Dwi_preprocessing.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "Convert_dwi",
            "mia_processes.bricks.preprocess.mrtrix.processes.MRConvert",
        )
        self.add_process(
            "Convert_b0",
            "mia_processes.bricks.preprocess.mrtrix.processes.MRConvert",
        )
        self.add_process(
            "Denoising",
            "mia_processes.bricks.preprocess.mrtrix.processes.DWIDenoise",
        )
        self.add_process(
            "Unringing",
            "mia_processes.bricks.preprocess.mrtrix.processes.MRDeGibbs",
        )
        self.add_process(
            "MotionDistortionCorrection",
            "mia_processes.bricks.preprocess.mrtrix.processes.DWIPreproc",
        )
        self.nodes["MotionDistortionCorrection"].process.rpe_options = (
            "rpe_pair"
        )
        self.nodes["MotionDistortionCorrection"].process.eddy_options = (
            " --slm=linear --data_is_shelled"
        )
        self.add_process(
            "Extractb0",
            "mia_processes.bricks.preprocess.mrtrix.processes.DWIExtract",
        )
        self.add_process(
            "Meanb0", "mia_processes.bricks.preprocess.mrtrix.processes.MRMath"
        )
        self.add_process(
            "Meanb0Corr",
            "mia_processes.bricks.preprocess.mrtrix.processes.MRMath",
        )
        self.add_process(
            "FilesToList", "mia_processes.bricks.tools.tools.Files_To_List"
        )
        self.add_process(
            "Concatenateb0",
            "mia_processes.bricks.preprocess.mrtrix.processes.MRCat",
        )
        self.add_process(
            "BiasCorrection",
            "mia_processes.bricks.preprocess.mrtrix.processes.DWIBiasCorrect",
        )
        self.add_process(
            "BrainMask",
            "mia_processes.bricks.preprocess.mrtrix.processes.DWIBrainMask",
        )

        # links
        self.export_parameter(
            "Convert_dwi", "in_file", "in_dwi", is_optional=False
        )
        self.export_parameter(
            "MotionDistortionCorrection",
            "pe_dir",
            "in_dwi_pe_dir",
            is_optional=True,
        )
        self.export_parameter(
            "MotionDistortionCorrection",
            "ro_time",
            "in_dwi_ro_time",
            is_optional=True,
        )
        self.export_parameter(
            "Convert_b0", "in_file", "in_b0_reverse", is_optional=False
        )
        self.add_link("Convert_dwi.out_file->Denoising.in_file")
        self.add_link("Convert_b0.out_file->Meanb0Corr.in_file")
        self.add_link("Denoising.out_file->Unringing.in_file")
        self.add_link("Unringing.out_file->MotionDistortionCorrection.in_file")
        self.add_link("Unringing.out_file->Extractb0.in_file")
        self.add_link(
            "MotionDistortionCorrection.out_file->BiasCorrection.in_file"
        )
        self.export_parameter(
            "BiasCorrection",
            "out_file",
            "preproc_dwi",
            is_optional=False,
        )
        self.add_link("Extractb0.out_file->Meanb0.in_file")
        self.add_link("Meanb0.out_file->FilesToList.file1")
        self.add_link("Meanb0Corr.out_file->FilesToList.file2")
        self.add_link("FilesToList.file_list->Concatenateb0.in_files")
        self.add_link(
            "Concatenateb0.out_file->MotionDistortionCorrection.se_epi_corr"
        )
        self.add_link("BiasCorrection.out_file->BrainMask.in_file")
        self.export_parameter(
            "BrainMask", "out_file", "brain_mask", is_optional=False
        )

        # parameters order

        self.reorder_traits(
            (
                "in_dwi",
                "in_dwi_pe_dir",
                "in_dwi_ro_time",
                "in_b0_reverse",
                "brain_mask",
                "preproc_dwi",
            )
        )

        # nodes positions
        self.node_position = {
            "Convert_dwi": (-1099.340793426031, -330.99389652344223),
            "Convert_b0": (-1079.6856105700367, 223.28226001559904),
            "inputs": (-1291.3993894677506, 48.14750084237923),
            "Denoising": (-817.46266529709, -344.2197355246564),
            "Unringing": (-642.6747160934633, -204.01988374578806),
            "MotionDistortionCorrection": (
                57.996932626461074,
                66.21320749222957,
            ),
            "Extractb0": (-451.7878747547151, -322.9089989467982),
            "Meanb0": (-259.73028467525273, -482.85714587016565),
            "Meanb0Corr": (-417.43894315496556, 209.10350773523913),
            "Concatenateb0": (-42.48238548560684, -207.5889444784554),
            "BiasCorrection": (220.9237026548114, -252.48423160549882),
            "BrainMask": (477.35300037914624, -69.03865707962856),
            "outputs": (526.6342301182356, 134.13224804042125),
            "FilesToList": (-199.17352902019883, -158.55231307808265),
        }

        # nodes dimensions
        self.node_dimension = {
            "Convert_dwi": (222.28125, 440.0),
            "Convert_b0": (222.28125, 440.0),
            "inputs": (147.046875, 250.0),
            "Denoising": (156.234375, 180.0),
            "Unringing": (135.75, 195.0),
            "MotionDistortionCorrection": (180.265625, 460.0),
            "Extractb0": (163.375, 195.0),
            "Meanb0": (186.96875, 160.0),
            "Meanb0Corr": (186.96875, 180.0),
            "Concatenateb0": (188.703125, 125.27466940582039),
            "BiasCorrection": (201.328125, 195.0),
            "BrainMask": (130.6875, 55.0),
            "outputs": (111.046875, 110.0),
            "FilesToList": (118.8125, 145.0),
        }

        self.do_autoexport_nodes_parameters = False
