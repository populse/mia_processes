# -*- coding: utf-8 -*-
"""The pipelines library of the mia_processes package.

The objective of this module is to propose pipelines built from atomic bricks
or other pipelines provided in the mia_processes library.

:Contains:
    :Class:
        - Dwi_tissue_boundaries

"""

from capsul.api import Pipeline


class Dwi_tissue_boundaries(Pipeline):

    """
    *DWI create tissue boudadiries with a T1w*

    Please, see the complete documentation for the
    `Dwi_tissue_boundaries pipeline in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/preprocess/Dwi_tissue_boundaries.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "ConvertT1w",
            "mia_processes.bricks.preprocess.mrtrix.processes.MRConvert",
        )
        self.add_process(
            "Generate5tt",
            "mia_processes.bricks.preprocess.mrtrix.processes.Generate5ttfsl",
        )
        self.add_process(
            "Extractb0",
            "mia_processes.bricks.preprocess.mrtrix.processes.DWIExtract",
        )
        self.add_process(
            "Meanb0", "mia_processes.bricks.preprocess.mrtrix.processes.MRMath"
        )
        self.add_process(
            "CoregistrationFlirt",
            "mia_processes.bricks.preprocess.fsl.processes.Flirt",
        )
        self.nodes["CoregistrationFlirt"].process.dof = 6
        self.nodes["CoregistrationFlirt"].process.interp = "nearestneighbour"
        self.add_process(
            "Transformfslconvert",
            "mia_processes.bricks.preprocess."
            "mrtrix.processes.TransformFSLConvert",
        )
        self.add_process(
            "ConvertMeanb0NIfTI",
            "mia_processes.bricks.preprocess.mrtrix.processes.MRConvert",
        )
        self.nodes["ConvertMeanb0NIfTI"].process.out_file_format = "NIFTI"
        self.add_process(
            "Transform_1",
            "mia_processes.bricks.preprocess.mrtrix.processes.MRTransform",
        )
        self.nodes["Transform_1"].process.inverse = True
        self.add_process(
            "Transform_2",
            "mia_processes.bricks.preprocess.mrtrix.processes.MRTransform",
        )
        self.nodes["Transform_2"].process.inverse = True
        self.add_process(
            "Generate5ttseed",
            "mia_processes.bricks.preprocess."
            "mrtrix.processes.Generate5tt2gmwmi",
        )

        # links
        self.export_parameter(
            "ConvertT1w", "in_file", "in_T1w", is_optional=False
        )
        self.add_link("in_T1w->CoregistrationFlirt.in_reference_file")
        self.export_parameter(
            "Extractb0", "in_file", "in_dwi", is_optional=False
        )
        self.export_parameter(
            "Generate5ttseed", "out_file", "gmwmSeed_coreg", is_optional=False
        )
        self.add_link("Extractb0.out_file->Meanb0.in_file")
        self.add_link("Meanb0.out_file->Transformfslconvert.in_file")
        self.add_link("Meanb0.out_file->ConvertMeanb0NIfTI.in_file")
        self.add_link("ConvertT1w.out_file->Generate5tt.in_file")
        self.add_link("ConvertT1w.out_file->Transform_1.in_file")
        self.add_link("ConvertT1w.out_file->Transformfslconvert.reference")
        self.add_link(
            "CoregistrationFlirt.out_matrix_file->"
            "Transformfslconvert.in_transform"
        )
        self.add_link(
            "Transformfslconvert.out_transform->Transform_2.linear_transform"
        )
        self.add_link(
            "Transformfslconvert.out_transform->Transform_1.linear_transform"
        )
        self.add_link(
            "ConvertMeanb0NIfTI.out_file->CoregistrationFlirt.in_file"
        )
        self.export_parameter(
            "Transform_1", "out_file", "T1w_coreg", is_optional=True
        )
        self.export_parameter(
            "Transform_2", "out_file", "5tt_coreg", is_optional=True
        )
        self.add_link("Transform_2.out_file->Generate5ttseed.in_file")
        self.add_link("Generate5tt.out_file->Transform_2.in_file")

        # parameters order

        self.reorder_traits(
            (
                "in_T1w",
                "in_dwi",
                "T1w_coreg",
                "5tt_coreg",
                "gmwmSeed_coreg",
            )
        )

        # nodes positions
        self.node_position = {
            "Generate5ttseed": (572.7923970930416, 563.1519654602814),
            "Extractb0": (-1005.4957981923569, 417.64372749745445),
            "Meanb0": (-806.8555420941282, 463.51008970621),
            "ConvertT1w": (-832.3275799578836, -255.8711697726233),
            "CoregistrationFlirt": (-290.10564078264713, 750.5788923698742),
            "Transformfslconvert": (-249.62126489965024, 497.1893676200373),
            "ConvertMeanb0NIfTI": (-585.5133725344976, 705.5775392569157),
            "Transform_1": (158.362335819779, -327.46138729475047),
            "Transform_2": (97.9136867940988, 468.4686461045751),
            "inputs": (-1161.5447656473912, 66.84405832847222),
            "Generate5tt": (-371.1225123162598, -137.1539719429656),
            "outputs": (791.427723268604, 42.77366014726516),
        }

        # nodes dimensions
        self.node_dimension = {
            "ConvertT1w": (228.28125, 460.0),
            "dwidenoise_1": (156.234375, 180.0),
            "Generate5ttseed": (153.359375, 110.0),
            "Extractb0": (166.375, 215.0),
            "Meanb0": (186.96875, 180.0),
            "CoregistrationFlirt": (257.890625, 215.0),
            "Transformfslconvert": (226.1875, 145.0),
            "ConvertMeanb0NIfTI": (222.28125, 460.0),
            "Transform_1": (220.34375, 600.0),
            "Transform_2": (220.34375, 600.0),
            "inputs": (106.984375, 145.0),
            "Generate5tt": (199.46875, 250.0),
            "outputs": (148.671875, 145.0),
        }

        self.do_autoexport_nodes_parameters = False
