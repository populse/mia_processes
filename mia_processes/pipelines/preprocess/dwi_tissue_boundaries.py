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
    `Dwi_tissue_boundaries pipeline in the mia_processes website
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
            "Convert5ttNIfTI",
            "mia_processes.bricks.preprocess.mrtrix.processes.MRConvert",
        )
        self.nodes["Convert5ttNIfTI"].process.out_file_format = "NIFTI"
        self.add_process(
            "ExtractROI",
            "mia_processes.bricks.preprocess.fsl.processes.ExtractROI",
        )
        self.nodes["ExtractROI"].process.t_min = 0
        self.nodes["ExtractROI"].process.t_size = 1
        self.nodes["ExtractROI"].process.suffix = "gm"
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
        self.add_link("Generate5tt.out_file->Convert5ttNIfTI.in_file")
        self.add_link("Convert5ttNIfTI.out_file->ExtractROI.in_file")
        self.add_link(
            "ExtractROI.roi_file->CoregistrationFlirt.in_reference_file"
        )
        self.export_parameter(
            "Extractb0", "in_file", "in_dwi", is_optional=False
        )
        self.export_parameter(
            "Generate5ttseed", "out_file", "gmwmSeed_coreg", is_optional=False
        )
        self.add_link("Extractb0.out_file->Meanb0.in_file")
        self.add_link("Meanb0.out_file->ConvertMeanb0NIfTI.in_file")
        self.add_link("ConvertT1w.out_file->Generate5tt.in_file")
        self.add_link("ConvertT1w.out_file->Transform_1.in_file")
        self.add_link("ExtractROI.roi_file->Transformfslconvert.reference")
        self.add_link(
            "CoregistrationFlirt.out_matrix_file->"
            "Transformfslconvert.in_transform"
        )
        self.add_link(
            "Transformfslconvert.out_transform->Transform_2.linear_transform"
        )
        self.add_link(
            "ConvertMeanb0NIfTI.out_file->Transformfslconvert.in_file"
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
        # nodes positions
        self.node_position = {
            "Generate5ttseed": (733.5035929823755, 496.18896717305904),
            "Extractb0": (-1013.5313579868234, 229.70091230465),
            "Meanb0": (-804.1770221626393, 294.3169140338279),
            "ConvertT1w": (-989.0209959499841, -274.6208092930456),
            "CoregistrationFlirt": (-182.9648435230911, 517.5476583303403),
            "Transformfslconvert": (128.05004544028458, 371.2989308400588),
            "ConvertMeanb0NIfTI": (-574.7992928085422, 453.79666569695905),
            "Transform_1": (348.53725095549083, -292.6406281853948),
            "Transform_2": (371.1227198059666, 653.2865213773091),
            "inputs": (-1161.5447656473912, 66.84405832847222),
            "Generate5tt": (-717.0980034669092, -99.20827291353953),
            "outputs": (791.427723268604, 42.77366014726516),
            "Convert5ttNIfTI": (-409.79240138257603, -146.41518218711377),
            "ExtractROI": (-110.24074392890529, -180.34117876155892),
        }

        # nodes dimensions
        self.node_dimension = {
            "ConvertT1w": (222.28125, 460.0),
            "dwidenoise_1": (156.234375, 180.0),
            "Generate5ttseed": (146.65625, 110.0),
            "Extractb0": (163.375, 215.0),
            "Meanb0": (186.96875, 180.0),
            "CoregistrationFlirt": (278.046875, 1580.0),
            "Transformfslconvert": (226.1875, 145.0),
            "ConvertMeanb0NIfTI": (222.28125, 460.0),
            "Transform_1": (220.34375, 600.0),
            "Transform_2": (220.34375, 600.0),
            "inputs": (95.26626223929541, 90.0),
            "Generate5tt": (199.46875, 250.0),
            "outputs": (151.69617003946317, 115.0),
            "Convert5ttNIfTI": (222.28125, 460.0),
            "ExtractROI": (171.046875, 425.0),
        }

        self.do_autoexport_nodes_parameters = False
