# -*- coding: utf-8 -*-

"""
A pipeline to get average signals in severals ROI for functional MRI data.
It used Assemblynet segmentation

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################


from capsul.api import Pipeline


class Bold_signals_roi(Pipeline):
    """
    *Get plots of the BOLD average signal in severals ROI
    (defined by AssemblyNet)*

    Please, see the complete documentation for the `Bold_signals_roi brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/qualiTyControl/Bold_signals_roi.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "assemblynetdocker_1",
            "mia_processes.bricks.preprocess."
            "volbrain.processes.AssemblyNetDocker",
        )
        self.add_process(
            "epireg_1", "mia_processes.bricks.preprocess.fsl.processes.EpiReg"
        )
        self.add_process(
            "inverse_transfo",
            "mia_processes.bricks.preprocess.fsl.processes.ConvertXFM",
        )
        self.nodes["inverse_transfo"].process.invert_xfm = True
        self.add_process(
            "flirt_native_tissues",
            "mia_processes.bricks.preprocess.fsl.processes.Flirt",
        )
        self.nodes["flirt_native_tissues"].process.apply_xfm = True
        self.nodes["flirt_native_tissues"].process.interp = "nearestneighbour"
        self.add_process(
            "flirt_native_structures",
            "mia_processes.bricks.preprocess.fsl.processes.Flirt",
        )
        self.nodes["flirt_native_structures"].process.apply_xfm = True
        self.nodes["flirt_native_structures"].process.interp = (
            "nearestneighbour"
        )
        self.add_process(
            "flirt_native_macrostructures",
            "mia_processes.bricks.preprocess.fsl.processes.Flirt",
        )
        self.nodes["flirt_native_macrostructures"].process.apply_xfm = True
        self.nodes["flirt_native_macrostructures"].process.interp = (
            "nearestneighbour"
        )
        self.add_process(
            "flirt_native_lobes",
            "mia_processes.bricks.preprocess.fsl.processes.Flirt",
        )
        self.nodes["flirt_native_lobes"].process.apply_xfm = True
        self.nodes["flirt_native_lobes"].process.interp = "nearestneighbour"
        self.add_process(
            "extractsignalroi_tissues",
            "mia_processes.bricks.preprocess."
            "others.processing.ExtractSignalROI",
        )
        self.nodes["extractsignalroi_tissues"].process.suffix = "tissues"
        self.add_process(
            "plotsignalroi_tissues",
            "mia_processes.bricks.reports.processes.PlotSignalROI",
        )
        self.nodes["plotsignalroi_tissues"].process.suffix = "tissues"
        self.add_process(
            "getlabels_tissues",
            "mia_processes.bricks.preprocess.volbrain.processes.GetLabels",
        )
        self.nodes["getlabels_tissues"].process.tissues = True
        self.add_process(
            "getlabels_macrostructures",
            "mia_processes.bricks.preprocess.volbrain.processes.GetLabels",
        )
        self.nodes["getlabels_macrostructures"].process.macrostructures = True
        self.add_process(
            "extractsignalroi_macrostructures",
            "mia_processes.bricks.preprocess."
            "others.processing.ExtractSignalROI",
        )
        self.nodes["extractsignalroi_macrostructures"].process.suffix = (
            "macrostructures"
        )
        self.add_process(
            "plotsignalroi_macrostructures",
            "mia_processes.bricks.reports.processes.PlotSignalROI",
        )
        self.nodes["plotsignalroi_macrostructures"].process.suffix = (
            "macrostructures"
        )
        self.add_process(
            "plotsignalroi_lobes",
            "mia_processes.bricks.reports.processes.PlotSignalROI",
        )
        self.nodes["plotsignalroi_lobes"].process.suffix = "lobes"
        self.add_process(
            "extractsignalroi_lobes",
            "mia_processes.bricks.preprocess."
            "others.processing.ExtractSignalROI",
        )
        self.nodes["extractsignalroi_lobes"].process.suffix = "lobes"
        self.add_process(
            "getlabels_lobes",
            "mia_processes.bricks.preprocess.volbrain.processes.GetLabels",
        )
        self.nodes["getlabels_lobes"].process.lobes = True
        self.add_process(
            "extractroibylabel_structures",
            "mia_processes.bricks.preprocess."
            "others.processing.ExtractROIbyLabel",
        )
        self.add_process(
            "extractsignalroi_structures",
            "mia_processes.bricks.preprocess."
            "others.processing.ExtractSignalROI",
        )
        self.add_process(
            "plotsignalroi_structures",
            "mia_processes.bricks.reports.processes.PlotSignalROI",
        )
        self.add_process(
            "labelscorrespondence_structures",
            "mia_processes.bricks.preprocess."
            "volbrain.processes.LabelsCorrespondence",
        )
        self.nodes["labelscorrespondence_structures"].process.structures = True
        self.add_process(
            "extractsignalroi_1",
            "mia_processes.bricks.preprocess."
            "others.processing.ExtractSignalROI",
        )
        self.add_process(
            "automask_bold",
            "mia_processes.bricks.preprocess.afni.processes.Automask",
        )
        self.add_process(
            "betsurfacesextraction_1",
            "mia_processes.bricks.preprocess."
            "fsl.processes.BetSurfacesExtraction",
        )

        # links
        self.export_parameter(
            "assemblynetdocker_1", "in_file", "anat", is_optional=False
        )
        self.add_link("anat->epireg_1.in_t1")
        self.add_link("anat->betsurfacesextraction_1.in_file")
        self.export_parameter(
            "extractsignalroi_structures", "in_file", "bold", is_optional=False
        )
        self.add_link("bold->flirt_native_tissues.in_reference_file")
        self.add_link("bold->flirt_native_macrostructures.in_reference_file")
        self.add_link("bold->flirt_native_structures.in_reference_file")
        self.add_link("bold->extractsignalroi_tissues.in_file")
        self.add_link("bold->epireg_1.in_epi")
        self.add_link("bold->plotsignalroi_tissues.in_file")
        self.add_link("bold->extractsignalroi_1.in_file")
        self.add_link("bold->plotsignalroi_structures.in_file")
        self.add_link("bold->extractsignalroi_macrostructures.in_file")
        self.add_link("bold->flirt_native_lobes.in_reference_file")
        self.add_link("bold->automask_bold.in_file")
        self.add_link("bold->plotsignalroi_lobes.in_file")
        self.add_link("bold->plotsignalroi_macrostructures.in_file")
        self.add_link("bold->extractsignalroi_lobes.in_file")
        self.export_parameter(
            "extractsignalroi_structures",
            "labels",
            "labels_structures",
            is_optional=False,
        )
        self.add_link("labels_structures->extractroibylabel_structures.labels")
        self.add_link(
            "labels_structures->labelscorrespondence_structures.labels_names"
        )
        self.add_link("labels_structures->plotsignalroi_structures.labels")
        self.add_link(
            "assemblynetdocker_1.native_structures->"
            "flirt_native_structures.in_file"
        )
        self.add_link(
            "assemblynetdocker_1.native_lobes->flirt_native_lobes.in_file"
        )
        self.add_link(
            "assemblynetdocker_1.native_macrostructures->"
            "flirt_native_macrostructures.in_file"
        )
        self.add_link(
            "assemblynetdocker_1.native_tissues->flirt_native_tissues.in_file"
        )
        self.export_parameter(
            "epireg_1", "out_file", "bold_struc", is_optional=False
        )
        self.add_link("epireg_1.epi2str_mat->inverse_transfo.in_transfo")
        self.add_link(
            "inverse_transfo.out_file->flirt_native_lobes.in_matrix_file"
        )
        self.add_link(
            "inverse_transfo.out_file->flirt_native_structures.in_matrix_file"
        )
        self.add_link(
            "inverse_transfo.out_file->flirt_native_tissues.in_matrix_file"
        )
        self.add_link(
            "inverse_transfo.out_file->"
            "flirt_native_macrostructures.in_matrix_file"
        )
        self.add_link(
            "flirt_native_tissues.out_file->plotsignalroi_tissues.rois_files"
        )
        self.add_link(
            "flirt_native_tissues.out_file->extractsignalroi_tissues.in_seg"
        )
        self.add_link(
            "flirt_native_structures.out_file->"
            "extractroibylabel_structures.in_file"
        )
        self.add_link(
            "flirt_native_macrostructures.out_file->"
            "plotsignalroi_macrostructures.rois_files"
        )
        self.add_link(
            "flirt_native_macrostructures.out_file->"
            "extractsignalroi_macrostructures.in_seg"
        )
        self.add_link(
            "flirt_native_lobes.out_file->plotsignalroi_lobes.rois_files"
        )
        self.add_link(
            "flirt_native_lobes.out_file->extractsignalroi_lobes.in_seg"
        )
        self.add_link(
            "extractsignalroi_tissues.signals->plotsignalroi_tissues.signals"
        )
        self.export_parameter(
            "plotsignalroi_tissues",
            "out_png",
            "out_png_tissues",
            is_optional=False,
        )
        self.add_link("getlabels_tissues.labels->plotsignalroi_tissues.labels")
        self.add_link(
            "getlabels_tissues.labels->extractsignalroi_tissues.labels"
        )
        self.add_link(
            "getlabels_tissues.names->plotsignalroi_tissues.labels_names"
        )
        self.add_link(
            "getlabels_macrostructures.labels->"
            "extractsignalroi_macrostructures.labels"
        )
        self.add_link(
            "getlabels_macrostructures.labels->"
            "plotsignalroi_macrostructures.labels"
        )
        self.add_link(
            "getlabels_macrostructures.names->"
            "plotsignalroi_macrostructures.labels_names"
        )
        self.add_link(
            "extractsignalroi_macrostructures.signals->"
            "plotsignalroi_macrostructures.signals"
        )
        self.export_parameter(
            "plotsignalroi_macrostructures",
            "out_png",
            "out_png_macrostructures",
            is_optional=False,
        )
        self.export_parameter(
            "plotsignalroi_lobes",
            "out_png",
            "out_png_lobes",
            is_optional=False,
        )
        self.add_link(
            "extractsignalroi_lobes.signals->plotsignalroi_lobes.signals"
        )
        self.add_link("getlabels_lobes.labels->plotsignalroi_lobes.labels")
        self.add_link("getlabels_lobes.labels->extractsignalroi_lobes.labels")
        self.add_link(
            "getlabels_lobes.names->plotsignalroi_lobes.labels_names"
        )
        self.add_link(
            "extractroibylabel_structures.out_concate->"
            "plotsignalroi_structures.rois_files"
        )
        self.add_link(
            "extractroibylabel_structures.out_concate->"
            "extractsignalroi_structures.in_seg"
        )
        self.add_link(
            "extractsignalroi_structures.signals->"
            "plotsignalroi_structures.signals"
        )
        self.export_parameter(
            "plotsignalroi_structures",
            "out_png",
            "out_png_structures",
            is_optional=False,
        )
        self.add_link(
            "labelscorrespondence_structures.correspondence->"
            "plotsignalroi_structures.labels_names"
        )
        self.add_link(
            "extractsignalroi_1.signals->"
            "plotsignalroi_structures.signals_whole_brain"
        )
        self.add_link(
            "extractsignalroi_1.signals->"
            "plotsignalroi_macrostructures.signals_whole_brain"
        )
        self.add_link(
            "extractsignalroi_1.signals->"
            "plotsignalroi_lobes.signals_whole_brain"
        )
        self.add_link(
            "extractsignalroi_1.signals->"
            "plotsignalroi_tissues.signals_whole_brain"
        )
        self.add_link("automask_bold.out_file->extractsignalroi_1.in_seg")
        self.add_link("betsurfacesextraction_1.out_file->epireg_1.in_t1_brain")

        # parameters order
        self.reorder_traits(
            (
                "anat",
                "bold",
                "out_png_tissues",
                "out_png_macrostructures",
                "out_png_lobes",
                "labels_structures",
                "out_png_structures",
                "bold_struc",
            )
        )

        # default and initial values
        self.labels_structures = [47, 48]

        # nodes positions
        self.node_position = {
            "assemblynetdocker_1": (-708.9103365730818, -1499.4914764535704),
            "inputs": (-905.9274258262157, -852.7262749946995),
            "epireg_1": (-542.6382655050916, -733.2514625397064),
            "inverse_transfo": (-255.60423887638504, -566.3713190537865),
            "flirt_native_tissues": (86.93184000000002, -1670.112776083243),
            "flirt_native_structures": (
                93.74415922646938,
                -1201.0890483997437,
            ),
            "flirt_native_macrostructures": (
                115.71027370446416,
                -759.6737487651787,
            ),
            "flirt_native_lobes": (134.35676364535732, -279.0722046052042),
            "outputs": (1242.0700112287686, -1575.5108974205662),
            "extractsignalroi_tissues": (
                571.8801574771767,
                -1670.193042698244,
            ),
            "plotsignalroi_tissues": (880.7829559300474, -1571.119436668161),
            "getlabels_tissues": (503.557355746152, -1392.5463280323697),
            "getlabels_macrostructures": (
                576.0170895497354,
                -478.59695305381115,
            ),
            "extractsignalroi_macrostructures": (
                587.3238713199722,
                -663.7653326408645,
            ),
            "plotsignalroi_macrostructures": (
                930.2447436556096,
                -706.685926228697,
            ),
            "plotsignalroi_lobes": (975.8567439896908, 102.6270007516834),
            "extractsignalroi_lobes": (712.9875841695895, 149.43931688403023),
            "getlabels_lobes": (554.5458987985697, 379.89995015096827),
            "extractroibylabel_structures": (
                574.5850017156421,
                -1166.2423517372858,
            ),
            "extractsignalroi_structures": (
                956.5161374949022,
                -1226.5533343129082,
            ),
            "plotsignalroi_structures": (
                1258.6736893354432,
                -1053.6577219194712,
            ),
            "labelscorrespondence_structures": (
                602.1988665783529,
                -919.9620458713234,
            ),
            "extractsignalroi_1": (-310.70828880160025, -15.389216879722653),
            "automask_bold": (-789.115951784914, -55.278067031963815),
            "betsurfacesextraction_1": (-931.041544325601, -632.5396751525075),
        }

        # nodes dimensions
        self.node_dimension = {
            "assemblynetdocker_1": (246.125, 530.0),
            "inputs": (171.1256372392954, 145.0),
            "epireg_1": (238.359375, 530.0),
            "convertxfm_1": (192.203125, 250.0),
            "inverse_transfo": (193.125, 250.0),
            "flirt_1": (278.046875, 1580.0),
            "flirt_native_tissues": (279.390625, 1580.0),
            "flirt_native_structures": (279.390625, 1580.0),
            "flirt_2": (278.046875, 1580.0),
            "flirt_native_macrostructures": (279.390625, 1580.0),
            "flirt_native_lobes": (279.390625, 1580.0),
            "plotsignalroi_1": (232.3125, 215.0),
            "extractsignalroi_1": (126.1875, 180.0),
            "outputs": (210.60242003946317, 215.0),
            "labelscorrespondence_1": (261.828125, 215.0),
            "extractsignalroi_tissues": (157.828125, 180.0),
            "plotsignalroi_tissues": (234.28125, 285.0),
            "getlabels_1": (195.265625, 180.0),
            "labelscorrespondence_tissues": (258.828125, 215.0),
            "getlabels_tissues": (195.265625, 180.0),
            "getlabels_macrostructures": (195.265625, 180.0),
            "extractsignalroi_macrostructures": (221.921875, 180.0),
            "plotsignalroi_macrostructures": (234.28125, 285.0),
            "plotsignalroi_lobes": (234.28125, 285.0),
            "extractsignalroi_lobes": (148.6875, 180.0),
            "getlabels_lobes": (195.265625, 180.0),
            "extractroibylabel_1": (238.96875, 180.0),
            "extractroibylabel_structures": (240.953125, 215.0),
            "extractsignalroi_structures": (181.125, 180.0),
            "plotsignalroi_structures": (234.28125, 285.0),
            "labelscorrespondence_structures": (263.421875, 215.0),
            "automask_1": (217.765625, 285.0),
            "automask_bold": (217.765625, 285.0),
            "betsurfacesextraction_1": (253.96875, 320.0),
        }

        self.do_autoexport_nodes_parameters = False
