# -*- coding: utf-8 -*-

"""
The atomic calculations from afni.

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from capsul.api import Pipeline


class Anat_mriqc(Pipeline):
    """
    *Get no-reference IQMs (image quality metrics) from structural (T1w and
    T2w) data using mriqc anatomical workflow (mriqc v22.06)*

    Please, see the complete documentation for the `Anat_mriqc_pipeline brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/qualiTyControl/Anat_mriqc.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "conformimage",
            "mia_processes.bricks.preprocess."
            "others.processing.ConformImage",
        )
        self.add_process(
            "anat_skullstrip_pipeline",
            "mia_processes.pipelines.preprocess."
            "anat_skullstrip_synthstrip."
            "Anat_skullstrip_synthstrip",
        )
        self.add_process(
            "segment",
            "mia_processes.bricks.preprocess." "fsl.processes.FastSegment",
        )
        self.add_process(
            "anat_spatial_norm",
            "mia_processes.pipelines.preprocess."
            "anat_spatial_norm.Anat_spatial_norm",
        )
        self.add_process(
            "anat_headmask",
            "mia_processes.bricks.preprocess" ".fsl.BetSurfacesExtraction",
        )
        self.add_process(
            "anat_airmask",
            "mia_processes.pipelines.preprocess." "anat_airmask.Anat_airmask",
        )
        self.add_process(
            "anat_mni_tpms_pipeline",
            "mia_processes.pipelines.preprocess."
            "anat_mni_tpms.Anat_mni_tpms",
        )
        self.add_process(
            "anatiqms", "mia_processes.bricks.reports.processes.AnatIQMs"
        )
        self.add_process(
            "harmonize",
            "mia_processes.bricks.preprocess." "others.processing.Harmonize",
        )
        self.add_process("fwhmx", "mia_processes.bricks.preprocess.afni.FWHMx")
        self.nodes["fwhmx"].plugs["detrend"].optional = True
        self.add_process(
            "list_to_file", "mia_processes.bricks.tools.tools.List_To_File"
        )
        self.add_process(
            "mriqc_anat_report",
            "mia_processes.bricks.reports." "reporting.ReportAnatMriqc",
        )

        # links
        self.export_parameter(
            "conformimage", "in_file", "anat_file", is_optional=False
        )
        self.add_link("anat_file->mriqc_anat_report.anat")
        self.add_link(
            "conformimage.out_file->" "anat_skullstrip_pipeline.in_file"
        )
        self.add_link("conformimage.out_file->anatiqms.in_ras")
        self.add_link("conformimage.out_file->fwhmx.in_file")
        self.add_link("conformimage.out_file->anat_mni_tpms_pipeline.in_ras")
        self.add_link("conformimage.out_file->anat_airmask.in_file")
        self.add_link("anat_skullstrip_pipeline.out_brain->segment.in_file")
        self.add_link(
            "anat_skullstrip_pipeline.out_mask_synthstrip->" "fwhmx.mask_file"
        )
        self.add_link(
            "anat_skullstrip_pipeline.out_mask_synthstrip->"
            "anat_spatial_norm.moving_mask"
        )
        self.add_link(
            "anat_skullstrip_pipeline.out_mask_synthstrip->"
            "anat_airmask.in_mask"
        )
        self.add_link(
            "anat_skullstrip_pipeline.out_mask_synthstrip->"
            "mriqc_anat_report.brain_mask"
        )
        self.add_link(
            "anat_skullstrip_pipeline.out_corrected" "->harmonize.in_file"
        )
        self.add_link(
            "anat_skullstrip_pipeline.out_corrected->" "anat_headmask.in_file"
        )
        self.add_link(
            "anat_skullstrip_pipeline.out_corrected->"
            "anat_spatial_norm.moving_image"
        )
        self.add_link("anat_skullstrip_pipeline.bias_image->anatiqms.in_inu")
        self.add_link("segment.tissue_class_map->anatiqms.segmentation")
        self.add_link(
            "segment.tissue_class_map->" "mriqc_anat_report.segmentation"
        )
        self.add_link("segment.partial_volume_files->list_to_file.file_list")
        self.add_link("segment.partial_volume_files->anatiqms.pvms")
        self.add_link(
            "anat_spatial_norm.inverse_composite_transform"
            "->anat_mni_tpms_pipeline.inverse_composite_transform"
        )
        self.add_link(
            "anat_spatial_norm.inverse_composite_transform->"
            "anat_airmask.inverse_composite_transform"
        )
        self.add_link(
            "anat_spatial_norm.warped_image->" "mriqc_anat_report.norm_anat"
        )
        self.add_link(
            "anat_headmask.outskin_mask_file->" "anat_airmask.head_mask"
        )
        self.add_link("anat_headmask.outskin_mask_file->anatiqms.headmask")
        self.add_link("anat_airmask.out_hat_mask->anatiqms.hatmask")
        self.add_link("anat_airmask.out_art_mask->anatiqms.artmask")
        self.add_link("anat_airmask.out_air_mask->anatiqms.airmask")
        self.add_link("anat_airmask.out_rot_mask->anatiqms.rotmask")
        self.add_link(
            "anat_headmask.outskin_mask_file->" "mriqc_anat_report.head_mask"
        )
        self.add_link("anat_airmask.out_art_mask->mriqc_anat_report.art_mask")
        self.add_link("anat_airmask.out_air_mask->mriqc_anat_report.air_mask")
        self.add_link("anat_mni_tpms_pipeline.mni_tpms->anatiqms.mni_tpms")
        self.add_link("anatiqms.out_file->mriqc_anat_report.IQMs_file")
        self.add_link("harmonize.out_file->anatiqms.in_noinu")
        self.add_link("fwhmx.out_file->anatiqms.in_fwhm")
        self.add_link("list_to_file.file->harmonize.wm_mask")
        self.export_parameter(
            "mriqc_anat_report",
            "report",
            pipeline_parameter="anat_report",
            is_optional=False,
        )

        # parameters order
        self.reorder_traits(("anat_file", "anat_report"))

        # default and initial values
        self.nodes["fwhmx"].process.args = "-ShowMeClassicFWHM"
        self.nodes["conformimage"].process.suffix = " "
        self.nodes["conformimage"].process.prefix = " "
        self.nodes["list_to_file"].process.index_filter = [3]
        self.nodes["harmonize"].process.prefix = " "

        # nodes positions
        self.node_position = {
            "conformimage": (-266.294964300051, 665.2984755165728),
            "inputs": (-760.1707537491396, 306.60591462332025),
            "anat_skullstrip_pipeline": (-318.30359196600546, 123.0195907),
            "segment": (-82.8528, -26.389119999999984),
            "anat_spatial_norm": (122.1775299994049, 474.3922141023147),
            "anat_headmask": (239.62046769206347, 245.85271538452383),
            "anat_airmask": (512.8535576922617, -39.083212461309614),
            "anat_mni_tpms_pipeline": (23.313110769047682, 1180.4801720512571),
            "anatiqms": (878.0130611266952, 880.9613282881348),
            "harmonize": (-274.74608025667965, 1129.801461025529),
            "fwhmx": (565.5987063585978, 1172.0964170767852),
            "list_to_file": (-520.7743035899467, 1036.3461573847217),
            "outputs": (1266.8414282725082, 887.7104479139109),
            "mriqc_anat_report": (1081.212516565563, 151.23336869359628),
        }

        # nodes dimensions
        self.node_dimension = {
            "conformimage": (108.3125, 145.0),
            "inputs": (221.953125, 600.0),
            "anat_skullstrip_pipeline": (179.21875, 180.0),
            "segment": (223.25, 145.0),
            "anat_spatial_norm": (377.796875, 635.0),
            "anat_headmask": (176.640625, 110.0),
            "anat_airmask": (276.828125, 250.0),
            "anat_mni_tpms_pipeline": (256.109375, 215.0),
            "anatiqms": (152.203125, 460.0),
            "harmonize": (141.0, 215.0),
            "fwhmx": (133.734375, 215.0),
            "list_to_file": (117.75, 110.0),
            "outputs": (184.18054503946317, 250.0),
            "mriqc_anat_report": (218.484375, 425.0),
        }

        self.do_autoexport_nodes_parameters = False
