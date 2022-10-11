from capsul.api import Pipeline
import traits.api as traits


class Anat_mriqc_pipeline(Pipeline):

    def pipeline_definition(self):
        # nodes
        self.add_process("conformimage",
                         "mia_processes.bricks.preprocess.others.processing.ConformImage")
        self.add_process("anat_skullstrip_pipeline",
                         "mia_processes.pipelines.preprocess.anat_skullstrip_pipeline.Anat_skullstrip_pipeline")
        self.add_process("segment",
                         "mia_processes.bricks.preprocess.fsl.processes.Segment")
        self.add_process("anat_spatial_norm",
                         "mia_processes.pipelines.preprocess.anat_spatial_norm.Anat_spatial_norm")
        self.add_process("anat_headmask_pipeline",
                         "mia_processes.pipelines.preprocess.anat_headmask_pipeline.Anat_headmask_pipeline")
        self.add_process("anat_airmask_pipeline",
                         "mia_processes.pipelines.preprocess.anat_airmask_pipeline.Anat_airmask_pipeline")
        self.add_process("anat_mni_tpms_pipeline",
                         "mia_processes.pipelines.preprocess.anat_mni_tpms_pipeline.Anat_mni_tpms_pipeline")
        self.add_process("anatiqms",
                         "mia_processes.bricks.reports.processes.AnatIQMs")
        self.add_process("harmonize",
                         "mia_processes.bricks.preprocess.others.processing.Harmonize")
        self.add_process("fwhmx",
                         "mia_processes.bricks.reports.processes.FWHMx")
        self.nodes["fwhmx"].plugs["detrend"].optional = True
        self.add_process("list_to_file",
                         "mia_processes.bricks.tools.tools.List_To_File")
        self.add_process("mriqc_anat_report",
                         "mia_processes.bricks.reports.reporting.MRIQC_anat_report")

        # links
        self.export_parameter("conformimage", "in_file", "anat_file",
                              is_optional=False)
        self.add_link("anat_file->mriqc_anat_report.anat")
        self.add_link(
            "conformimage.out_file->anat_skullstrip_pipeline.in_file")
        self.add_link("conformimage.out_file->anatiqms.in_ras")
        self.add_link("conformimage.out_file->fwhmx.in_file")
        self.add_link("conformimage.out_file->anat_mni_tpms_pipeline.in_ras")
        self.add_link("conformimage.out_file->anat_airmask_pipeline.in_file")
        self.add_link("anat_skullstrip_pipeline.out_file->segment.in_file")
        self.add_link("anat_skullstrip_pipeline.out_mask->fwhmx.mask_file")
        self.add_link(
            "anat_skullstrip_pipeline.out_mask->anat_spatial_norm.moving_mask")
        self.add_link(
            "anat_skullstrip_pipeline.out_mask->anat_airmask_pipeline.in_mask")
        self.add_link(
            "anat_skullstrip_pipeline.bias_corrected->harmonize.in_file")
        self.add_link(
            "anat_skullstrip_pipeline.bias_corrected->anat_headmask_pipeline.in_file")
        self.add_link(
            "anat_skullstrip_pipeline.bias_corrected->anat_spatial_norm.moving_image")
        self.add_link("anat_skullstrip_pipeline.bias_image->anatiqms.in_inu")
        self.add_link(
            "segment.tissue_class_map->anat_headmask_pipeline.seg_file")
        self.add_link("segment.tissue_class_map->anatiqms.segmentation")
        self.add_link("segment.partial_volume_files->list_to_file.file_list")
        self.add_link("segment.partial_volume_files->anatiqms.pvms")
        self.add_link(
            "anat_spatial_norm.inverse_composite_transform->anat_mni_tpms_pipeline.inverse_composite_transform")
        self.add_link(
            "anat_spatial_norm.inverse_composite_transform->anat_airmask_pipeline.inverse_composite_transform")
        self.add_link(
            "anat_spatial_norm.warped_image->mriqc_anat_report.norm_anat")
        self.add_link(
            "anat_headmask_pipeline.out_file->anat_airmask_pipeline.head_mask")
        self.add_link("anat_headmask_pipeline.out_file->anatiqms.headmask")
        self.add_link("anat_airmask_pipeline.out_hat_mask->anatiqms.hatmask")
        self.add_link("anat_airmask_pipeline.out_art_mask->anatiqms.artmask")
        self.add_link("anat_airmask_pipeline.out_air_mask->anatiqms.airmask")
        self.add_link("anat_airmask_pipeline.out_rot_mask->anatiqms.rotmask")
        self.add_link("anat_mni_tpms_pipeline.mni_tpms->anatiqms.mni_tpms")
        self.add_link("anatiqms.out_file->mriqc_anat_report.IQMs_file")
        self.add_link("harmonize.out_file->anatiqms.in_noinu")
        self.add_link("fwhmx.out_file->anatiqms.in_fwhm")
        self.add_link("list_to_file.file->harmonize.wm_mask")
        self.export_parameter("mriqc_anat_report", "report", is_optional=False)

        # parameters order

        self.reorder_traits(("anat_file", "report"))

        # default and initial values
        self.nodes["conformimage"].process.suffix = ' '
        self.nodes["conformimage"].process.prefix = ' '
        self.nodes["list_to_file"].process.index_filter = [2]
        self.nodes["harmonize"].process.prefix = ' '

        # nodes positions
        self.node_position = {
            "conformimage": (-266.294964300051, 665.2984755165728),
            "inputs": (-760.1707537491396, 306.60591462332025),
            "anat_skullstrip_pipeline": (
            -318.30359196600546, 123.01959076904758),
            "segment": (-82.8528, -26.389119999999984),
            "anat_spatial_norm": (122.1775299994049, 474.3922141023147),
            "anat_headmask_pipeline": (239.62046769206347, 245.85271538452383),
            "anat_airmask_pipeline": (512.8535576922617, -39.083212461309614),
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
            "anat_headmask_pipeline": (176.640625, 110.0),
            "anat_airmask_pipeline": (276.828125, 250.0),
            "anat_mni_tpms_pipeline": (256.109375, 215.0),
            "anatiqms": (152.203125, 460.0),
            "harmonize": (141.0, 215.0),
            "fwhmx": (133.734375, 215.0),
            "list_to_file": (117.75, 110.0),
            "outputs": (184.18054503946317, 250.0),
            "mriqc_anat_report": (218.484375, 425.0),
        }

        self.do_autoexport_nodes_parameters = False
