# -*- coding: utf-8 -*-

from capsul.api import Pipeline
import traits.api as traits


class Bold_mriqc_pipeline(Pipeline):

    def pipeline_definition(self):
        # nodes
        self.add_process("bold_iqms_pipeline",
                         "mia_processes.pipelines.reports.bold_iqms_pipeline."
                         "Bold_iqms_pipeline")
        self.nodes["bold_iqms_pipeline"].process.nodes_activation = {
            'outliercount': True,
            'boldiqms': True,
            'carpetparcellation': True,
            'framewisedisplacement': True,
            'spikes': True,
            'gcor': True,
            'computedvars': True,
            'fwhmx': True,
            'qualityindex': True}
        self.add_process("nonsteadystatedetector",
                         "mia_processes.bricks.preprocess.others.processing."
                         "NonSteadyStateDetector")
        self.add_process("sanitize",
                         "mia_processes.bricks.preprocess.others.processing."
                         "Sanitize")
        self.add_process("tsnr",
                         "mia_processes.bricks.preprocess.others.processing."
                         "TSNR")
        self.add_process("mean",
                         "mia_processes.bricks.preprocess.afni.processes."
                         "Mean")
        self.add_process("automask",
                         "mia_processes.bricks.preprocess.afni.processes."
                         "Automask")
        self.add_process("bold_hmc_pipeline",
                         "mia_processes.pipelines.preprocess.bold_hmc_pipeline."
                         "Bold_hmc_pipeline")
        self.nodes["bold_hmc_pipeline"].process.nodes_activation = {
            'despike': False,
            'deoblique': False,
            'volreg': True}
        self.add_process("bold_mni_align",
                         "mia_processes.pipelines.preprocess.bold_mni_align."
                         "Bold_mni_align")
        self.nodes["bold_mni_align"].set_plug_value("epi_mask",
                                                    traits.Undefined)
        self.nodes["bold_mni_align"].process.nodes[
            "registration"].set_plug_value("moving_image_masks",
                                           traits.Undefined)
        self.nodes["bold_mni_align"].process.nodes_activation = {
            'affineinitializer': True,
            'registration': True,
            'n4biasfieldcorrection': True,
            'applytransforms': True,
            'template': True,
            'template_mask': True,
            'seg_template': True}
        self.nodes["bold_mni_align"].process.epi_mask = traits.Undefined
        self.add_process("mriqc_func_report",
                         "mia_processes.bricks.reports.reporting."
                         "MRIQC_func_report")

        # links
        self.export_parameter("bold_hmc_pipeline", "despike",
                              pipeline_parameter="hmc_despike")
        self.export_parameter("bold_hmc_pipeline", "deoblique",
                              pipeline_parameter="hmc_deoblique")
        self.export_parameter("sanitize", "in_file",
                              pipeline_parameter="func_file",
                              is_optional=False)
        self.add_link("func_file->nonsteadystatedetector.in_file")
        self.add_link("func_file->mriqc_func_report.func")
        self.add_link("bold_iqms_pipeline.BoldQC_out_file->"
                      "mriqc_func_report.IQMs_file")
        self.add_link("nonsteadystatedetector.n_volumes_to_discard->"
                      "sanitize.n_volumes_to_discard")
        self.add_link("nonsteadystatedetector.n_volumes_to_discard->"
                      "bold_iqms_pipeline.dummy_TRs")
        self.add_link("sanitize.out_file->bold_iqms_pipeline.ras_epi")
        self.add_link("sanitize.out_file->bold_hmc_pipeline.in_file")
        self.add_link("tsnr.out_tsnr_file->bold_iqms_pipeline.epi_tsnr")
        self.add_link("mean.out_file->bold_mni_align.epi_mean")
        self.add_link("mean.out_file->automask.in_file")
        self.add_link("mean.out_file->bold_iqms_pipeline.epi_mean")
        self.add_link("automask.out_file->bold_mni_align.epi_mask")
        self.add_link("automask.out_file->bold_iqms_pipeline.brainmask")
        self.add_link("bold_hmc_pipeline.out_file->tsnr.in_file")
        self.add_link("bold_hmc_pipeline.out_file->mean.in_file")
        self.add_link("bold_hmc_pipeline.out_file->bold_iqms_pipeline.hmc_epi")
        self.add_link("bold_hmc_pipeline.oned_file->"
                      "bold_iqms_pipeline.hmc_motion")
        self.add_link("bold_mni_align.epi_parc->bold_iqms_pipeline.epi_parc")
        self.add_link("bold_mni_align.epi_mni->mriqc_func_report.norm_func")
        self.export_parameter("mriqc_func_report", "report",
                              pipeline_parameter="func_report",
                              is_optional=True)
        self.export_parameter("bold_iqms_pipeline", "carpet_seg",
                              is_optional=True)

        # parameters order

        self.reorder_traits(("func_file", "carpet_seg", "func_report"))

        # default and initial values

        # nodes positions
        self.node_position = {
            "bold_iqms_pipeline": (842.4679035022216, 159.92981034666673),
            "nonsteadystatedetector": (-717.1107544384864, 113.59195862760868),
            "sanitize": (-234.51652984508183, 346.43391671594577),
            "tsnr": (-127.7391993962027, 60.050718720000134),
            "mean": (-173.11641600000007, -181.6966963199999),
            "automask": (144.1320191999998, -189.81009407999983),
            "bold_hmc_pipeline": (-527.8792115199997, -149.2405196799998),
            "bold_mni_align": (283.31445431594545, 116.08725895240536),
            "inputs": (-890.4021152762024, 569.7810895696218),
            "outputs": (1530.1459691893501, 147.22476292910784),
            "mriqc_func_report": (1282.5825850836268, 1058.7313069644374),
        }

        # nodes dimensions
        self.node_dimension = {
            "bold_iqms_pipeline": (299.734375, 810.0),
            "nonsteadystatedetector": (194.09375, 75.0),
            "sanitize": (202.453125, 215.0),
            "tsnr": (194.03125, 215.0),
            "mean": (146.0, 145.0),
            "automask": (146.0, 145.0),
            "bold_hmc_pipeline": (154.453125, 215.0),
            "bold_mni_align": (377.796875, 670.0),
            "inputs": (243.9375, 1265.0),
            "outputs": (203.16492003946317, 355.0),
            "mriqc_func_report": (210.8125, 425.0),
        }

        self.do_autoexport_nodes_parameters = False