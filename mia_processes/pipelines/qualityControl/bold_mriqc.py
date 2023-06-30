# -*- coding: utf-8 -*-

"""
The MRIQC pipeline for functional MRI data.

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

import traits.api as traits
from capsul.api import Pipeline


class Bold_mriqc(Pipeline):
    """
    *Get no-reference IQMs (image quality metrics) from functional MRI data
    using mriqc functional workflow (mriqc v22.06)*

    Please, see the complete documentation for the `Bold_mriqc brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/pipelines/qualiTyControl/Bold_mriqc.html>`_

    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "bold_iqms_pipeline",
            "mia_processes.pipelines.reports.bold_iqms." "Bold_iqms",
        )
        self.nodes["bold_iqms_pipeline"].process.nodes_activation = {
            "outliercount": True,
            "boldiqms": True,
            "carpetparcellation": True,
            "framewisedisplacement": True,
            "spikes": True,
            "gcor": True,
            "computedvars": True,
            "fwhmx": True,
            "qualityindex": True,
        }
        self.add_process(
            "nonsteadystatedetector",
            "mia_processes.bricks.preprocess.others.processing."
            "NonSteadyStateDetector",
        )
        self.add_process(
            "sanitize",
            "mia_processes.bricks.preprocess.others.processing." "Sanitize",
        )
        self.add_process(
            "tsnr", "mia_processes.bricks.preprocess.others.processing." "TSNR"
        )
        self.add_process(
            "mean",
            "mia_processes.bricks.preprocess.afni.processes." "TStatMean",
        )
        self.add_process(
            "automask",
            "mia_processes.bricks.preprocess.afni.processes." "Automask",
        )
        self.add_process(
            "volreg",
            "mia_processes.bricks.preprocess.afni.processes." "Volreg",
        )
        self.nodes["volreg"].process.twopass = True
        self.nodes["volreg"].process.interpolation = "Fourier"
        self.add_process(
            "bold_mni_align",
            "mia_processes.pipelines.preprocess.bold_mni_align."
            "Bold_mni_align",
        )
        self.nodes["bold_mni_align"].set_plug_value(
            "epi_mask", traits.Undefined
        )
        self.nodes["bold_mni_align"].process.nodes[
            "registration"
        ].set_plug_value("moving_image_masks", traits.Undefined)
        self.nodes["bold_mni_align"].process.nodes_activation = {
            "affineinitializer": True,
            "registration": True,
            "n4biasfieldcorrection": True,
            "applytransforms": True,
            "template": True,
            "template_mask": True,
            "seg_template": True,
        }
        self.nodes["bold_mni_align"].process.epi_mask = traits.Undefined
        self.add_process(
            "mriqc_func_report",
            "mia_processes.bricks.reports.reporting." "ReportFuncMriqc",
        )
        self.add_process(
            "plot_iqms", "mia_processes.bricks.reports.processes.BoldIQMsPlot"
        )

        # links
        self.export_parameter(
            "sanitize",
            "in_file",
            pipeline_parameter="func_file",
            is_optional=False,
        )
        self.add_link("func_file->nonsteadystatedetector.in_file")
        self.add_link("func_file->mriqc_func_report.func")
        self.add_link(
            "bold_iqms_pipeline.BoldQC_out_file->"
            "mriqc_func_report.IQMs_file"
        )
        self.add_link(
            "nonsteadystatedetector.n_volumes_to_discard->"
            "sanitize.n_volumes_to_discard"
        )
        self.add_link(
            "nonsteadystatedetector.n_volumes_to_discard->"
            "bold_iqms_pipeline.dummy_TRs"
        )
        self.add_link("sanitize.out_file->bold_iqms_pipeline.ras_epi")
        self.add_link("sanitize.out_file->volreg.in_file")
        self.add_link("tsnr.out_tsnr_file->bold_iqms_pipeline.epi_tsnr")
        self.add_link("mean.out_file->bold_mni_align.epi_mean")
        self.add_link("mean.out_file->automask.in_file")
        self.add_link("mean.out_file->bold_iqms_pipeline.epi_mean")
        self.add_link("automask.out_file->bold_mni_align.epi_mask")
        self.add_link("automask.out_file->bold_iqms_pipeline.brainmask")
        self.add_link("volreg.out_file->tsnr.in_file")
        self.add_link("volreg.out_file->mean.in_file")
        self.add_link("volreg.out_file->bold_iqms_pipeline.hmc_epi")
        self.add_link("volreg.oned_file->" "bold_iqms_pipeline.hmc_motion")
        self.add_link("bold_mni_align.epi_parc->bold_iqms_pipeline.epi_parc")
        self.add_link("bold_mni_align.epi_mni->mriqc_func_report.norm_func")
        self.add_link(
            "bold_iqms_pipeline.outliercount_file->"
            "plot_iqms.in_outliers_file"
        )
        self.add_link("bold_iqms_pipeline.fd_file->plot_iqms.in_fd_file")
        self.add_link("bold_iqms_pipeline.dvars_file->plot_iqms.in_dvars_file")
        self.add_link(
            "bold_iqms_pipeline.spike_file->" "plot_iqms.in_spikes_file"
        )
        self.add_link(
            "nonsteadystatedetector.n_volumes_to_discard->"
            "plot_iqms.drop_trs"
        )
        self.add_link("bold_iqms_pipeline.carpet_seg->plot_iqms.carpet_seg")
        self.add_link("volreg.out_file->plot_iqms.in_func")
        self.add_link("plot_iqms.out_file->mriqc_func_report.IQMs_plot")
        self.add_link("tsnr.out_stddev_file->mriqc_func_report.stddev_func")
        self.add_link("mean.out_file->mriqc_func_report.func_mean")
        self.add_link("automask.out_file->mriqc_func_report.brain_mask")
        self.export_parameter(
            "mriqc_func_report",
            "report",
            pipeline_parameter="func_report",
            is_optional=True,
        )

        # parameters order
        self.reorder_traits(("func_file", "func_report"))

        # nodes positions
        self.node_position = {
            "bold_iqms_pipeline": (842.4679035022216, 159.92981034666673),
            "nonsteadystatedetector": (-717.1107544384864, 113.59195862760868),
            "sanitize": (-234.51652984508183, 346.43391671594577),
            "tsnr": (-127.7391993962027, 60.050718720000134),
            "mean": (-173.11641600000007, -181.6966963199999),
            "automask": (144.1320191999998, -189.81009407999983),
            "volreg": (-527.8792115199997, -149.2405196799998),
            "bold_mni_align": (283.31445431594545, 116.08725895240536),
            "inputs": (-890.4021152762024, 569.7810895696218),
            "outputs": (1530.1459691893501, 147.22476292910784),
            "mriqc_func_report": (1549.0303706025636, 592.4476823062976),
            "plot_iqms": (1300.6566434431888, 640.3365297646346),
        }

        # nodes dimensions
        self.node_dimension = {
            "bold_iqms_pipeline": (364.546875, 810.0),
            "nonsteadystatedetector": (233.125, 75.0),
            "sanitize": (243.9375, 215.0),
            "tsnr": (239.265625, 215.0),
            "mean": (175.671875, 145.0),
            "automask": (217.53125, 285.0),
            "volreg": (251.46875, 425.0),
            "bold_mni_align": (459.453125, 670.0),
            "inputs": (106.00063723929541, 65.0),
            "outputs": (107.25, 110.0),
            "mriqc_func_report": (260.484375, 460.0),
            "plot_iqms": (192.265625, 355.0),
        }

        self.do_autoexport_nodes_parameters = False
