# -*- coding: utf-8 -*-

from capsul.api import Pipeline
import traits.api as traits


class Bold_iqms_pipeline(Pipeline):

    def pipeline_definition(self):
        # nodes
        self.add_process("outliercount", "mia_processes.bricks.reports.processes.OutlierCount")
        self.add_process("boldiqms", "mia_processes.bricks.reports.processes.BoldIQMs")
        self.add_process("carpetparcellation", "mia_processes.bricks.reports.processes.CarpetParcellation")
        self.add_process("framewisedisplacement", "mia_processes.bricks.reports.processes.FramewiseDisplacement")
        self.nodes["framewisedisplacement"].process.parameter_source = 'AFNI'
        self.add_process("spikes", "mia_processes.bricks.reports.processes.Spikes")
        self.add_process("gcor", "mia_processes.bricks.reports.processes.GCOR")
        self.add_process("computedvars", "mia_processes.bricks.reports.processes.ComputeDVARS")
        self.add_process("fwhmx", "mia_processes.bricks.reports.processes.FWHMx")
        self.add_process("qualityindex", "mia_processes.bricks.reports.processes.QualityIndex")

        # links
        self.export_parameter("fwhmx", "in_file", "epi_mean", is_optional=False)
        self.add_link("epi_mean->boldiqms.in_epi")
        self.export_parameter("boldiqms", "in_hmc", "hmc_epi", is_optional=True)
        self.add_link("hmc_epi->computedvars.in_file")
        self.add_link("hmc_epi->qualityindex.in_file")
        self.add_link("hmc_epi->outliercount.in_file")
        self.add_link("hmc_epi->gcor.in_file")
        self.export_parameter("boldiqms", "in_tsnr", "epi_tsnr", is_optional=True)
        self.export_parameter("gcor", "mask_file", "brainmask", is_optional=True)
        self.add_link("brainmask->fwhmx.mask_file")
        self.add_link("brainmask->carpetparcellation.brainmask")
        self.add_link("brainmask->computedvars.in_mask")
        self.add_link("brainmask->boldiqms.in_mask")
        self.add_link("brainmask->outliercount.mask_file")
        self.export_parameter("carpetparcellation", "segmentation", "epi_parc", is_optional=False)
        self.export_parameter("framewisedisplacement", "in_file", "hmc_motion", is_optional=False)
        self.export_parameter("spikes", "in_file", "ras_epi", is_optional=False)
        self.export_parameter("outliercount", "fraction", "outlier_fraction", is_optional=True)
        self.export_parameter("qualityindex", "automask", "quality_index_automask", is_optional=False)
        self.export_parameter("fwhmx", "combine", "fwhm_combine", is_optional=True)
        self.export_parameter("fwhmx", "detrend", "fwhm_detrend", is_optional=False)
        self.export_parameter("computedvars", "remove_zero_variance", "dvars_remove_zero_variance", is_optional=True)
        self.export_parameter("computedvars", "intensity_normalization", "dvars_intensity_normalization", is_optional=True)
        self.export_parameter("framewisedisplacement", "parameter_source", "fd_parameter_source", is_optional=False)
        self.export_parameter("framewisedisplacement", "radius", "fd_radius", is_optional=True)
        self.export_parameter("framewisedisplacement", "normalize", "fd_normalize", is_optional=True)
        self.export_parameter("spikes", "no_zscore", "spikes_no_zscore", is_optional=True)
        self.export_parameter("spikes", "detrend", "spikes_detrend", is_optional=True)
        self.export_parameter("spikes", "spike_thresh", "spikes_spike_thresh", is_optional=True)
        self.export_parameter("spikes", "skip_frames", "spikes_skip_frames", is_optional=True)
        self.add_link("outliercount.out_file->boldiqms.in_outliers_file")
        self.export_parameter("boldiqms", "out_file", "BoldQC_out_file", is_optional=False)
        self.export_parameter("carpetparcellation", "out_file", "carpet_seg", is_optional=False)
        self.add_link("framewisedisplacement.out_file->boldiqms.in_fd_file")
        self.add_link("spikes.out_file->boldiqms.in_spikes_file")
        self.add_link("gcor.out->boldiqms.in_gcor")
        self.add_link("computedvars.out_file->boldiqms.in_dvars_file")
        self.add_link("fwhmx.out_file->boldiqms.in_fwhm_file")
        self.add_link("qualityindex.out_file->boldiqms.in_QI_file")

        # parameters order

        self.reorder_traits(("epi_mean", "hmc_epi", "epi_tsnr", "brainmask", "epi_parc", "hmc_motion", "ras_epi", "carpet_seg", "BoldQC_out_file", "outlier_fraction", "quality_index_automask", "fwhm_combine", "fwhm_detrend", "dvars_remove_zero_variance", "dvars_intensity_normalization", "fd_parameter_source", "fd_radius", "fd_normalize", "spikes_no_zscore", "spikes_detrend", "spikes_spike_thresh", "spikes_skip_frames"))

        # default and initial values
        self.fd_parameter_source = 'AFNI'

        # nodes positions
        self.node_position = {
            "outliercount": (-274.30400000000003, -168.95199999999994),
            "boldiqms": (151.84, -78.44799999999998),
            "carpetparcellation": (-80.77920000000003, -260.15679999999986),
            "inputs": (-714.5310000000001, -105.81600000000006),
            "outputs": (412.246875, -94.0),
            "framewisedisplacement": (-431.24, 342.48),
            "spikes": (-166.744, 452.03999999999996),
            "gcor": (32.63839999999999, 592.8847999999999),
            "computedvars": (-161.52800000000022, 183.92799999999983),
            "fwhmx": (-357.17600000000004, 38.98400000000001),
            "qualityindex": (-85.24000000000002, -28.67999999999995),
        }

        # nodes dimensions
        self.node_dimension = {
            "boldiqms_1": (185.15625, 495.0),
            "spikes_1": (169.953125, 250.0),
            "qualityindex_1": (152.234375, 145.0),
            "outliercount_1": (152.234375, 180.0),
            "gcor_1": (119.765625, 110.0),
            "framewisedisplacement_1": (239.28125, 215.0),
            "fwhmx_1": (152.234375, 215.0),
            "computedvars_1": (243.0625, 215.0),
            "carpetparcellation_1": (180.9375, 145.0),
            "outliercount": (152.234375, 180.0),
            "boldiqms": (185.15625, 495.0),
            "carpetparcellation": (179.5, 145.0),
            "inputs": (238.03125, 497.0),
            "outputs": (134.109375, 83.0),
            "framewisedisplacement": (206.359375, 215.0),
            "spikes": (169.953125, 250.0),
            "gcor": (119.765625, 110.0),
            "computedvars": (243.0625, 215.0),
            "fwhmx": (152.234375, 215.0),
            "qualityindex": (152.234375, 145.0),
        }

        self.do_autoexport_nodes_parameters = False