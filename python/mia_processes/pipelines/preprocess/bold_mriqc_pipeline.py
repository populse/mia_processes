# -*- coding: utf-8 -*-

from capsul.api import Pipeline
import traits.api as traits


class Bold_mriqc_pipeline(Pipeline):

    def pipeline_definition(self):
        # nodes
        self.add_process("bold_iqms_pipeline_1", "mia_processes.pipelines.reports.bold_iqms_pipeline.Bold_iqms_pipeline")
        self.nodes["bold_iqms_pipeline_1"].process.nodes_activation = {'outliercount': True, 'boldiqms': True, 'carpetparcellation': True, 'framewisedisplacement': True, 'spikes': True, 'gcor': True, 'computedvars': True, 'fwhmx': True, 'qualityindex': True}
        self.add_process("nonsteadystatedetector_1", "mia_processes.bricks.preprocess.others.processing.NonSteadyStateDetector")
        self.add_process("sanitize_1", "mia_processes.bricks.preprocess.others.processing.Sanitize")
        self.add_process("tsnr_1", "mia_processes.bricks.preprocess.others.processing.TSNR")
        self.add_process("mean_1", "mia_processes.bricks.preprocess.afni.processes.Mean")
        self.add_process("automask_1", "mia_processes.bricks.preprocess.afni.processes.Automask")
        self.add_process("bold_hmc_pipeline_1", "mia_processes.pipelines.preprocess.bold_hmc_pipeline.Bold_hmc_pipeline")
        self.nodes["bold_hmc_pipeline_1"].process.nodes_activation = {'droptrs': True, 'tshift': True, 'despike': True, 'deoblique': True, 'volreg': True}
        self.add_process("bold_mni_align_1", "mia_processes.pipelines.preprocess.bold_mni_align.Bold_mni_align")
        self.nodes["bold_mni_align_1"].process.nodes_activation = {'affineinitializer': True, 'registration': True, 'n4biasfieldcorrection': True, 'applytransforms': True, 'template': True, 'template_mask': True, 'seg_template': True}
        self.nodes["bold_mni_align_1"].process.epi_mask = traits.Undefined

        # links
        self.export_parameter("sanitize_1", "in_file", is_optional=False)
        self.add_link("in_file->nonsteadystatedetector_1.in_file")
        self.export_parameter("bold_mni_align_1", "template", is_optional=False)
        self.export_parameter("bold_mni_align_1", "template_res", is_optional=False)
        self.export_parameter("bold_mni_align_1", "seg_template_res", is_optional=False)
        self.export_parameter("bold_mni_align_1", "reg_transforms", is_optional=False)
        self.export_parameter("bold_mni_align_1", "reg_transform_parameters", is_optional=False)
        self.export_parameter("bold_mni_align_1", "reg_metric", is_optional=False)
        self.export_parameter("bold_mni_align_1", "reg_metric_weight", is_optional=False)
        self.export_parameter("bold_mni_align_1", "reg_shrink_factors", is_optional=False)
        self.export_parameter("bold_mni_align_1", "reg_smoothing_sigmas", is_optional=False)
        self.export_parameter("bold_mni_align_1", "reg_number_of_iterations", is_optional=True)
        self.export_parameter("bold_mni_align_1", "reg_radius_or_number_of_bins", is_optional=True)
        self.export_parameter("bold_mni_align_1", "reg_convergence_threshold", is_optional=True)
        self.export_parameter("bold_mni_align_1", "reg_convergence_window_size", is_optional=True)
        self.export_parameter("bold_mni_align_1", "reg_sampling_percentage", is_optional=True)
        self.export_parameter("bold_mni_align_1", "reg_sampling_strategy", is_optional=True)
        self.export_parameter("bold_mni_align_1", "reg_interpolation", is_optional=True)
        self.export_parameter("bold_iqms_pipeline_1", "outlier_fraction", "iqm_outlier_fraction", is_optional=True)
        self.export_parameter("bold_iqms_pipeline_1", "quality_index_automask", "iqm_quality_index_automask", is_optional=False)
        self.export_parameter("bold_iqms_pipeline_1", "fwhm_combine", "iqm_fwhm_combine", is_optional=True)
        self.export_parameter("bold_iqms_pipeline_1", "fwhm_detrend", "iqm_fwhm_detrend", is_optional=False)
        self.export_parameter("bold_iqms_pipeline_1", "dvars_remove_zero_variance", "iqm_dvars_remove_zero_variance", is_optional=True)
        self.export_parameter("bold_iqms_pipeline_1", "dvars_intensity_normalization", "iqm_dvars_intensity_normalization", is_optional=True)
        self.export_parameter("bold_iqms_pipeline_1", "fd_parameter_source", "iqm_fd_parameter_source", is_optional=False)
        self.export_parameter("bold_iqms_pipeline_1", "fd_radius", "iqm_fd_radius", is_optional=True)
        self.export_parameter("bold_iqms_pipeline_1", "fd_normalize", "iqm_fd_normalize", is_optional=True)
        self.export_parameter("bold_iqms_pipeline_1", "fd_thresh", "iqm_fd_thresh", is_optional=True)
        self.export_parameter("bold_iqms_pipeline_1", "spikes_no_zscore", "iqm_spikes_no_zscore", is_optional=True)
        self.export_parameter("bold_iqms_pipeline_1", "spikes_detrend", "iqm_spikes_detrend", is_optional=True)
        self.export_parameter("bold_iqms_pipeline_1", "spikes_spike_thresh", "iqm_spikes_spike_thresh", is_optional=True)
        self.export_parameter("bold_iqms_pipeline_1", "spikes_skip_frames", "iqm_spikes_skip_frames", is_optional=True)
        self.export_parameter("bold_hmc_pipeline_1", "start_idx", "hmc_start_idx", is_optional=False)
        self.export_parameter("bold_hmc_pipeline_1", "stop_idx", "hmc_stop_idx", is_optional=False)
        self.export_parameter("bold_hmc_pipeline_1", "despike", "hmc_despike", is_optional=False)
        self.export_parameter("bold_hmc_pipeline_1", "deoblique", "hmc_deoblique", is_optional=False)
        self.export_parameter("bold_iqms_pipeline_1", "carpet_seg", is_optional=False)
        self.export_parameter("bold_iqms_pipeline_1", "BoldQC_out_file", is_optional=False)
        self.add_link("nonsteadystatedetector_1.n_volumes_to_discard->sanitize_1.n_volumes_to_discard")
        self.add_link("nonsteadystatedetector_1.n_volumes_to_discard->bold_iqms_pipeline_1.dummy_TRs")
        self.add_link("sanitize_1.out_file->bold_iqms_pipeline_1.ras_epi")
        self.add_link("sanitize_1.out_file->bold_hmc_pipeline_1.in_file")
        self.export_parameter("tsnr_1", "out_tsnr_file", "tsnr_file", is_optional=False)
        self.add_link("tsnr_1.out_tsnr_file->bold_iqms_pipeline_1.epi_tsnr")
        self.export_parameter("tsnr_1", "out_stddev_file", "stddev_file", is_optional=False)
        self.add_link("mean_1.out_file->bold_mni_align_1.epi_mean")
        self.add_link("mean_1.out_file->bold_iqms_pipeline_1.epi_mean")
        self.export_parameter("mean_1", "out_file", "epi_mean_file", is_optional=False)
        self.add_link("mean_1.out_file->automask_1.in_file")
        self.add_link("automask_1.out_file->bold_iqms_pipeline_1.brainmask")
        self.add_link("automask_1.out_file->bold_mni_align_1.epi_mask")
        self.export_parameter("automask_1", "out_file", "brainmask_file", is_optional=False)
        self.add_link("bold_hmc_pipeline_1.out_file->bold_iqms_pipeline_1.hmc_epi")
        self.add_link("bold_hmc_pipeline_1.out_file->mean_1.in_file")
        self.add_link("bold_hmc_pipeline_1.out_file->tsnr_1.in_file")
        self.add_link("bold_hmc_pipeline_1.oned_file->bold_iqms_pipeline_1.hmc_motion")
        self.add_link("bold_mni_align_1.epi_parc->bold_iqms_pipeline_1.epi_parc")
        self.export_parameter("bold_mni_align_1", "epi_mni", "epi_mni_file", is_optional=False)
        self.export_parameter("bold_mni_align_1", "composite_transform", "_composite_transform", is_optional=False)
        self.export_parameter("bold_mni_align_1", "inverse_composite_transform", "_inverse_composite_transform", is_optional=False)
        self.export_parameter("bold_mni_align_1", "bias_image", "_bias_image", is_optional=False)

        # parameters order

        self.reorder_traits(("in_file", "tsnr_file", "stddev_file", "epi_mni_file", "epi_mean_file", "brainmask_file", "BoldQC_out_file", "carpet_seg", "_composite_transform", "_inverse_composite_transform", "_bias_image", "template", "template_res", "seg_template_res", "reg_transforms", "reg_transform_parameters", "reg_metric", "reg_metric_weight", "reg_shrink_factors", "reg_smoothing_sigmas", "reg_number_of_iterations", "reg_radius_or_number_of_bins", "reg_convergence_threshold", "reg_convergence_window_size", "reg_sampling_percentage", "reg_sampling_strategy", "reg_interpolation", "iqm_outlier_fraction", "iqm_quality_index_automask", "iqm_fwhm_combine", "iqm_fwhm_detrend", "iqm_dvars_remove_zero_variance", "iqm_dvars_intensity_normalization", "iqm_fd_parameter_source", "iqm_fd_radius", "iqm_fd_normalize", "iqm_fd_thresh", "iqm_spikes_no_zscore", "iqm_spikes_detrend", "iqm_spikes_spike_thresh", "iqm_spikes_skip_frames", "hmc_start_idx", "hmc_stop_idx", "hmc_despike", "hmc_deoblique"))

        # default and initial values
        self.template = 'MNI152NLin2009cAsym'
        self.seg_template_res = 1
        self.reg_transforms = ['Rigid', 'Affine', 'SyN']
        self.reg_transform_parameters = [(0.05,), (0.08,), (0.1, 3.0, 0.0)]
        self.reg_metric = ['Mattes', 'Mattes', 'CC']
        self.reg_metric_weight = [1, 1, 1]
        self.reg_shrink_factors = [[4, 2, 1], [8, 4, 2], [8, 4, 2]]
        self.reg_smoothing_sigmas = [[4.0, 2.0, 1.0], [4.0, 2.0, 0.0], [3.0, 2.0, 1.0]]
        self.reg_number_of_iterations = [[10000, 1000, 100], [500, 250, 100], [100, 30, 20]]
        self.reg_radius_or_number_of_bins = [56, 56, 4]
        self.reg_convergence_threshold = [1e-06, 1e-06, 1e-06]
        self.reg_convergence_window_size = [20, 20, 10]
        self.reg_sampling_percentage = [0.25, 0.25, 1.0]
        self.reg_sampling_strategy = ['Regular', 'Regular', 'None']
        self.reg_interpolation = 'LanczosWindowedSinc'
        self.iqm_fd_parameter_source = 'AFNI'

        # nodes positions
        self.node_position = {
            "bold_iqms_pipeline_1": (842.4679035022216, 159.92981034666673),
            "nonsteadystatedetector_1": (-812.0104243199997, 341.9442892799999),
            "sanitize_1": (-542.9404569599999, 287.1216230400001),
            "tsnr_1": (-130.70481408, 60.050718720000134),
            "mean_1": (-173.11641600000007, -181.6966963199999),
            "automask_1": (144.1320191999998, -189.81009407999983),
            "bold_hmc_pipeline_1": (-527.8792115199997, -149.2405196799998),
            "bold_mni_align_1": (224.00216063999977, 122.01848831999996),
            "inputs": (-893.3677299599997, 593.50600704),
            "outputs": (1257.3094182799998, 46.39386368000004),
        }

        # nodes dimensions
        self.node_dimension = {
            "bold_iqms_pipeline_1": (383.640625, 810.0),
            "nonsteadystatedetector_1": (226.4375, 75.0),
            "sanitize_1": (236.53125, 215.0),
            "tsnr_1": (230.671875, 215.0),
            "mean_1": (165.65625, 145.0),
            "automask_1": (165.65625, 145.0),
            "bold_hmc_pipeline_1": (210.625, 215.0),
            "bold_mni_align_1": (481.546875, 670.0),
            "inputs": (256.71875, 1265.0),
            "outputs": (229.78125, 390.0),
        }

        self.do_autoexport_nodes_parameters = False
