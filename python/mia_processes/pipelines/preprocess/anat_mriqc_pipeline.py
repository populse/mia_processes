from capsul.api import Pipeline
import traits.api as traits


class Anat_mriqc_pipeline(Pipeline):

    def pipeline_definition(self):
        # nodes
        self.add_process("conformimage", "mia_processes.bricks.preprocess.others.processing.ConformImage")
        self.add_process("anat_skullstrip_pipeline", "mia_processes.pipelines.preprocess.anat_skullstrip_pipeline.Anat_skullstrip_pipeline")
        self.nodes["anat_skullstrip_pipeline"].process.nodes_activation = {'n4_bias_field_correction': True, 'skull_stripping': True, 'calc': True, 'binarize': True}
        self.add_process("segment", "mia_processes.bricks.preprocess.fsl.processes.Segment")
        self.add_process("anat_spatial_norm", "mia_processes.pipelines.preprocess.anat_spatial_norm.Anat_spatial_norm")
        self.nodes["anat_spatial_norm"].process.nodes_activation = {'mask_moving_image': True, 'mask_fixed_image': True, 'affine_initializer': True, 'registration': True, 'template_mask': True, 'template': True}
        self.add_process("anat_headmask_pipeline", "mia_processes.pipelines.preprocess.anat_headmask_pipeline.Anat_headmask_pipeline")
        self.nodes["anat_headmask_pipeline"].process.nodes_activation = {'denoise': True, 'enhance': True, 'gradient_threshold': True}
        self.add_process("anat_airmask_pipeline", "mia_processes.pipelines.preprocess.anat_airmask_pipeline.Anat_airmask_pipeline")
        self.nodes["anat_airmask_pipeline"].process.nodes_activation = {'rotationmask': True, 'applytransforms': True, 'artifactmask': True, 'template': True}
        self.add_process("anat_mni_tpms_pipeline", "mia_processes.pipelines.preprocess.anat_mni_tpms_pipeline.Anat_mni_tpms_pipeline")
        self.nodes["anat_mni_tpms_pipeline"].process.nodes_activation = {'template_CSF': True, 'applytransforms_CSF': True, 'template_GM': True, 'template_WM': True, 'applytransforms_WM': True, 'applytransforms_GM': True, 'files_to_list': True}
        self.add_process("anatiqms", "mia_processes.bricks.reports.processes.AnatIQMs")
        self.add_process("harmonize", "mia_processes.bricks.preprocess.others.processing.Harmonize")
        self.add_process("fwhmx", "mia_processes.bricks.preprocess.afni.processes.FWHMx")
        self.add_process("list_to_file", "mia_processes.bricks.tools.tools.List_To_File")

        # links
        self.export_parameter("conformimage", "in_file", "anat_file", is_optional=False)
        self.export_parameter("anat_skullstrip_pipeline", "expr_mask", "ss_expr_mask", is_optional=False)
        self.export_parameter("anat_spatial_norm", "template", is_optional=False)
        self.add_link("template->anat_airmask_pipeline.template")
        self.add_link("template->anat_mni_tpms_pipeline.template")
        self.export_parameter("anat_spatial_norm", "template_res", "reg_template_res", is_optional=False)
        self.export_parameter("anat_airmask_pipeline", "tpl_res", "mask_template_res", is_optional=False)
        self.add_link("mask_template_res->anat_mni_tpms_pipeline.template_res")
        self.export_parameter("anat_spatial_norm", "reg_transforms", is_optional=False)
        self.export_parameter("anat_spatial_norm", "reg_transform_parameters", is_optional=False)
        self.export_parameter("anat_spatial_norm", "reg_metric", is_optional=False)
        self.export_parameter("anat_spatial_norm", "reg_metric_weight", is_optional=False)
        self.export_parameter("anat_spatial_norm", "reg_shrink_factors", is_optional=False)
        self.export_parameter("anat_spatial_norm", "reg_smoothing_sigmas", is_optional=False)
        self.export_parameter("anat_spatial_norm", "reg_number_of_iterations", is_optional=True)
        self.export_parameter("anat_spatial_norm", "reg_radius_or_number_of_bins", is_optional=True)
        self.export_parameter("anat_spatial_norm", "reg_convergence_threshold", is_optional=True)
        self.export_parameter("anat_spatial_norm", "reg_convergence_window_size", is_optional=True)
        self.export_parameter("anat_spatial_norm", "reg_sampling_percentage", is_optional=True)
        self.export_parameter("anat_spatial_norm", "reg_sampling_strategy", is_optional=True)
        self.export_parameter("anat_spatial_norm", "reg_interpolation", is_optional=True)
        self.export_parameter("fwhmx", "detrend", "fwhm_detrend", is_optional=False)
        self.export_parameter("conformimage", "suffix", "conform_suffix", is_optional=True)
        self.export_parameter("conformimage", "prefix", "conform_prefix", is_optional=True)
        self.export_parameter("segment", "segments", "reg_segments", is_optional=True)
        self.export_parameter("segment", "output_type", "reg_output_type", is_optional=True)
        self.export_parameter("anat_mni_tpms_pipeline", "template_suffix", "tpms_template_suffix", is_optional=True)
        self.export_parameter("harmonize", "erodemask", "harmonize_erodemask", is_optional=True)
        self.export_parameter("harmonize", "suffix", "harmonize_suffix", is_optional=True)
        self.export_parameter("harmonize", "prefix", "harmonize_prefix", is_optional=True)
        self.export_parameter("fwhmx", "combine", "fwhm_combine", is_optional=True)
        self.export_parameter("fwhmx", "out_prefix", "fwhm_out_prefix", is_optional=True)
        self.export_parameter("list_to_file", "index_filter", "wm_index_filter", is_optional=True)
        self.add_link("conformimage.out_file->anatiqms.in_ras")
        self.add_link("conformimage.out_file->anat_airmask_pipeline.in_file")
        self.add_link("conformimage.out_file->fwhmx.in_file")
        self.add_link("conformimage.out_file->anat_skullstrip_pipeline.in_file")
        self.add_link("conformimage.out_file->anat_mni_tpms_pipeline.in_ras")
        self.add_link("anat_skullstrip_pipeline.out_file->segment.in_file")
        self.add_link("anat_skullstrip_pipeline.out_mask->anat_airmask_pipeline.in_mask")
        self.add_link("anat_skullstrip_pipeline.out_mask->fwhmx.mask_file")
        self.add_link("anat_skullstrip_pipeline.out_mask->anat_spatial_norm.moving_mask")
        self.add_link("anat_skullstrip_pipeline.bias_corrected->anat_headmask_pipeline.in_file")
        self.add_link("anat_skullstrip_pipeline.bias_corrected->anat_spatial_norm.moving_image")
        self.add_link("anat_skullstrip_pipeline.bias_corrected->harmonize.in_file")
        self.add_link("anat_skullstrip_pipeline.bias_image->anatiqms.in_inu")
        self.add_link("segment.tissue_class_map->anat_headmask_pipeline.seg_file")
        self.add_link("segment.tissue_class_map->anatiqms.segmentation")
        self.add_link("segment.partial_volume_files->list_to_file.file_list")
        self.add_link("segment.partial_volume_files->anatiqms.pvms")
        self.export_parameter("anat_spatial_norm", "composite_transform", "_composite_transform", is_optional=False)
        self.add_link("anat_spatial_norm.inverse_composite_transform->anat_airmask_pipeline.inverse_composite_transform")
        self.add_link("anat_spatial_norm.inverse_composite_transform->anat_mni_tpms_pipeline.inverse_composite_transform")
        self.export_parameter("anat_spatial_norm", "warped_image", "_warped_image", is_optional=False)
        self.add_link("anat_headmask_pipeline.out_file->anat_airmask_pipeline.head_mask")
        self.add_link("anat_headmask_pipeline.out_file->anatiqms.headmask")
        self.add_link("anat_airmask_pipeline.out_hat_mask->anatiqms.hatmask")
        self.add_link("anat_airmask_pipeline.out_art_mask->anatiqms.artmask")
        self.add_link("anat_airmask_pipeline.out_air_mask->anatiqms.airmask")
        self.add_link("anat_airmask_pipeline.out_rot_mask->anatiqms.rotmask")
        self.add_link("anat_mni_tpms_pipeline.mni_tpms->anatiqms.mni_tpms")
        self.export_parameter("anatiqms", "out_file", "AnatQC_out_file", is_optional=False)
        self.add_link("harmonize.out_file->anatiqms.in_noinu")
        self.add_link("fwhmx.out_file->anatiqms.in_fwhm")
        self.add_link("list_to_file.file->harmonize.wm_mask")

        # parameters order

        self.reorder_traits(("anat_file", "ss_expr_mask", "template", "reg_template_res", "mask_template_res", "AnatQC_out_file", "_composite_transform", "_warped_image", "reg_transforms", "reg_transform_parameters", "reg_metric", "reg_metric_weight", "reg_shrink_factors", "reg_smoothing_sigmas", "reg_number_of_iterations", "reg_radius_or_number_of_bins", "reg_convergence_threshold", "reg_convergence_window_size", "reg_sampling_percentage", "reg_sampling_strategy", "reg_interpolation", "fwhm_detrend", "conform_suffix", "conform_prefix", "reg_segments", "reg_output_type", "tpms_template_suffix", "harmonize_erodemask", "harmonize_suffix", "harmonize_prefix", "fwhm_combine", "fwhm_out_prefix", "wm_index_filter"))

        # default and initial values
        self.ss_expr_mask = 'a*step(b)'
        self.template = 'MNI152NLin2009cAsym'
        self.mask_template_res = 1
        self.reg_transforms = ['Rigid', 'Affine', 'SyN']
        self.reg_transform_parameters = [(0.01,), (0.08,), (0.1, 3.0, 0.0)]
        self.reg_metric = ['Mattes', 'Mattes', 'Mattes']
        self.reg_metric_weight = [1, 1, 1]
        self.reg_shrink_factors = [[4], [4, 2, 1], [2, 1]]
        self.reg_smoothing_sigmas = [[4.0], [4.0, 2.0, 0.0], [1.0, 0.0]]
        self.reg_number_of_iterations = [[1000], [500, 250, 100], [50, 20]]
        self.reg_radius_or_number_of_bins = [32, 32, 56]
        self.reg_convergence_threshold = [1e-06, 1e-06, 1e-06]
        self.reg_convergence_window_size = [20, 20, 10]
        self.reg_sampling_percentage = [0.15, 0.15, 0.25]
        self.reg_sampling_strategy = ['Random', 'Regular', 'Regular']
        self.reg_interpolation = 'LanczosWindowedSinc'
        self.conform_suffix = ' '
        self.conform_prefix = ' '
        self.tpms_template_suffix = 'probseg'
        self.harmonize_prefix = ' '
        self.wm_index_filter = [2]

        # nodes positions
        self.node_position = {
            "conformimage": (-590.006976, 216.9358719999999),
            "inputs": (-951.4481057999996, 184.46495488000005),
            "anat_skullstrip_pipeline": (-405.8763555555557, 40.05591999999988),
            "segment": (-82.8528, -26.389119999999984),
            "anat_spatial_norm": (-147.45440000000005, 361.46943999999974),
            "anat_headmask_pipeline": (129.00223999999986, 204.37087999999997),
            "anat_airmask_pipeline": (492.1126399999997, 64.621376),
            "anat_mni_tpms_pipeline": (-59.65056000000003, 1168.9574400000004),
            "anatiqms": (538.1351423999997, 483.30971135999994),
            "harmonize": (-396.8870399999999, 1079.1014399999997),
            "fwhmx": (395.0622719999998, 1109.8736639999993),
            "list_to_file": (-612.9561599999998, 1084.7416319999998),
            "outputs": (815.9185139999993, 607.3239679999998),
        }

        # nodes dimensions
        self.node_dimension = {
            "conformimage": (150.40625, 145.0),
            "inputs": (227.765625, 1090.0),
            "anat_skullstrip_pipeline": (249.859375, 180.0),
            "segment": (250.234375, 145.0),
            "anat_spatial_norm": (481.546875, 635.0),
            "anat_headmask_pipeline": (248.375, 110.0),
            "anat_airmask_pipeline": (332.34375, 250.0),
            "anat_mni_tpms_pipeline": (302.875, 215.0),
            "anatiqms": (179.5, 460.0),
            "harmonize": (159.734375, 215.0),
            "fwhmx": (155.234375, 215.0),
            "list_to_file": (131.84375, 110.0),
            "outputs": (174.03125, 145.0),
        }

        self.do_autoexport_nodes_parameters = False
