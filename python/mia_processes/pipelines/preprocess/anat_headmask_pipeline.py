from capsul.api import Pipeline
import traits.api as traits


class Anat_headmask_pipeline(Pipeline):

    def pipeline_definition(self):
        # nodes
        self.add_process("denoise",
                         "mia_processes.bricks.preprocess."
                         "dipy.processes.Denoise")
        self.add_process("enhance",
                         "mia_processes.bricks.preprocess."
                         "others.processing.Enhance")
        self.add_process("gradient_threshold",
                         "mia_processes.bricks.preprocess."
                         "others.processing.GradientThreshold")

        # links
        self.export_parameter("denoise", "in_file_snr", "in_file",
                              is_optional=False)
        self.add_link("in_file->enhance.in_files")
        self.export_parameter("gradient_threshold", "seg_file",
                              is_optional=False)
        self.add_link("seg_file->denoise.seg_file")
        self.add_link("denoise.out_file->gradient_threshold.in_file")
        self.add_link("enhance.out_files->denoise.in_file")
        self.export_parameter("gradient_threshold", "out_file",
                              is_optional=False)

        # parameters order

        self.reorder_traits(("in_file", "seg_file", "out_file"))

        # nodes positions
        self.node_position = {
            "denoise": (1.0, -70.0),
            "inputs": (-346.415625, 4.0),
            "enhance": (-180.0, 182.0),
            "gradient_threshold": (216.0, 125.0),
            "outputs": (487.028125, 125.0),
        }

        # nodes dimensions
        self.node_dimension = {
            "denoise": (157.84375, 180.0),
            "inputs": (86.503125, 110.0),
            "enhance": (141.546875, 145.0),
            "gradient_threshold": (178.265625, 180.0),
            "outputs": (79.0625, 60.0),
        }

        self.do_autoexport_nodes_parameters = False
