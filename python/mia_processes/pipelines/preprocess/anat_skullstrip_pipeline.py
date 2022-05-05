from capsul.api import Pipeline
import traits.api as traits


class Anat_skullstrip_pipeline(Pipeline):

    def pipeline_definition(self):
        # nodes
        self.add_process("n4_bias_field_correction", "mia_processes.bricks.preprocess.ants.processes.N4BiasFieldCorrection")
        self.nodes["n4_bias_field_correction"].process.dimension = 3
        self.add_process("skull_stripping", "mia_processes.bricks.preprocess.afni.processes.SkullStripping")
        self.add_process("calc", "mia_processes.bricks.preprocess.afni.processes.Calc")
        self.nodes["calc"].process.expr = 'a*step(b)'
        self.nodes["calc"].process.out_prefix = 'ss_orig_'
        self.add_process("binarize", "mia_processes.bricks.preprocess.others.processing.Binarize")

        # links
        self.export_parameter("calc", "in_file_a", "in_file", is_optional=False)
        self.add_link("in_file->n4_bias_field_correction.in_file")
        self.export_parameter("calc", "expr", "expr_mask", is_optional=False)
        self.add_link("n4_bias_field_correction.out_file->skull_stripping.in_file")
        self.export_parameter("n4_bias_field_correction", "out_file", "bias_corrected", is_optional=False)
        self.export_parameter("n4_bias_field_correction", "bias_image", is_optional=False)
        self.add_link("skull_stripping.out_file->calc.in_file_b")
        self.export_parameter("calc", "out_file", is_optional=False)
        self.add_link("calc.out_file->binarize.in_files")
        self.export_parameter("binarize", "out_files", "out_mask", is_optional=False)

        # parameters order

        self.reorder_traits(("out_file", "out_mask", "bias_corrected", "bias_image", "in_file", "expr_mask"))

        # default and initial values
        self.expr_mask = 'a*step(b)'

        # nodes positions
        self.node_position = {
            "n4_bias_field_correction": (-656.0, 91.0),
            "skull_stripping": (-363.84000000000003, -41.879999999999995),
            "calc": (-188.39200000000002, -325.20799999999997),
            "binarize": (34.55199999999999, -79.86399999999996),
            "outputs": (277.71875, 29.0),
            "inputs": (-825.5356250000002, -74.59999999999998),
        }

        # nodes dimensions
        self.node_dimension = {
            "n4_bias_field_correction": (209.78125, 145.0),
            "skull_stripping": (165.65625, 145.0),
            "calc": (165.65625, 250.0),
            "binarize": (164.59375, 180.0),
            "outputs": (127.75, 129.0),
            "inputs": (116.346875, 83.0),
        }

        self.do_autoexport_nodes_parameters = False