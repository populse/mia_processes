from capsul.api import Pipeline
import traits.api as traits


class Bold_hmc_pipeline(Pipeline):

    def pipeline_definition(self):
        # nodes
        self.add_process("droptrs", "mia_processes.bricks.preprocess.afni.processes.DropTRs")
        self.add_process("tshift", "mia_processes.bricks.preprocess.afni.processes.TShift")
        self.add_process("despike", "mia_processes.bricks.preprocess.afni.processes.Despike")
        self.add_process("deoblique", "mia_processes.bricks.preprocess.afni.processes.Deoblique")
        self.add_process("volreg", "mia_processes.bricks.preprocess.afni.processes.Volreg")
        self.nodes["volreg"].process.twopass = True

        # links
        self.export_parameter("droptrs", "in_file", is_optional=False)
        self.export_parameter("droptrs", "start_idx", is_optional=False)
        self.export_parameter("droptrs", "stop_idx", is_optional=False)
        self.export_parameter("despike", "despike", is_optional=False)
        self.export_parameter("deoblique", "deoblique", is_optional=False)
        self.add_link("droptrs.out_file->tshift.in_file")
        self.add_link("tshift.out_file->despike.in_file")
        self.add_link("despike.out_file->deoblique.in_file")
        self.add_link("deoblique.out_file->volreg.in_file")
        self.export_parameter("volreg", "out_file", is_optional=False)
        self.export_parameter("volreg", "oned_file", is_optional=False)

        # parameters order

        self.reorder_traits(("in_file", "start_idx", "stop_idx", "despike", "deoblique", "out_file", "oned_file"))

        # nodes positions
        self.node_position = {
            "droptrs": (-364.0, -61.0),
            "tshift": (-189.79999999999998, 113.80000000000004),
            "despike": (21.799999999999983, 294.59999999999997),
            "deoblique": (186.20000000000005, 439.9999999999999),
            "volreg": (353.6, 514.3999999999999),
            "inputs": (-575.73325, 294.9679999999999),
            "outputs": (633.9562500000001, 514.3999999999999),
        }

        # nodes dimensions
        self.node_dimension = {
            "droptrs": (165.65625, 215.0),
            "tshift": (207.046875, 215.0),
            "despike": (165.65625, 180.0),
            "deoblique": (151.03125, 110.0),
            "volreg": (176.625, 250.0),
            "inputs": (110.26875, 152.0),
            "outputs": (90.03125, 83.0),
        }

        self.do_autoexport_nodes_parameters = False
