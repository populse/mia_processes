from capsul.api import Pipeline
import os
import traits.api as traits


class Anat_airmask(Pipeline):
    """
    * | Compute rotation mask, air mask, artifact mask and hat mask for
      | structural data.

    following step 1 from
    (`[Mortamet2009] <https://onlinelibrary.wiley.com/doi/10.1002/mrm.21992>`_).

    """

    def pipeline_definition(self):
        # nodes
        self.add_process("rotationmask", "mia_processes.bricks.preprocess."
                         "others.processing.RotationMask")
        self.add_process("applytransforms", "mia_processes.bricks.preprocess."
                         "ants.processes.ApplyTransforms")
        self.add_process("artifactmask", "mia_processes.bricks.preprocess."
                         "others.processing.ArtifactMask")
        self.add_process("template", "mia_processes.bricks.preprocess."
                         "others.processing.TemplateFromTemplateFlow")
        self.nodes["template"].process.in_template = 'MNI152NLin2009cAsym'
        self.nodes["template"].process.resolution = 1
        self.nodes["template"].process.suffix = 'mask'
        self.nodes["template"].process.desc = 'head'
        self.nodes["applytransforms"].process.interpolation = 'MultiLabel'
        self.nodes["applytransforms"].process.float = True
        self.nodes["applytransforms"].process.num_threads = int(
            os.getenv("OMP_NUM_THREADS", os.cpu_count()))

        # links
        self.export_parameter("artifactmask", "in_file", is_optional=False)
        self.add_link("in_file->rotationmask.in_file")
        self.export_parameter("applytransforms", "reference_image",
                              "in_mask", is_optional=False)
        self.export_parameter("applytransforms", "transforms",
                              "inverse_composite_transform", is_optional=False)
        self.export_parameter("artifactmask", "head_mask", is_optional=False)
        self.export_parameter("rotationmask", "out_file",
                              "out_rot_mask", is_optional=False)
        self.add_link("rotationmask.out_file->"
                      "artifactmask.rot_mask")
        self.add_link("applytransforms.output_image->"
                      "artifactmask.nasion_post_mask")
        self.export_parameter("artifactmask", "out_hat_mask",
                              is_optional=False)
        self.export_parameter("artifactmask", "out_art_mask",
                              is_optional=False)
        self.export_parameter("artifactmask", "out_air_mask",
                              is_optional=False)
        self.add_link("template.template->applytransforms.input_image")

        # parameters order
        self.reorder_traits(("in_file", "in_mask",
                             "inverse_composite_transform", "head_mask",
                             "out_hat_mask", "out_art_mask", "out_air_mask"))

        # default and initial values
        self.template = 'MNI152NLin2009cAsym'
        self.tpl_res = 1

        # nodes positions
        self.node_position = {
            "rotationmask": (-477.0, -307.0),
            "applytransforms": (-407.33600000000007, -12.78400000000002),
            "artifactmask": (-69.0, -133.0),
            "template": (-694.6239999999999, 126.00000000000011),
            "inputs": (-1007.2956249999999, -146.15999999999985),
            "outputs": (326.85625, -133.0),
        }

        # nodes dimensions
        self.node_dimension = {
            "rotationmask": (141.359375, 145.0),
            "applytransforms": (248.328125, 215.0),
            "artifactmask": (258.59375, 215.0),
            "template": (183.640625, 215.0),
            "inputs": (223.78125, 250.0),
            "outputs": (123.5625, 106.0),
        }

        self.do_autoexport_nodes_parameters = False
