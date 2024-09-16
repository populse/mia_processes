# -*- coding: utf-8 -*-
"""The matlab script preprocess library of the mia_processes package.

The purpose of this module is to add bricks that wrap MATLAB
scripts using Nipype.
Thoses bricks requires an exclusive Matlab license.

The matlab script should be added in the folder "matlab_wrap/sripts"
"""

import os
import re
import tempfile

import nibabel as nb
import nipype.interfaces.matlab as Matlab
from nipype.interfaces.base import File

# populse_mia import
from populse_mia.software_properties import Config
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA
from scipy.io import savemat
from traits.api import Int

from mia_processes.utils import checkFileExt

EXT = {"NIFTI_GZ": "nii.gz", "NIFTI": "nii"}


class ComputeBrainVolume(ProcessMIA):
    """
    *Compute brain volume*

    Please, see the complete documentation for the `ComputeBrainVolume brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/matlab_wrap/ComputeBrainVolume.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        "self.requirement" attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ComputeBrainVolume, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["nipype", "matlab"]
        # Matlab standalone is not working to launch Matlab script
        # TODO: add a way to check that matlab with
        # license (not standalone) is used

        # Mandatory inputs description
        in_file_desc = (
            "Input file (a pathlike object"
            "string representing an existing file)"
        )
        # Outputs description
        volume_desc = "Brain volume (int)"

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Outputs traits
        self.add_trait("volume", Int(output=True, desc=volume_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None, iteration=False):
        """Dedicated to the initialisation step of the brick.
        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key "notInDb" of the
        self.outputs dictionary. To work properly this method must return
        self.make_initResult() object.
        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(ComputeBrainVolume, self).list_outputs()
        config = Config()
        use_matlab = config.get_use_matlab()

        if not use_matlab:
            # This is a patch pending modification of the requirements check
            # at initialisation time. Currently, initialisation will fail,
            # but the user will not be informed why
            # (unless he looks at the stdout ....)
            print(
                "\nThis brick requires Matlab (with license). "
                "Please modify the configuration this way ...\n"
            )
            return

        if self.in_file:
            valid_ext, in_ext, file_name = checkFileExt(self.in_file, EXT)
            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return

            # Outputs definition and tags inheritance (optional)
            self.outputs["volume"] = 0

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ComputeBrainVolume, self).run_process_mia()

        # Matlab path configuration
        config = Config()
        use_matlab = config.get_use_matlab()
        if use_matlab:
            matlab_path = config.get_matlab_path()
            if "MATLABCMD" not in os.environ:
                os.environ["MATLABCMD"] = matlab_path
        tmp_dir = tempfile.TemporaryDirectory()

        # Specif pre-processing for this script
        valid_ext, in_ext, file_name = checkFileExt(self.in_file, EXT)
        image_mat = os.path.join(tmp_dir.name, file_name + ".mat")
        data = nb.load(self.in_file).get_fdata()
        savemat(image_mat, {"data": data}, do_compression=False)

        # Matlab script configuration
        # Get matlab script
        matlab_script = os.path.join(
            os.path.realpath(os.path.dirname(__file__)),
            "scripts",
            "compute_brainvolume.m",
        )
        # Create a new script with the good inputs / outputs
        with open(matlab_script, encoding="utf-8") as script_file:
            script_content = script_file.read()

        tmp_script_name = "tmp_script"
        tmp_script = os.path.join(tmp_dir.name, tmp_script_name + ".m")

        with open(tmp_script, "w", encoding="utf-8") as script_file:
            script_file.write(
                script_content.replace("input_image.mat", image_mat)
            )

        # Change working directory and launch script using matlab
        os.chdir(tmp_dir.name)
        Matlab.MatlabCommand.set_default_matlab_cmd(
            Matlab.get_matlab_command()
        )
        mlab = Matlab.MatlabCommand(script=tmp_script_name)
        result = mlab.run()

        # Get outputs using result.runtime.stdout
        expr_reg = r"total\s*=\s*(\d+)"
        correspondance = re.search(expr_reg, result.runtime.stdout)
        self.volume = int(correspondance.group(1))

        # Clean tmp directory
        tmp_dir.cleanup()
