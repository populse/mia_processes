# -*- coding: utf-8 -*-

"""The volBrain (https://github.com/volBrain/AssemblyNet)
preprocess library of the mia_processes package.

The purpose of this module is to launch volbrain docker in mia_processes

:Contains:
    :Class:
        - Assemblynet

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# Other import
import os
import subprocess

from nipype.interfaces.base import File

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

from mia_processes.utils import checkFileExt

EXT = {"NIFTI_GZ": "nii.gz", "NIFTI": "nii"}


class AssemblyNet(ProcessMIA):
    """
    *3D Whole Brain MRI Segmentation using AssemblyNet (volBrain)*

    Please, see the complete documentation for the `AssemblyNet brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/volbrain/AssemblyNet.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(AssemblyNet, self).__init__()

        # Third party softwares required for the execution of the brick
        # TODO: add reuqirement (docker needed)

        # Mandatory inputs description
        in_file_desc = (
            "Input T1 file (a pathlike object"
            "string representing an existing file)"
        )

        # Outputs description
        native_t1_desc = (
            "Filtered and normalized T1 image in native space "
            "(a pathlike object or string representing a file)"
        )
        native_structures_desc = (
            "Structure segmentation in native space"
            "(a pathlike object or string representing a file)"
        )
        native_mask_desc = (
            "Intracranial Cavity mask image in native space"
            "(a pathlike object or string representing a file)"
        )
        native_lobes_desc = (
            "Lobe segmentation in native space"
            "(a pathlike object or string representing a file)"
        )
        native_macrosctrutures_desc = (
            "Macrostructure segmentation in native space"
            "(a pathlike object or string representing a file)"
        )
        native_tissues_desc = (
            "Intracranial Cavity mask image in native space"
            "(a pathlike object or string representing a file)"
        )
        mni_t1_desc = (
            "Filtered and normalized T1 image in MNI space"
            "(a pathlike object or string representing a file)"
        )
        mni_mask_desc = (
            "IIntracranial Cavity mask image in MNI space"
            "(a pathlike object or string representing a file)"
        )
        mni_tissues_desc = (
            "Tissues segmentation in MNI space"
            "(a pathlike object or string representing a file)"
        )
        mni_macrostructures_desc = (
            "Macrostructure segmentation in MNI space"
            "(a pathlike object or string representing a file)"
        )
        mni_lobes_desc = (
            "Lobe segmentation in MNI space"
            "(a pathlike object or string representing a file)"
        )
        mni_structures_desc = (
            "Structure segmentation in MNI space"
            "(a pathlike object or string representing a file)"
        )
        report_pdf_desc = (
            "PDF format volumetry report"
            "(a pathlike object or string representing a file)"
        )
        report_csv_desc = (
            "CSV format volumetry report."
            "(a pathlike object or string representing a file)"
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits

        # Outputs traits
        self.add_trait(
            "native_t1", File(output=True, optional=True, desc=native_t1_desc)
        )
        self.add_trait(
            "native_structures",
            File(output=True, optional=True, desc=native_structures_desc),
        )
        self.add_trait(
            "native_mask",
            File(output=True, optional=True, desc=native_mask_desc),
        )
        self.add_trait(
            "native_lobes",
            File(output=True, optional=True, desc=native_lobes_desc),
        )
        self.add_trait(
            "native_macrostructures",
            File(output=True, optional=True, desc=native_macrosctrutures_desc),
        )
        self.add_trait(
            "native_tissues",
            File(output=True, optional=True, desc=native_tissues_desc),
        )
        self.add_trait(
            "mni_t1", File(output=True, optional=True, desc=mni_t1_desc)
        )
        self.add_trait(
            "mni_mask", File(output=True, optional=True, desc=mni_mask_desc)
        )
        self.add_trait(
            "mni_tissues",
            File(output=True, optional=True, desc=mni_tissues_desc),
        )
        self.add_trait(
            "mni_macrostructures",
            File(output=True, optional=True, desc=mni_macrostructures_desc),
        )
        self.add_trait(
            "mni_lobes", File(output=True, optional=True, desc=mni_lobes_desc)
        )
        self.add_trait(
            "mni_structures",
            File(output=True, optional=True, desc=mni_structures_desc),
        )
        self.add_trait(
            "report_pdf",
            File(output=True, optional=True, desc=report_pdf_desc),
        )
        self.add_trait(
            "report_csv",
            File(output=True, optional=True, desc=report_csv_desc),
        )

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.
        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key 'notInDb' of the
        self.outputs dictionary. To work properly this method must return
        self.make_initResult() object.
        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(AssemblyNet, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            valid_ext, in_ext, file_name = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print("\nThe input image format is not recognized...!")
                return

            if self.output_directory:
                self.outputs["native_t1"] = os.path.join(
                    self.output_directory, "native_t1_" + file_name + ".nii.gz"
                )
                self.outputs["native_structures"] = os.path.join(
                    self.output_directory,
                    "native_structures_" + file_name + ".nii.gz",
                )
                self.outputs["native_mask"] = os.path.join(
                    self.output_directory,
                    "native_mask_" + file_name + ".nii.gz",
                )
                self.outputs["native_lobes"] = os.path.join(
                    self.output_directory,
                    "native_lobes_" + file_name + ".nii.gz",
                )
                self.outputs["native_macrostructures"] = os.path.join(
                    self.output_directory,
                    "native_macrostructures_" + file_name + ".nii.gz",
                )
                self.outputs["native_tissues"] = os.path.join(
                    self.output_directory,
                    "native_tissues_" + file_name + ".nii.gz",
                )
                self.outputs["mni_t1"] = os.path.join(
                    self.output_directory, "mni_t1_" + file_name + ".nii.gz"
                )
                self.outputs["mni_structures"] = os.path.join(
                    self.output_directory,
                    "mni_structures_" + file_name + ".nii.gz",
                )
                self.outputs["mni_mask"] = os.path.join(
                    self.output_directory, "mni_mask_" + file_name + ".nii.gz"
                )
                self.outputs["mni_lobes"] = os.path.join(
                    self.output_directory, "mni_lobes_" + file_name + ".nii.gz"
                )
                self.outputs["mni_macrostructures"] = os.path.join(
                    self.output_directory,
                    "mni_macrostructures_" + file_name + ".nii.gz",
                )
                self.outputs["mni_tissues"] = os.path.join(
                    self.output_directory,
                    "mni_tissues_" + file_name + ".nii.gz",
                )
                self.outputs["report_csv"] = os.path.join(
                    self.output_directory, "report_" + file_name + ".csv"
                )
                self.outputs["report_pdf"] = os.path.join(
                    self.output_directory, "report_" + file_name + ".pdf"
                )

        if self.outputs:
            # TODO change inheritance for report et mni
            for k in (
                "native_t1",
                "native_structures",
                "native_mask",
                "native_lobes",
                "native_macrostructures",
                "native_tissues",
                "mni_t1",
                "mni_structures",
                "mni_mask",
                "mni_lobes",
                "mni_macrostructures",
                "mni_tissues",
                "report_csv",
                "report_pdf",
            ):
                self.tags_inheritance(
                    in_file=self.in_file,
                    out_file=self.outputs[k],
                )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(AssemblyNet, self).run_process_mia()

        in_directory = os.path.dirname(os.path.realpath(self.in_file))
        file_name = os.path.basename(self.in_file)

        id_u = (
            subprocess.check_output(["id", "-u"])
            .decode("utf-8")
            .replace("\n", "")
        )
        id_g = (
            subprocess.check_output(["id", "-g"])
            .decode("utf-8")
            .replace("\n", "")
        )
        cmd = [
            "docker",
            "run",
            "--rm",
            "--user",
            id_u + ":" + id_g,
            "-v",
            in_directory + ":/data",
            "-v",
            self.output_directory + ":/data_out",
            "volbrain/assemblynet:1.0.0",
            "/data/" + file_name,
            "/data_out/",
        ]
        p = subprocess.Popen(
            cmd,
            shell=False,
            bufsize=-1,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
        )

        print("--------->PID:", p.pid)

        (sdtoutl, stderrl) = p.communicate()
        if str(sdtoutl) != "":
            print("sdtoutl: ", sdtoutl.decode())
        if str(stderrl) != "":
            print("stderrl: ", stderrl.decode())
