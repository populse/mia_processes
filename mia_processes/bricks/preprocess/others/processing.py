# -*- coding: utf-8 -*-
"""The other preprocess library of the mia_processes package.

The purpose of this module is to provide bricks generally necessary for the
pre-processing steps, which are not found in nipype.

:Contains:
    :Class:
        - ApplyBiasCorrection
        - ArtifactMask
        - Binarize
        - ConformImage
        - ConvROI
        - Enhance
        - EstimateSNR
        - GradientThreshold
        - Harmonize
        - IntensityClip
        - Mask
        - NonSteadyStateDetector
        - Resample1
        - Resample2
        - RotationMask
        - Sanitize
        - TSNR
        - TemplateFromTemplateFlow
        - Threshold

    :Function:
        - artifact_mask
        - is_outlier
        - threshold
"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

import os
import shutil
import tempfile

# nibabel import
import nibabel as nib
import nibabel.processing as nibp
import numpy as np

# nipype import
from nipype.interfaces.base import (
    File,
    InputMultiPath,
    OutputMultiPath,
    TraitListObject,
    Undefined,
    traits,
)
from nipype.interfaces.spm.base import ImageFileSPM

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA
from scipy import ndimage as sim
from skimage.morphology import ball
from skimage.transform import resize

# soma-base imports
from soma.qt_gui.qt_backend.Qt import QMessageBox
from statsmodels.robust.scale import mad

# import templateflow to get anatomical templates
from templateflow.api import get as get_template

# mia_processes import
from mia_processes.utils import (
    checkFileExt,
    get_dbFieldValue,
    set_dbFieldValue,
)

# Other import
# from distutils.dir_util import copy_tree


# from populse_mia.software_properties import Config


EXT = {"NIFTI_GZ": "nii.gz", "NIFTI": "nii"}


class ApplyBiasCorrection(ProcessMIA):
    """
    *Apply bias field correction to an input file using the bias image*

    Please, see the complete documentation for the 'ApplyBiasCorrection brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/ApplyMask.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ApplyBiasCorrection, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = "A file"
        bias_image_desc = "Bias image"

        # Outputs description
        out_file_desc = "Out file"

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "bias_image",
            File(output=False, optional=False, desc=bias_image_desc),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(ApplyBiasCorrection, self).list_outputs()

        if self.in_file:
            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print(
                        "\nApplyBiasCorrection brick: The input image "
                        "format is not recognized...!"
                    )
                    return

                self.outputs["out_file"] = os.path.join(
                    self.output_directory, fileName + "_inu." + in_ext
                )
            else:
                print(
                    "ApplyBiasCorrection brick: No output_directory was "
                    "found...!\n"
                )
                return

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""

        super(ApplyBiasCorrection, self).run_process_mia()

        in_file = self.in_file
        out_file = self.out_file
        bias_image = self.bias_image

        img = nib.load(in_file)
        data = np.clip(
            img.get_fdata() * nib.load(bias_image).get_fdata(),
            a_min=0,
            a_max=None,
        )
        out_img = img.__class__(
            data.astype(img.get_data_dtype()),
            img.affine,
            img.header,
        )

        out_img.to_filename(out_file)


class ArtifactMask(ProcessMIA):
    """
    *Computes the artifact mask*

    Please, see the complete documentation for the `ArtifactMask brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/ArtifactMask.html>`_

    adapted from
    https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/anatomical.py#L301
    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ArtifactMask, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = "An existing path file."
        head_mask_desc = "An existing path file."
        rot_mask_desc = "An existing path file."
        nasion_post_mask_desc = "An existing path file."
        suffix_desc = "Suffix of the output image (a string)."

        # Outputs description
        out_hat_mask_desc = (
            "Path of the outputted air mask "
            "(a pathlike object or string"
            "representing a file)."
        )
        out_art_mask_desc = (
            "Path of the outputted artifact mask "
            "(a pathlike object or string"
            "representing a file)."
        )
        out_air_mask_desc = (
            "Path of the outputted air mask "
            "(a pathlike object or string"
            "representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "head_mask",
            File(output=False, optional=False, desc=head_mask_desc),
        )

        self.add_trait(
            "nasion_post_mask",
            File(output=False, optional=False, desc=nasion_post_mask_desc),
        )

        self.add_trait(
            "rot_mask",
            File(
                default=None, output=False, optional=True, desc=rot_mask_desc
            ),
        )

        self.add_trait(
            "suffix",
            traits.String(
                "_mask", output=False, optional=True, desc=suffix_desc
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_hat_mask", File(output=True, desc=out_hat_mask_desc)
        )

        self.add_trait(
            "out_art_mask", File(output=True, desc=out_art_mask_desc)
        )

        self.add_trait(
            "out_air_mask", File(output=True, desc=out_air_mask_desc)
        )

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(ArtifactMask, self).list_outputs()

        # Outputs definition
        if self.in_file:
            if (
                (not self.suffix)
                or (self.suffix.isspace())
                or self.suffix in [Undefined, "<undefined>"]
            ):
                self.suffix = " "

            file_hat = ""
            file_art = ""
            file_air = ""

            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print(
                    "\nArtifactMask brick: The input image format is not "
                    "recognized...!"
                )
                return

            file_hat = os.path.join(
                self.output_directory,
                ("hat_" + fileName + self.suffix.strip() + "." + in_ext),
            )

            file_art = os.path.join(
                self.output_directory,
                ("art_" + fileName + self.suffix.strip() + "." + in_ext),
            )

            file_air = os.path.join(
                self.output_directory,
                ("air_" + fileName + self.suffix.strip() + "." + in_ext),
            )

            self.outputs["out_hat_mask"] = file_hat
            self.outputs["out_art_mask"] = file_art
            self.outputs["out_air_mask"] = file_air

            # tags inheritance (optional)
            if self.outputs:
                self.inheritance_dict[
                    self.outputs["out_hat_mask"]
                ] = self.in_file
                self.inheritance_dict[
                    self.outputs["out_art_mask"]
                ] = self.in_file
                self.inheritance_dict[
                    self.outputs["out_air_mask"]
                ] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ArtifactMask, self).run_process_mia()

        file_name = self.in_file
        head_mask_name = self.head_mask
        rot_mask_name = self.rot_mask
        nasion_post_mask_name = self.nasion_post_mask

        try:
            imnii = nib.load(file_name)
            hmnii = nib.load(head_mask_name)
            if rot_mask_name:
                rmnii = nib.load(rot_mask_name)
            npmnii = nib.load(nasion_post_mask_name)

        except (
            nib.filebasedimages.ImageFileError,
            FileNotFoundError,
            TypeError,
        ) as e:
            print(
                "\nArtifactMask brick: Error with files to mask, during "
                "initialisation: ",
                e,
            )
            imnii = None
            hmnii = None
            rmnii = None
            npmnii = None

        if (
            (imnii is not None)
            and (hmnii is not None)
            and (rmnii is not None)
            and (npmnii is not None)
        ):
            imdata = np.nan_to_num(imnii.get_fdata().astype(np.float32))

            # Remove negative values
            imdata[imdata < 0] = 0

            hmdata = np.asanyarray(nib.load(head_mask_name).dataobj)
            npdata = np.asanyarray(nib.load(nasion_post_mask_name).dataobj)

            # Invert head mask
            airdata = np.ones_like(hmdata, dtype=np.uint8)
            airdata[hmdata == 1] = 0

            # Calculate distance to border
            dist = sim.morphology.distance_transform_edt(airdata)

            # Apply nasion-to-posterior mask
            airdata[npdata == 1] = 0
            dist[npdata == 1] = 0
            dist /= dist.max()

            # Apply rotation mask (if supplied)
            if rot_mask_name != Undefined:
                rotmskdata = np.asanyarray(nib.load(rot_mask_name).dataobj)
                airdata[rotmskdata == 1] = 0

            # Run the artifact detection
            qi1_img = artifact_mask(imdata, airdata, dist)
            hdr = imnii.header.copy()
            hdr.set_data_dtype(np.uint8)

            # Image save
            _, file_name = os.path.split(file_name)
            file_name_no_ext, file_extension = os.path.splitext(file_name)
            if file_extension == ".gz":
                (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                    file_name_no_ext
                )
                if file_extension_2 == ".nii":
                    file_name_no_ext = file_name_no_ext_2
                    file_extension = ".nii.gz"

            art_file_out = os.path.join(
                self.output_directory,
                (
                    "art_"
                    + file_name_no_ext
                    + self.suffix.strip()
                    + file_extension
                ),
            )
            nib.Nifti1Image(qi1_img, imnii.affine, hdr).to_filename(
                art_file_out
            )

            hat_file_out = os.path.join(
                self.output_directory,
                (
                    "hat_"
                    + file_name_no_ext
                    + self.suffix.strip()
                    + file_extension
                ),
            )
            nib.Nifti1Image(airdata, imnii.affine, hdr).to_filename(
                hat_file_out
            )

            air_file_out = os.path.join(
                self.output_directory,
                (
                    "air_"
                    + file_name_no_ext
                    + self.suffix.strip()
                    + file_extension
                ),
            )
            airdata[qi1_img > 0] = 0
            nib.Nifti1Image(airdata, imnii.affine, hdr).to_filename(
                air_file_out
            )


class Binarize(ProcessMIA):
    """
    *Image binarization*

    Please, see the complete documentation for the `Binarize brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/Binarize.html>`_

    adapted from:
    https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/nibabel.py#L92

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Binarize, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_files_desc = (
            "A list of items with string elements corresponding "
            "to existing path files."
        )
        thresh_low_desc = "Lower threshold for binarization (float)."
        prefix_desc = "Prefix of the output images (a string)."
        suffix_desc = "Suffix of the output images (a string)."

        # Outputs description
        out_files_desc = (
            "Path of the scan after application of the"
            "binarization (a pathlike object or string"
            "representing a file, or a list of pathlike"
            "objects or strings representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_files",
            InputMultiPath(
                traits.Either(File(), traits.List(File())),
                output=False,
                desc=in_files_desc,
            ),
        )

        self.add_trait(
            "prefix",
            traits.String("", output=False, optional=True, desc=prefix_desc),
        )
        self.add_trait(
            "suffix",
            traits.String(
                "_bin", output=False, optional=True, desc=suffix_desc
            ),
        )

        self.add_trait(
            "thresh_low",
            traits.Float(
                default=0.0, output=False, optional=True, desc=thresh_low_desc
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_files",
            OutputMultiPath(File(), output=True, desc=out_files_desc),
        )

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Binarize, self).list_outputs()

        # Outputs definition
        if self.in_files:
            files_name = self.in_files

            if (
                (not self.suffix)
                or (self.suffix.isspace())
                or self.suffix in [Undefined, "<undefined>"]
            ):
                self.suffix = " "

            if (
                (not self.prefix)
                or (self.prefix.isspace())
                or (self.prefix in [Undefined, "<undefined>"])
            ):
                self.prefix = " "

            files = []
            flag = True  # If False, suf/pref check will not be performed later

            for file_name1 in files_name:
                valid_ext, in_ext, fileName = checkFileExt(file_name1, EXT)
                if not valid_ext:
                    print(
                        "\nBinarize brick: The input image format is not "
                        "recognized...!"
                    )
                    return

                if self.suffix == " " and self.prefix == " " and flag is True:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Warning)
                    msg.setWindowTitle(
                        "mia_processes - " "Binarize brick Warning!"
                    )
                    msg.setText(
                        "Suffix and prefix input parameters are not "
                        "defined or consist only of one or more white "
                        "spaces.\nThe {0} input parameter will be "
                        "overwritten ...\n Yes or "
                        "Abort?".format(file_name1)
                    )
                    msg.setStandardButtons(
                        QMessageBox.Yes
                        | QMessageBox.YesToAll
                        | QMessageBox.Abort
                    )
                    retval = msg.exec_()

                    if retval != QMessageBox.Abort:
                        print(
                            "\nBinarize brick warning: the out_files output "
                            "parameter is the same as the in_files input "
                            "parameter (suffix and prefix are not defined):"
                            "\n{0} will be overwrited ...".format(file_name1)
                        )

                        if retval == QMessageBox.YesToAll:
                            flag = False
                            print(
                                "\nBinarize brick: YesToAll selected ; end of "
                                "overwrite checks on input images ..."
                            )
                    else:
                        files_name = []
                        print(
                            "\nBinarize brick: Initialization Aborted. Please "
                            "check your input parameters ..."
                        )
                        return

                files.append(
                    os.path.join(
                        self.output_directory,
                        (
                            self.prefix.strip()
                            + fileName
                            + self.suffix.strip()
                            + "."
                            + in_ext
                        ),
                    )
                )

            if files:
                self.outputs["out_files"] = files

            else:
                print(
                    "- Binarize brick: There was no output file deducted "
                    "during initialization. Please check the "
                    "input parameters...!"
                )

        # tags inheritance (optional)
        if self.outputs:
            for key, val in self.outputs.items():
                if key == "out_files":
                    for in_val, out_val in zip(files_name, val):
                        _, fileOval = os.path.split(out_val)
                        fileOval_no_ext, file_extension = os.path.splitext(
                            fileOval
                        )

                        if file_extension == ".gz":
                            (
                                fileOval_no_ext_2,
                                file_extension_2,
                            ) = os.path.splitext(fileOval_no_ext)
                            if file_extension_2 == ".nii":
                                fileOval_no_ext = fileOval_no_ext_2

                        _, fileIval = os.path.split(in_val)
                        fileIval_no_ext, file_extension = os.path.splitext(
                            fileIval
                        )

                        if file_extension == ".gz":
                            (
                                fileIval_no_ext_2,
                                file_extension_2,
                            ) = os.path.splitext(fileIval_no_ext)
                            if file_extension_2 == ".nii":
                                fileIval_no_ext = fileIval_no_ext_2

                        if (self.prefix) and (
                            fileOval_no_ext.startswith(self.prefix)
                        ):
                            # fmt: off
                            fileOval_no_ext = fileOval_no_ext[len(
                                                              self.prefix):]
                            # fmt: on

                        if (self.suffix) and (
                            fileOval_no_ext.endswith(self.suffix)
                        ):
                            fileOval_no_ext = fileOval_no_ext[
                                : -len(self.suffix)
                            ]

                        if fileOval_no_ext == fileIval_no_ext:
                            self.inheritance_dict[out_val] = in_val

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Binarize, self).run_process_mia()

        files_name = self.in_files

        for file_name in files_name:
            # Image processing
            try:
                img = nib.load(file_name)

            except (
                nib.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                print(
                    "\nBinarize brick: Error with files to binarize, during "
                    "the run: ",
                    e,
                )
                img = None

            if img is not None:
                data = img.get_fdata()
                mask = data > self.thresh_low
                data[~mask] = 0.0
                img.header.set_data_dtype("uint8")
                maskimg = img.__class__(
                    mask.astype("uint8"), img.affine, img.header
                )

                # Image save
                _, file_name = os.path.split(file_name)
                file_name_no_ext, file_extension = os.path.splitext(file_name)
                if file_extension == ".gz":
                    (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                        file_name_no_ext
                    )
                    if file_extension_2 == ".nii":
                        file_name_no_ext = file_name_no_ext_2
                        file_extension = ".nii.gz"

                file_out = os.path.join(
                    self.output_directory,
                    (
                        self.prefix.strip()
                        + file_name_no_ext
                        + self.suffix.strip()
                        + file_extension
                    ),
                )
                nib.save(maskimg, file_out)


class ConformImage(ProcessMIA):
    """
    *Conform T1w image to standard*

    Please, see the complete documentation for the `ConformImage brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/ConformImage.html>`_

    adapted from:
    https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/common/conform_image.py#L75

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ConformImage, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = "An existing path file."
        prefix_desc = "Prefix of the output image (a string)."
        suffix_desc = "Suffix of the output image (a string)."

        # Outputs description
        out_file_desc = (
            "Path of the conformed scan "
            "(a pathlike object or string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "prefix",
            traits.String("", output=False, optional=True, desc=prefix_desc),
        )

        self.add_trait(
            "suffix",
            traits.String("", output=False, optional=True, desc=suffix_desc),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(ConformImage, self).list_outputs()

        # Outputs definition
        if self.in_file:
            if (
                (not self.suffix)
                or (self.suffix.isspace())
                or self.suffix in [Undefined, "<undefined>"]
            ):
                self.suffix = " "

            if (
                (not self.prefix)
                or (self.prefix.isspace())
                or (self.prefix in [Undefined, "<undefined>"])
            ):
                self.prefix = " "

            valid_ext, in_ext, file_name = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print(
                    "\nConformImage brick: The input image format is not "
                    "recognized...!"
                )
                return

            path, _ = os.path.split(self.in_file)

            if (
                self.suffix == " "
                and self.prefix == " "
                and path == self.output_directory
            ):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle(
                    "mia_processes - ConformImage brick Warning!"
                )
                msg.setText(
                    "Suffix and prefix input parameters are not "
                    "defined or consist only of one or more white "
                    "spaces.\nThe {0} input parameter will be "
                    "overwritten ...\n Yes or "
                    "Abort?".format(file_name)
                )
                msg.setStandardButtons(QMessageBox.Yes | QMessageBox.Abort)
                retval = msg.exec_()

                if retval != QMessageBox.Abort:
                    print(
                        "\nConformImage brick warning: the out_file output "
                        "parameter is the same as the in_files input "
                        "parameter (suffix and prefix are not defined):"
                        "\n{0} will be overwritten ...".format(file_name)
                    )
                else:
                    print("\nAborted. Please check your input parameters ...")
                    return

            file = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + file_name
                    + self.suffix.strip()
                    + "."
                    + in_ext
                ),
            )
            if file:
                self.outputs["out_file"] = file

            else:
                print(
                    "- ConformImage brick: There was no output file deducted "
                    "during initialization. Please check the input "
                    "parameters...!"
                )

            # tags inheritance (optional)
            if self.outputs["out_file"]:
                self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(ConformImage, self).run_process_mia()

        file_name = self.in_file

        # Image processing
        try:
            img = nib.load(file_name)

        except (
            nib.filebasedimages.ImageFileError,
            FileNotFoundError,
            TypeError,
        ) as e:
            print(
                "\nConformImage brick: Error with file to conform, during the "
                "run: ",
                e,
            )
            img = None

        if img is not None:
            out_img = nib.squeeze_image(img)
            out_img = nib.as_closest_canonical(out_img)

            # Image save
            _, file_name = os.path.split(file_name)
            file_name_no_ext, file_extension = os.path.splitext(file_name)
            if file_extension == ".gz":
                (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                    file_name_no_ext
                )
                if file_extension_2 == ".nii":
                    file_name_no_ext = file_name_no_ext_2
                    file_extension = ".nii.gz"

            file_out = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + file_name_no_ext
                    + self.suffix.strip()
                    + file_extension
                ),
            )
            nib.save(out_img, file_out)


class ConvROI(ProcessMIA):
    """
    *Image convolution with one image*

    - Resampling the convolve_with to the size of images_to_convolve.
    - Then convolve each element of images_to_convolve with resized
      convolve_with.
    - The output_directory/PatientName_data/ROI_data/convROI_BOLD directory is
      created to receive the convolved images. If this directory exists at
      runtime it is deleted.
    - To work correctly, the database entry for the convolve_with
      parameter must have the "PatientName" tag filled in.

    Please, see the complete documentation for the `ConvROI brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/ConvROI.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(ConvROI, self).__init__()

        # Inputs description
        images_to_convolve_desc = (
            "A list of images to convolve with convolve_with image"
        )
        convolve_with_desc = "Aa string representing a path to the file"
        prefix_desc = "The prefix for the out_images (a string)"

        # Outputs description
        out_images_desc = "The convoluted images"

        # Inputs traits
        self.add_trait(
            "images_to_convolve",
            InputMultiPath(
                ImageFileSPM(),
                output=False,
                optional=False,
                desc=images_to_convolve_desc,
            ),
        )
        self.add_trait(
            "convolve_with",
            ImageFileSPM(
                output=False, optional=False, desc=convolve_with_desc
            ),
        )

        self.add_trait(
            "prefix",
            traits.String(
                "conv", output=False, optional=True, desc=prefix_desc
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_images",
            OutputMultiPath(File(), output=True, desc=out_images_desc),
        )

        # Special parameter used as a messenger for the run_process_mia method
        self.add_trait(
            "dict4runtime",
            traits.Dict(output=False, optional=True, userlevel=1),
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
        super(ConvROI, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if (
            self.images_to_convolve != Undefined
            and self.convolve_with != Undefined
        ):
            patient_name = get_dbFieldValue(
                self.project, self.convolve_with, "PatientName"
            )

            if patient_name is None:
                print(
                    '\nConvROI brick:\n The "PatientName" tag is not filled '
                    "in the database for the {} file ...\n The calculation"
                    "is aborted...".format(self.convolve_with)
                )
                return self.make_initResult()

            self.dict4runtime["patient_name"] = patient_name
            conv_dir = os.path.join(
                self.output_directory,
                patient_name + "_data",
                "ROI_data",
                "convROI_BOLD",
            )
            list_out = []

            if self.prefix.isspace():
                self.prefix = ""

            for roi_file in self.images_to_convolve:
                list_out.append(
                    os.path.join(
                        conv_dir, self.prefix + os.path.basename(roi_file)
                    )
                )

            self.outputs["out_images"] = list_out

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""

        # No need the next line (we don't use self.process et SPM)
        # super(ConvROI, self).run_process_mia()

        pat_name_dir = os.path.join(
            self.output_directory, self.dict4runtime["patient_name"] + "_data"
        )

        if not os.path.exists(pat_name_dir):
            os.mkdir(pat_name_dir)

        roi_data_dir = os.path.join(pat_name_dir, "ROI_data")

        if not os.path.exists(roi_data_dir):
            os.mkdir(roi_data_dir)

        conv_dir = os.path.join(roi_data_dir, "convROI_BOLD")
        tmp = "None"

        if os.path.isdir(conv_dir):
            tmp = tempfile.mktemp(dir=os.path.dirname(roi_data_dir))
            os.mkdir(tmp)
            shutil.move(conv_dir, os.path.join(tmp, "convROI_BOLD"))
            print(
                '\nConvROI brick:\nA "{}" folder already exists, '
                "it will be overwritten by this new "
                "calculation...".format(conv_dir)
            )
        os.mkdir(conv_dir)

        if os.path.isdir(tmp):
            shutil.rmtree(tmp)

        # Resample the convolve_with to the size of images_to_convolve, then
        # convolve each images_to_convolve with resized convolve_with.
        mask_thresh = threshold(self.convolve_with, 0.5).get_fdata()

        for roi_file in self.images_to_convolve:
            roi_img = nib.load(roi_file)
            roi_data = roi_img.get_fdata()
            roi_size = roi_data.shape[:3]
            resized_mask = resize(mask_thresh, roi_size)
            mult = (roi_data * resized_mask).astype(float)
            # TODO: Should we take info from images_to_convolve or from
            #       convolve_with ?
            #       Currently we take from images_to_convolve
            mult_img = nib.Nifti1Image(mult, roi_img.affine, roi_img.header)

            # Image save in conv_dir
            out_file = os.path.join(
                conv_dir, self.prefix + os.path.basename(roi_file)
            )
            nib.save(mult_img, out_file)
            print("{0} saved".format(out_file))


class Enhance(ProcessMIA):
    """
    *Image enhancing*

    Please, see the complete documentation for the `Enhance brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/Enhance.html>`_

    adapted from:
    https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/workflows/anatomical.py#L974

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Enhance, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_files_desc = (
            "A list of items with string elements corresponding "
            "to existing path files."
        )
        prefix_desc = "Prefix of the output image (a string)."
        suffix_desc = "Suffix of the output image (a string)."

        # Outputs description
        out_files_desc = (
            "Path of the scan after enhancing "
            "(a pathlike object or string representing a file, or "
            "a list of pathlike objects or strings representing a "
            "file)."
        )

        # Inputs traits
        self.add_trait(
            "in_files",
            InputMultiPath(
                traits.Either(File(), traits.List(File())),
                output=False,
                desc=in_files_desc,
            ),
        )

        self.add_trait(
            "prefix",
            traits.String("", output=False, optional=True, desc=prefix_desc),
        )

        self.add_trait(
            "suffix",
            traits.String(
                "_enh", output=False, optional=True, desc=suffix_desc
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_files",
            OutputMultiPath(File(), output=True, desc=out_files_desc),
        )

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Enhance, self).list_outputs()

        # Outputs definition
        if self.in_files:
            files_name = self.in_files

            if (
                (not self.suffix)
                or (self.suffix.isspace())
                or self.suffix in [Undefined, "<undefined>"]
            ):
                self.suffix = " "

            if (
                (not self.prefix)
                or (self.prefix.isspace())
                or (self.prefix in [Undefined, "<undefined>"])
            ):
                self.prefix = " "

            files = []
            flag = True  # If False, suf/pref check will not be performed later

            for file_name1 in files_name:
                valid_ext, in_ext, fileName = checkFileExt(file_name1, EXT)

                if not valid_ext:
                    print(
                        "\nEnhance brick: The input image format is not "
                        "recognized...!"
                    )
                    return
                path, _ = os.path.split(file_name1)

                if (
                    self.suffix == " "
                    and self.prefix == " "
                    and path == self.output_directory
                    and flag is True
                ):
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Warning)
                    msg.setWindowTitle(
                        "mia_processes - Enhance brick Warning!"
                    )
                    msg.setText(
                        "Suffix and prefix input parameters are not "
                        "defined or consist only of one or more white "
                        "spaces.\nThe {0} input parameter will be "
                        "overwritten ...\n Yes or "
                        "Abort?".format(file_name1)
                    )
                    msg.setStandardButtons(
                        QMessageBox.Yes
                        | QMessageBox.YesToAll
                        | QMessageBox.Abort
                    )
                    retval = msg.exec_()

                    if retval != QMessageBox.Abort:
                        print(
                            "\nEnhance brick warning: the out_files output "
                            "parameter is the same as the in_files input "
                            "parameter (suffix and prefix are not defined):"
                            "\n{0} will be overwrited ...".format(file_name1)
                        )

                        if retval == QMessageBox.YesToAll:
                            flag = False
                            print(
                                "\nEnhance brick brick: YesToAll selected ; "
                                "end of overwrite checks on input images ..."
                            )
                    else:
                        files_name = []
                        print(
                            "\nEnhance brick brick: Aborted ; Please check "
                            "your input parameters ..."
                        )
                        return

                files.append(
                    os.path.join(
                        self.output_directory,
                        (
                            self.prefix.strip()
                            + fileName
                            + self.suffix.strip()
                            + "."
                            + in_ext
                        ),
                    )
                )

            if files_name:
                self.outputs["out_files"] = files

            else:
                print(
                    "- Enhance brick: There was no output file deducted "
                    "during initialisation. Please check the input "
                    "parameters...!"
                )

        # tags inheritance (optional)
        if self.outputs:
            for key, val in self.outputs.items():
                if key == "out_files":
                    for in_val, out_val in zip(files_name, val):
                        _, fileOval = os.path.split(out_val)
                        fileOval_no_ext, file_extension = os.path.splitext(
                            fileOval
                        )
                        if file_extension == ".gz":
                            (
                                fileOval_no_ext_2,
                                file_extension_2,
                            ) = os.path.splitext(fileOval_no_ext)
                            if file_extension_2 == ".nii":
                                fileOval_no_ext = fileOval_no_ext_2

                        _, fileIval = os.path.split(in_val)
                        fileIval_no_ext, file_extension = os.path.splitext(
                            fileIval
                        )
                        if file_extension == ".gz":
                            (
                                fileIval_no_ext_2,
                                file_extension_2,
                            ) = os.path.splitext(fileIval_no_ext)
                            if file_extension_2 == ".nii":
                                fileIval_no_ext = fileIval_no_ext_2

                        if (self.prefix) and (
                            fileOval_no_ext.startswith(self.prefix)
                        ):
                            # fmt: off
                            fileOval_no_ext = fileOval_no_ext[len(
                                                              self.prefix):]
                            # fmt: on
                        if (self.suffix) and (
                            fileOval_no_ext.endswith(self.suffix)
                        ):
                            fileOval_no_ext = fileOval_no_ext[
                                : -len(self.suffix)
                            ]

                        if fileOval_no_ext == fileIval_no_ext:
                            self.inheritance_dict[out_val] = in_val

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Enhance, self).run_process_mia()

        files_name = self.in_files

        for file_name in files_name:
            # Image processing
            try:
                img = nib.load(file_name)

            except (
                nib.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                print(
                    "\nEnhance brick: Error with file to enhance, during the "
                    "run: ",
                    e,
                )
                img = None

            if img is not None:
                data = img.get_fdata().astype(np.float32)
                range_max = np.percentile(data[data > 0], 99.98)
                range_min = np.median(data[data > 0])

                # Resample signal excess pixels
                excess = np.where(data > range_max)
                data[excess] = 0
                data[excess] = np.random.choice(
                    data[data > range_min], size=len(excess[0])
                )

                out_image = img.__class__(data, img.affine, img.header)

                # Image save
                _, file_name = os.path.split(file_name)
                file_name_no_ext, file_extension = os.path.splitext(file_name)
                if file_extension == ".gz":
                    (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                        file_name_no_ext
                    )
                    if file_extension_2 == ".nii":
                        file_name_no_ext = file_name_no_ext_2
                        file_extension = ".nii.gz"

                file_out = os.path.join(
                    self.output_directory,
                    (
                        self.prefix.strip()
                        + file_name_no_ext
                        + self.suffix.strip()
                        + file_extension
                    ),
                )
                nib.save(out_image, file_out)


class EstimateSNR(ProcessMIA):
    """
    *Estimate SNR using a segmentation file*

    Please, see the complete documentation for the `EstimateSNR brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/EstimateSNR.html>`_

    adapted from:
    https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/workflows/anatomical.py#L970

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(EstimateSNR, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = (
            "Input file (a pathlike object or string " "representing a file)."
        )
        seg_file_desc = (
            "A segmentation file to calculate SNR (a pathlike "
            "object or string representing a file)."
        )
        # Outputs description
        out_snr_desc = (
            "Estimated SNR" "(a pathlike object or string representing a file"
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "seg_file", File(output=False, optional=False, desc=seg_file_desc)
        )

        # Outputs traits
        self.add_trait("out_snr", traits.Float(output=True, desc=out_snr_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(EstimateSNR, self).list_outputs()

        # Outputs definition
        if self.in_file:
            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print(
                        "\nEstimateSNR brick: The input image format is not "
                        "recognized...!"
                    )
                    return

            else:
                print("EstimateSNR brick: No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(EstimateSNR, self).run_process_mia()

        file_name = self.in_file
        seg_file_name = self.seg_file

        try:
            img = nib.load(file_name)
            seg_img = nib.load(seg_file_name)
        except (
            nib.filebasedimages.ImageFileError,
            FileNotFoundError,
            TypeError,
        ) as e:
            print("\nEstimateSNR brick: Error with files, during the run: ", e)
            return

        data = img.get_fdata()
        mask = seg_img.get_fdata() == 2  # WM label
        self.out_snr = float(
            np.mean(data[mask])
            / (data[mask].std() * np.sqrt(mask.sum() / (mask.sum() - 1)))
        )


class GradientThreshold(ProcessMIA):
    """
    *Computes a threshold from the histogram of the magnitude gradient image*

    Please, see the complete documentation for the `GradientThreshold brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/GradientThreshold.html>`_

    adapted from:
    https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/workflows/anatomical.py#L1039

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(GradientThreshold, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = "An existing path file."
        seg_file_desc = "An existing path of segmentation file."
        prefix_desc = "Prefix of the output image (a string)."
        suffix_desc = "Suffix of the output image (a string)."

        # Outputs description
        out_file_desc = (
            "Path of the thresholded scan "
            "(a pathlike object or string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "seg_file", File(output=False, optional=False, desc=seg_file_desc)
        )

        self.add_trait(
            "prefix",
            traits.String("", output=False, optional=True, desc=prefix_desc),
        )
        self.add_trait(
            "suffix",
            traits.String(
                "_grad", output=False, optional=True, desc=suffix_desc
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(GradientThreshold, self).list_outputs()

        # Outputs definition
        if self.in_file:
            if (
                (not self.suffix)
                or (self.suffix.isspace())
                or self.suffix in [Undefined, "<undefined>"]
            ):
                self.suffix = " "

            if (
                (not self.prefix)
                or (self.prefix.isspace())
                or (self.prefix in [Undefined, "<undefined>"])
            ):
                self.prefix = " "

            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print(
                    "\nGradientThreshold brick: The input image format is "
                    "not recognized...!"
                )
                return

            path, _ = os.path.split(self.in_file)

            if (
                self.suffix == " "
                and self.prefix == " "
                and path == self.output_directory
            ):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle(
                    "mia_processes - " "GradientThreshold brick Warning!"
                )
                msg.setText(
                    "Suffix and prefix input parameters are not "
                    "defined or consist only of one or more white "
                    "spaces.\nThe {0} input parameter will be "
                    "overwritten ...\n Yes or "
                    "Abort?".format(fileName)
                )
                msg.setStandardButtons(QMessageBox.Yes | QMessageBox.Abort)
                retval = msg.exec_()

                if retval != QMessageBox.Abort:
                    print(
                        "\nGradientThreshold brick warning: the out_file"
                        "output parameter is the same as the in_file input "
                        "parameter (suffix and prefix are not defined):"
                        "\n{0} will be overwritten ...".format(fileName)
                    )

                else:
                    print(
                        "\nGradientThreshold brick Aborted. Please check "
                        "your input parameters ..."
                    )
                    return

            file = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + fileName
                    + self.suffix.strip()
                    + "."
                    + in_ext
                ),
            )

            if file:
                self.outputs["out_file"] = file

            else:
                print(
                    "- GradientThreshold brick: There was no output file "
                    "deducted during initialization. Please check the input "
                    "parameters...!"
                )

            # tags inheritance (optional)
            if self.outputs["out_file"]:
                self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(GradientThreshold, self).run_process_mia()

        file_name = self.in_file
        seg_file_name = self.seg_file

        # Image processing
        try:
            img = nib.load(file_name)
            seg_img = nib.load(seg_file_name)

        except (
            nib.filebasedimages.ImageFileError,
            FileNotFoundError,
            TypeError,
        ) as e:
            print(
                "\nGradientThreshold brick: Error with file to enhance, "
                "during the run: ",
                e,
            )
            img = None
            seg_img = None

        if (img is not None) and (seg_img is not None):
            # Image gradient
            data = img.get_fdata().astype(np.float32)
            data_max = np.percentile(data.reshape(-1), 99.5)
            data *= 100 / data_max
            grad = sim.gaussian_gradient_magnitude(data, 3.0)
            grad_max = np.percentile(grad.reshape(-1), 99.5)
            grad *= 100.0
            grad /= grad_max

            # Gradient threshold
            struc = sim.iterate_structure(
                sim.generate_binary_structure(3, 2), 2
            )
            mask = np.zeros_like(grad, dtype=np.uint8)
            mask[grad > 15.0] = 1

            seg_data = seg_img.get_fdata().astype(np.uint8)
            seg_data[seg_data > 0] = 1
            seg_data = sim.binary_dilation(
                seg_data, struc, iterations=2, border_value=1
            ).astype(np.uint8)
            mask[seg_data > 0] = 1
            mask = sim.binary_closing(mask, struc, iterations=2).astype(
                np.uint8
            )

            label_im, nb_labels = sim.label(mask)
            artmsk = np.zeros_like(mask)
            if nb_labels > 2:
                sizes = sim.sum(mask, label_im, list(range(nb_labels + 1)))
                ordered = list(
                    reversed(sorted(zip(sizes, list(range(nb_labels + 1)))))
                )
                for _, label in ordered[2:]:
                    mask[label_im == label] = 0
                    artmsk[label_im == label] = 1

            mask = sim.binary_fill_holes(mask, struc).astype(np.uint8)

            out_image = img.__class__(mask, img.affine, img.header)

            # Image save
            _, file_name = os.path.split(file_name)
            file_name_no_ext, file_extension = os.path.splitext(file_name)
            if file_extension == ".gz":
                (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                    file_name_no_ext
                )
                if file_extension_2 == ".nii":
                    file_name_no_ext = file_name_no_ext_2
                    file_extension = ".nii.gz"

            file_out = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + file_name_no_ext
                    + self.suffix.strip()
                    + file_extension
                ),
            )
            nib.save(out_image, file_out)


class Harmonize(ProcessMIA):
    """
    *Harmonize input image using a white matter mask*

    Please, see the complete documentation for the `Harmonize brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/Harmonize.html>`_

    adapted from:
    https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/anatomical.py#L405

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Harmonize, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = "An existing path file."
        wm_mask_desc = "An existing path of a WM mask file."
        erodemask_desc = "Boolean (a string)."
        suffix_desc = "Suffix of the output image (a string)."
        prefix_desc = "Prefix of the output image (a string)."

        # Outputs description
        out_file_desc = (
            "Path of the harmonized image "
            "(a pathlike object or string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "wm_mask", File(output=False, optional=False, desc=wm_mask_desc)
        )

        self.add_trait(
            "erodemask",
            traits.Bool(
                True, output=False, optional=True, desc=erodemask_desc
            ),
        )
        self.add_trait(
            "prefix",
            traits.String("", output=False, optional=True, desc=prefix_desc),
        )
        self.add_trait(
            "suffix",
            traits.String(
                "_harmonized", output=False, optional=True, desc=suffix_desc
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Harmonize, self).list_outputs()

        # Outputs definition
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print(
                    "\nHarmonize brick: The input image format is "
                    "not recognized...!"
                )
                return

            if (
                (not self.suffix)
                or (self.suffix.isspace())
                or self.suffix in [Undefined, "<undefined>"]
            ):
                self.suffix = " "

            if (
                (not self.prefix)
                or (self.prefix.isspace())
                or (self.prefix in [Undefined, "<undefined>"])
            ):
                self.prefix = " "

            path, _ = os.path.split(self.in_file)

            if (
                self.suffix == " "
                and self.prefix == " "
                and path == self.output_directory
            ):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle(
                    "mia_processes - " "Harmonize brick Warning!"
                )
                msg.setText(
                    "Suffix and prefix input parameters are not "
                    "defined or consist only of one or more white "
                    "spaces.\nThe {0} input parameter will be "
                    "overwritten ...\n Yes or "
                    "Abort?".format(fileName)
                )
                msg.setStandardButtons(QMessageBox.Yes | QMessageBox.Abort)
                retval = msg.exec_()

                if retval != QMessageBox.Abort:
                    print(
                        "\nHarmonize brick warning: the out_file output "
                        "parameter is the same as the in_file input "
                        "parameter (suffix and prefix are not defined):"
                        "\n{0} will be overwrited ...".format(fileName)
                    )

                else:
                    print(
                        "\nHarmonize brick: Aborted. Please check your "
                        "input parameters ..."
                    )
                    return

            file = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + fileName
                    + self.suffix.strip()
                    + "."
                    + in_ext
                ),
            )

            if file:
                self.outputs["out_file"] = file

            else:
                print(
                    "- Harmonize brick: There was no output file deducted "
                    "during initialization. Please check the input "
                    "parameters...!"
                )

            # tags inheritance (optional)
            if self.outputs["out_file"]:
                self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Harmonize, self).run_process_mia()

        try:
            img = nib.load(self.in_file)
            wm_img = nib.load(self.wm_mask).get_fdata()
        except (
            nib.filebasedimages.ImageFileError,
            FileNotFoundError,
            TypeError,
        ) as e:
            print(
                "\nHarmonize brick: Error with file to enhance, during the "
                "run: ",
                e,
            )
            img = None
            wm_img = None

        if (img is not None) and (wm_img is not None):
            wm_img[wm_img < 0.9] = 0
            wm_img[wm_img > 0] = 1
            wm_mask = wm_img.astype(np.uint8)

            if self.erodemask:
                # Create a structural element to be used
                # in an opening operation.
                struc = sim.generate_binary_structure(3, 2)
                # Perform an opening operation on the background data.
                wm_mask = sim.binary_erosion(wm_mask, structure=struc).astype(
                    np.uint8
                )

            data = img.get_fdata()
            data = data * (1000.0 / np.median(data[wm_mask > 0]))

            out_img = img.__class__(data, img.affine, img.header)

            # Image save
            _, file_name = os.path.split(self.in_file)
            file_name_no_ext, file_extension = os.path.splitext(file_name)
            if file_extension == ".gz":
                (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                    file_name_no_ext
                )
                if file_extension_2 == ".nii":
                    file_name_no_ext = file_name_no_ext_2
                    file_extension = ".nii.gz"

            file_out = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + file_name_no_ext
                    + self.suffix.strip()
                    + file_extension
                ),
            )
            nib.save(out_img, file_out)


class IntensityClip(ProcessMIA):
    """
    *Clip the intensity range as prescribed by the percentiles*

    Please, see the complete documentation for the `IntensityClip brick in
    the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/IntensityClip.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(IntensityClip, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = "3D file which intensity will be clipped"
        p_min_desc = "Percentile for the lower bound"
        p_max_desc = "Percentile for the upper bound"
        dtype_desc = "Output datatype"
        nonnegative_desc = "Whether input intensities must be positive"
        invert_desc = "Finalize by inverting contrast"

        # Outputs description
        out_file_desc = "File after clipping"

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits
        self.add_trait(
            "dtype",
            traits.Enum(
                "int16",
                "float32",
                "uint8",
                output=False,
                optional=True,
                desc=dtype_desc,
            ),
        )

        self.add_trait(
            "invert",
            traits.Bool(
                default_value=False,
                output=False,
                optional=True,
                desc=invert_desc,
            ),
        )

        self.add_trait(
            "nonnegative",
            traits.Bool(
                default_value=True,
                output=False,
                optional=True,
                desc=nonnegative_desc,
            ),
        )

        self.add_trait(
            "p_max",
            traits.Float(
                default_value=99.9,
                output=False,
                optional=True,
                desc=p_max_desc,
            ),
        )

        self.add_trait(
            "p_min",
            traits.Float(
                default_value=10.0,
                output=False,
                optional=True,
                desc=p_min_desc,
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(IntensityClip, self).list_outputs()

        if self.in_file:
            if self.output_directory:
                self.outputs["out_file"] = os.path.join(
                    self.output_directory,
                    os.path.split(self.in_file)[1].replace(
                        ".nii", "_clipped.nii"
                    ),
                )
            else:
                print(
                    "IntensityClip brick: No output_directory was "
                    "found...!\n"
                )
                return

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""

        super(IntensityClip, self).run_process_mia()

        in_file = self.in_file
        out_file = self.out_file
        nonnegative = self.nonnegative
        p_min = self.p_min
        p_max = self.p_max
        invert = self.invert
        dtype = self.dtype

        # Load data
        img = nib.squeeze_image(nib.load(in_file))
        if len(img.shape) != 3:
            raise RuntimeError(f"<{in_file}> is not a 3D file.")
        data = img.get_fdata(dtype="float32")

        # Calculate stats on denoised version, to preempt outliers from biasing
        denoised = sim.median_filter(data, footprint=ball(3))

        a_min = np.percentile(
            denoised[denoised > 0] if nonnegative else denoised, p_min
        )
        a_max = np.percentile(
            denoised[denoised > 0] if nonnegative else denoised, p_max
        )

        # Clip and cast
        data = np.clip(data, a_min=a_min, a_max=a_max)
        data -= data.min()
        data /= data.max()

        if invert:
            data = 1.0 - data

        if dtype in ("uint8", "int16"):
            data = np.round(255 * data).astype(dtype)

        hdr = img.header.copy()
        hdr.set_data_dtype(dtype)
        img.__class__(data, img.affine, hdr).to_filename(out_file)


class Mask(ProcessMIA):
    """
    *Apply a binary mask to an image*

    Please, see the complete documentation for the `Mask brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/Mask.html>`_

    adapted from:
    https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/norm.py#L474

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Mask, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = "An existing path file."
        mask_file_desc = "An existing path file."
        prefix_desc = "Prefix of the output image (a string)."
        suffix_desc = "Suffix of the output image (a string)."

        # Outputs description
        out_file_desc = (
            "Path of the scan after masking "
            "(a pathlike object or string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "mask_file",
            File(output=False, optional=False, desc=mask_file_desc),
        )

        self.add_trait(
            "prefix",
            traits.String("", output=False, optional=True, desc=prefix_desc),
        )

        self.add_trait(
            "suffix",
            traits.String(
                "_masked", output=False, optional=True, desc=suffix_desc
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Mask, self).list_outputs()

        # Outputs definition
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print(
                    "\nMask brick: The input image format is not "
                    "recognized...!"
                )
                return

            if (
                (not self.suffix)
                or (self.suffix.isspace())
                or self.suffix in [Undefined, "<undefined>"]
            ):
                self.suffix = " "

            if (
                (not self.prefix)
                or (self.prefix.isspace())
                or (self.prefix in [Undefined, "<undefined>"])
            ):
                self.prefix = " "

            path, _ = os.path.split(self.in_file)

            if (
                self.suffix == " "
                and self.prefix == " "
                and path == self.output_directory
            ):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle("mia_processes - " "Mask brick Warning!")
                msg.setText(
                    "Suffix and prefix input parameters are not "
                    "defined or consist only of one or more white "
                    "spaces.\nThe {0} input parameter will be "
                    "overwritten ...\n Yes or "
                    "Abort?".format(fileName)
                )
                msg.setStandardButtons(QMessageBox.Yes | QMessageBox.Abort)
                retval = msg.exec_()

                if retval != QMessageBox.Abort:
                    print(
                        "\nMask brick warning: the out_file output "
                        "parameter is the same as the in_files input "
                        "parameter (suffix and prefix are not defined):"
                        "\n{0} will be overwrited ...".format(fileName)
                    )

                else:
                    print(
                        "\nMask brick: Aborted. Please check your input "
                        "parameters ..."
                    )
                    return

            file = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + fileName
                    + self.suffix.strip()
                    + "."
                    + in_ext
                ),
            )

            if file:
                self.outputs["out_file"] = file

            else:
                print(
                    "- Mask brick: There was no output file deducted during "
                    "initialization. Please check the input parameters...!"
                )

        # tags inheritance (optional)
        if self.outputs["out_file"]:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Mask, self).run_process_mia()

        file_name = self.in_file
        mask_name = self.mask_file

        # Image processing
        try:
            img = nib.load(file_name)
            mask_img = nib.load(mask_name)

        except (
            nib.filebasedimages.ImageFileError,
            FileNotFoundError,
            TypeError,
        ) as e:
            print(
                "\nMask brick: Error with files to mask, during " "the run: ",
                e,
            )
            img = None

        if (img is not None) and (mask_img is not None):
            data = img.get_fdata()
            data[np.asanyarray(mask_img.dataobj) == 0] = 0

            maskimg = img.__class__(data, img.affine, img.header)

            # Image save
            _, file_name = os.path.split(file_name)
            file_name_no_ext, file_extension = os.path.splitext(file_name)
            if file_extension == ".gz":
                (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                    file_name_no_ext
                )
                if file_extension_2 == ".nii":
                    file_name_no_ext = file_name_no_ext_2
                    file_extension = ".nii.gz"

            file_out = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + file_name_no_ext
                    + self.suffix.strip()
                    + file_extension
                ),
            )
            nib.save(maskimg, file_out)


class NonSteadyStateDetector(ProcessMIA):
    """
    *Detect non-steady-state at the beginning of a bold 4D image*

    Please, see the complete documentation for the
    `NonSteadyStateDetector brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/NonSteadyStateDetector.html>`_

    adapted from:
    https://github.com/nipy/nipype/blob/f662acfce8def4717e0c3414618f3a5de5913b31/nipype/algorithms/confounds.py#L974

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(NonSteadyStateDetector, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = "An existing path file."

        # Outputs description
        n_volumes_to_discard_desc = "Number of non steady" "state volumes"

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Outputs traits
        self.add_trait(
            "n_volumes_to_discard",
            traits.Int(output=True, desc=n_volumes_to_discard_desc),
        )

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(NonSteadyStateDetector, self).list_outputs()

        # Outputs definition
        if self.in_file:
            self.outputs["n_volumes_to_discard"] = 0

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(NonSteadyStateDetector, self).run_process_mia()

        file_name = self.in_file

        try:
            in_nii = nib.load(file_name)
        except (
            nib.filebasedimages.ImageFileError,
            FileNotFoundError,
            TypeError,
        ) as e:
            print(
                "\nNonSteadyStateDetector brick: Error with file to enhance, "
                "during the run: ",
                e,
            )
            in_nii = None

        if in_nii is not None:
            global_signal = (
                in_nii.dataobj[:, :, :, :50]
                .mean(axis=0)
                .mean(axis=0)
                .mean(axis=0)
            )

        self.n_volumes_to_discard = is_outlier(global_signal)


class Resample1(ProcessMIA):
    """
    *Resamples an image to the resolution of a reference image*

    - Uses nibabel.processing.resample_from_to().

    Please, see the complete documentation for the
    `Resample1 brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/Resample1.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Resample1, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        files_to_resample_desc = (
            "The 3D images that will be resampled (a "
            "list of pathlike object or string "
            "representing a file or a list of items "
            "which are a pathlike object or string "
            "representing a file, valid extensions: "
            "[.img, .nii, .hdr])."
        )
        reference_image_desc = (
            "A 3D or 4D image used as reference to "
            "resample the files_to_resample images (a "
            "pathlike object or string representing a file "
            "with extension in [.img, .nii, .hdr])."
        )
        interp_desc = (
            "The order of the spline interpolation (an integer "
            "between 0 and 5; trilinear == 3)."
        )
        prefix_desc = "The prefix for the out_files image (a string)."
        suffix_to_delete_desc = (
            "The suffix to delete from the "
            "files_to_resample, when creating the "
            "out_files (a string)."
        )
        suffix_desc = "The suffix for the out_files image (a string)."

        # Outputs description
        out_files_desc = (
            "The resulting images after resampling (a pathlike "
            "object or string representing a file, or a list of "
            "pathlike objects or strings representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "files_to_resample",
            InputMultiPath(
                ImageFileSPM(), output=False, desc=files_to_resample_desc
            ),
        )

        self.add_trait(
            "reference_image",
            ImageFileSPM(output=False, desc=reference_image_desc),
        )

        self.add_trait(
            "interp",
            traits.Range(
                value=3,
                low=0,
                high=5,
                output=False,
                optional=True,
                desc=interp_desc,
            ),
        )

        self.add_trait(
            "prefix",
            traits.String(" ", output=False, optional=True, desc=prefix_desc),
        )

        self.add_trait(
            "suffix_to_delete",
            traits.String(
                "_002", output=False, optional=True, desc=suffix_to_delete_desc
            ),
        )

        self.add_trait(
            "suffix",
            traits.String(
                "_003", output=False, optional=True, desc=suffix_desc
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_files",
            OutputMultiPath(File(), output=True, desc=out_files_desc),
        )

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Resample1, self).list_outputs()

        # Outputs definition
        if (
            self.reference_image != Undefined
            and self.files_to_resample != Undefined
            and self.interp != Undefined
        ):
            files = []
            files_name = self.files_to_resample

            try:
                refName = nib.load(self.reference_image)

            except (
                nib.filebasedimages.ImageFileError,
                FileNotFoundError,
                TypeError,
            ) as e:
                print(
                    "\nResample1 brick: Error with reference_image, during "
                    "initialisation: ",
                    e,
                )
                refName = None

            if (
                (not self.suffix_to_delete)
                or (self.suffix_to_delete.isspace())
                or (self.suffix_to_delete in [Undefined, "<undefined>"])
            ):
                self.suffix_to_delete = " "

            if (
                (not self.suffix)
                or (self.suffix.isspace())
                or (self.suffix in [Undefined, "<undefined>"])
            ):
                self.suffix = " "

            if (
                (not self.prefix)
                or (self.prefix.isspace())
                or (self.prefix in [Undefined, "<undefined>"])
            ):
                self.prefix = " "

            for file_name in files_name:
                try:
                    fileName = nib.load(file_name)

                except (
                    nib.filebasedimages.ImageFileError,
                    FileNotFoundError,
                    TypeError,
                ) as e:
                    print(
                        "\nResample1 brick: Error with files_to_resample, "
                        "during initialisation: ",
                        e,
                    )
                    fileName = None

                if (refName) and (fileName):
                    if (len(fileName.shape) != 3) or (
                        3 > len(refName.shape) > 4
                    ):
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Warning)
                        msg.setWindowTitle(
                            "mia_processes - Resample1 brick Warning!"
                        )
                        msg.setText(
                            "Currently, the Resample1 brick only "
                            "allows to resample 3D images from 3D or "
                            "4D reference images.\n\n However the "
                            "'{0}' reference image is a {1}D and the "
                            "'{2}' image is a {3}D ...\n\nPlease, "
                            "modify your input and initialise again "
                            "this brick.".format(
                                self.reference_image,
                                len(refName.shape),
                                file_name,
                                len(fileName.shape),
                            )
                        )
                        msg.setStandardButtons(QMessageBox.Close)
                        msg.buttonClicked.connect(msg.close)
                        msg.exec()
                        print(
                            "\nResample1 brick: Initialisation failed ... "
                            "Please check the input parameters!"
                        )
                        files_name = None
                        break

                else:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Warning)
                    msg.setWindowTitle(
                        "mia_processes - Resample1 brick Warning!"
                    )
                    msg.setText(
                        "files_to_resample '{0}' or/and "
                        "reference_image '{1}' is (are) empty ... Do "
                        "you want to continue? \n\n"
                        "- To correct the input parameters "
                        "click 'Abort'.\n\n"
                        "- If the Resample1 brick is located in a "
                        "pipeline, it may be usual that this(these) "
                        "parameter(s) is(are) empty at the "
                        "initialisation time. In this case, if the "
                        "files_to_resample and the reference_image "
                        "parameters correspond to a 3D image and a 3D "
                        "or a 4D, respectively, click 'Yes', "
                        "otherwise click 'Abort' and change the input "
                        "parameters.".format(file_name, self.reference_image)
                    )
                    msg.setStandardButtons(QMessageBox.Yes | QMessageBox.Abort)
                    retval = msg.exec_()

                    if retval == QMessageBox.Yes:
                        print(
                            "\nResample1 brick Warning!: Empty input "
                            "file(s). No check of the dimensions of the "
                            "input images is done ... (files_to_resample "
                            "must be a 3D and reference_image must be a 3D "
                            "or a 4D) ..."
                        )

                    else:
                        print(
                            "\nResample1 brick: Initialisation failed ... "
                            "Please check the input parameters!"
                        )
                        files_name = None
                        break

                path, file_name = os.path.split(file_name)
                file_name_no_ext, file_extension = os.path.splitext(file_name)

                # fmt: off
                if (self.suffix_to_delete != " ") and (
                    file_name_no_ext[-len(self.suffix_to_delete):]
                    == self.suffix_to_delete
                ):
                    file_name_no_ext = file_name_no_ext[
                        :-len(self.suffix_to_delete)
                    ]

                # fmt: on
                files.append(
                    os.path.join(
                        self.output_directory,
                        (
                            self.prefix.strip()
                            + file_name_no_ext
                            + self.suffix.strip()
                            + file_extension
                        ),
                    )
                )

            # May overwrite an input image with an output image
            if (self.suffix == " ") and (self.prefix == " ") and (files_name):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle("mia_processes - Resample1 brick Warning!")
                msg.setText(
                    "Suffix and prefix input parameters are not "
                    "defined or consist only of one or more white "
                    "spaces.\n\nThe files_to_resample input parameter "
                    "(if suffix_to_delete input parameter is not set) "
                    "or other files (if suffix_to_delete input "
                    "parameter is set) could be overwritten ... this "
                    "concerns the following images:\n{}\n\nDo you "
                    "agree to use these input "
                    "parameters?".format(set(files))
                )
                msg.setStandardButtons(QMessageBox.Yes | QMessageBox.Abort)
                retval = msg.exec_()

                if retval == QMessageBox.Yes:
                    print(
                        "\nResample1 brick warning: suffix and prefix input "
                        "parameters are not defined, the following files "
                        "could be overwrite!:\n{0} ...\n".format(set(files))
                    )

                else:
                    files_name = None

            # Two (at least) input images give the same output image
            if (len(files) != len(set(files))) and (files_name):
                dupes = set()
                seen = set()

                for x in files:
                    if x in seen:
                        dupes.add(x)

                    seen.add(x)

                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle("mia_processes - Resample1 brick Warning!")
                msg.setText(
                    "The suffix, prefix and suffix_to_delete input "
                    "parameters combination seems to create for at "
                    "least two different input images the same output "
                    "image:\n{}\n\nPlease change the suffix, prefix "
                    "and suffix_to_delete input parameters combination "
                    "and initialise again this brick.".format(dupes)
                )
                msg.setStandardButtons(QMessageBox.Close)
                msg.buttonClicked.connect(msg.close)
                msg.exec()
                files_name = None

            if files_name:
                self.outputs["out_files"] = files

        # Tags inheritance (optional)
        if self.outputs:
            for key, val in self.outputs.items():
                if key == "out_files":
                    for in_val, out_val in zip(files_name, val):
                        _, fileOval = os.path.split(out_val)
                        fileOval_no_ext, _ = os.path.splitext(fileOval)

                        # fmt: off
                        if (self.suffix != " ") and (
                                fileOval_no_ext[-len(
                                                self.suffix):] == self.suffix
                        ):
                            fileOval_no_ext = fileOval_no_ext[:-len(
                                self.suffix)]

                        if (self.prefix != " ") and (
                                fileOval_no_ext[0:len(
                                    self.prefix)] == self.prefix
                        ):
                            fileOval_no_ext = fileOval_no_ext[len(
                                self.prefix):]
                        # fmt: on
                        _, fileIval = os.path.split(in_val)
                        fileIval_no_ext, _ = os.path.splitext(fileIval)

                        # fmt: off
                        if (self.suffix_to_delete != " ") and (
                            fileIval_no_ext[-len(self.suffix_to_delete):]
                            == self.suffix_to_delete
                        ):
                            # fmt: on
                            fileOval_no_ext = (
                                fileOval_no_ext + self.suffix_to_delete
                            )

                        if fileOval_no_ext == fileIval_no_ext:
                            self.inheritance_dict[out_val] = in_val
                            # FIXME: In the latest version of mia, indexing of
                            #        the database with particular tags defined
                            #        in the processes is done only at the end
                            #        of the initialisation of the whole
                            #        pipeline. So we cannot use the value of
                            #        these tags in other processes of the
                            #        pipeline at the time of initialisation
                            #        (see populse_mia #290). Until better we
                            #        use a quick and dirty hack with the
                            #        set_dbFieldValue() function !
                            tag_to_add = dict()
                            tag_to_add["name"] = "PatientName"
                            tag_to_add["field_type"] = "string"
                            tag_to_add["description"] = ""
                            tag_to_add["visibility"] = True
                            tag_to_add["origin"] = "user"
                            tag_to_add["unit"] = None
                            tag_to_add["default_value"] = None
                            tag_to_add["value"] = get_dbFieldValue(
                                self.project, in_val, "PatientName"
                            )

                            if tag_to_add["value"] is not None:
                                set_dbFieldValue(
                                    self.project, out_val, tag_to_add
                                )

                            else:
                                print(
                                    "\nResample1 brick:\nThe 'PatientName' "
                                    "tag could not be added to the database "
                                    "for the '{}' parameter. This can lead to "
                                    "a subsequent issue during "
                                    "initialization!!\n".format(out_val)
                                )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def _check_interp(self):
        """Checks the order of the splines interpolation.

        :returns: a message according to the selected order
        """

        # skimage.transform.resize (does not seem to be the best solution for
        # resampling):
        # - B-splines of the order 0 and 1 correspond to nearest neighbour and
        #   linear interpolation, respectively.
        # - B-splines of a higher order can be defined by a repetitive
        #   convolution of the zeroth-order spline (the box function) with
        #   itself.
        #
        # nibabel.processing.resample_from_to (seems to be a good solution, but
        # can give non-zero intensity for the background -> needs
        # thresolding ?)
        #
        # nilearn.image.resample_to_img (seems to be a good solution, no
        # non-zero intensity for the background observed)

        if self.interp == 0:
            # skimage and nibabel: Nearest neighbour
            return "spline interpolation of order 0 (Nearest-neighbour)"

        elif self.interp == 1:
            # skimage: It seems that with dimensions=3 and B-spline order=1,
            #          the interpolation would be equivalent to "trilinear"
            # return ("Bi-linear (with dimensions=3, Bi-linear is equivalent "
            #        "to 'trilinear')")

            # nibabel: spline interpolation of order 1
            return "spline interpolation of order 1"

        elif self.interp == 2:
            # skimage:
            # return "Bi-quadratic"

            # nibabel: spline interpolation of order 2
            return "spline interpolation of order 2"

        elif self.interp == 3:
            # skimage
            # return "Bi-cubic"

            # nibabel: spline interpolation of order 3: trilinear
            return "spline interpolation of order 3 (trilinear)"

        elif self.interp == 4:
            # skimage
            # return "Bi-quartic"

            # nibabel: spline interpolation of order 4
            return "spline interpolation of order 4"

        elif self.interp == 5:
            # skimage
            # return "Bi-quintic"

            # nibabel: spline interpolation of order 5
            return "spline interpolation of order 5"

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Resample1, self).run_process_mia()
        files_name = self.files_to_resample
        refName = nib.load(self.reference_image)
        print(
            "\nResample1 process from mia_processes: \n",
            self._check_interp(),
            " \n",
        )

        for file_name in files_name:
            fileName = nib.load(file_name)
            # mask_data = mask.get_fdata()  # skimage

            # Currently, resampling is only done for 3D images using 3D
            # or 4D images for reference. It would be interesting to widen the
            # possibilities to more cases (at least a 4D using a 3D, and why
            # not, to larger dimensions?). Currently, we use
            # nibabel.processing.resample_from_to(), (see
            # https://mail.python.org/pipermail/neuroimaging/2019-January/001902.html).
            # It seems this method produces a little noise (can give non-zero
            # intensity for the background). In addition to the resolution
            # settings it seems the size of the reference_image are also
            # transferred to the resampled image but not the orientations and
            # positions (output space comes from the affine of
            # files_to_resample).

            if len(fileName.shape) == len(refName.shape) == 3:
                # nibabel:
                fileFinal = nibp.resample_from_to(
                    fileName, refName, order=self.interp
                )

            if len(fileName.shape) == 3 and len(refName.shape) == 4:
                # skimage
                # ref_size = ref_data.shape[:3]
                # resized_mask_data = resize(mask_data,
                #                           ref_size,
                #                           order=self.interp,
                #                           mode='reflect')
                # TODO: Taking info of mask's or ref's header?
                # mask_final = nib.Nifti1Image(resized_mask_data,
                #                             ref.affine,
                #                             ref.header)

                # nibabel:
                fileFinal = nibp.resample_from_to(
                    fileName, refName.slicer[:, :, :, 0], order=self.interp
                )

            # Image save
            path, file_name = os.path.split(file_name)
            file_name_no_ext, file_extension = os.path.splitext(file_name)

            # fmt: off
            if (self.suffix_to_delete != " ") and (
                file_name_no_ext[-len(self.suffix_to_delete):]
                == self.suffix_to_delete
            ):
                # fmt: on
                file_name_no_ext = file_name_no_ext[
                    : -len(self.suffix_to_delete)
                ]

            file_out = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + file_name_no_ext
                    + self.suffix.strip()
                    + file_extension
                ),
            )
            nib.save(fileFinal, file_out)


class Resample2(ProcessMIA):
    """
    *Resamples images to the resolution of a reference image*

    - Uses skimage.transform.resize()
    - The output_directory/PatientName_data/ROI_data/convROI_BOLD2 directory is
      created to receive the resampled images. If this directory exists at
      runtime it is deleted.
    - To work correctly, the database entry for the reference_image parameter
      must have the "PatientName" tag filled in.

    Please, see the complete documentation for the
    `Resample2 brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/Resample2.html>`_
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        super(Resample2, self).__init__()

        # Inputs description
        files_to_resample_desc = (
            "The images that will be resampled (a "
            "list of pathlike object or string "
            "representing a file or a list of items "
            "which are a pathlike object or string "
            "representing a file, valid extensions: "
            "[.img, .nii, .hdr])."
        )

        reference_image_desc = (
            "The reference image for resampling (an "
            "existing, uncompressed file)"
        )

        suffix_desc = "The suffix for the out_images (a string)"

        # Outputs description
        out_images_desc = "The resampled images"

        # Inputs traits
        self.add_trait(
            "files_to_resample",
            InputMultiPath(
                ImageFileSPM(), output=False, desc=files_to_resample_desc
            ),
        )

        self.add_trait(
            "reference_image",
            ImageFileSPM(output=False, desc=reference_image_desc),
        )

        self.add_trait(
            "suffix",
            traits.String("_2", output=False, optional=True, desc=suffix_desc),
        )

        # Outputs traits
        self.add_trait(
            "out_images",
            OutputMultiPath(File(), output=True, desc=out_images_desc),
        )

        # Special parameter used as a messenger for the run_process_mia method
        self.add_trait(
            "dict4runtime",
            traits.Dict(output=False, optional=True, userlevel=1),
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
        super(Resample2, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if (
            self.files_to_resample != Undefined
            and self.reference_image != Undefined
        ):
            patient_name = get_dbFieldValue(
                self.project, self.reference_image, "PatientName"
            )

            if patient_name is None:
                print(
                    "\nResample2 brick:\n The PatientName tag is not filled "
                    "in the database for the {} file ...\n The calculation"
                    "is aborted...".format(self.reference_image)
                )
                return self.make_initResult()

            self.dict4runtime["patient_name"] = patient_name
            conv_dir2 = os.path.join(
                self.output_directory,
                patient_name + "_data",
                "ROI_data",
                "convROI_BOLD2",
            )
            list_out = []

            for roi_file in self.files_to_resample:
                out_file = os.path.basename(roi_file)
                out_file_no_ext, file_extension = os.path.splitext(out_file)
                out_file = out_file_no_ext + self.suffix + file_extension
                out_file = os.path.join(conv_dir2, out_file)
                list_out.append(out_file)

            self.outputs["out_images"] = list_out
            # FIXME: What are we doing about the tags inheritance?

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""

        # No need the next line (we don't use self.process et SPM)
        # super(Resample2, self).run_process_mia()

        pat_name_dir = os.path.join(
            self.output_directory, self.dict4runtime["patient_name"] + "_data"
        )

        if not os.path.exists(pat_name_dir):
            os.mkdir(pat_name_dir)

        roi_data_dir = os.path.join(pat_name_dir, "ROI_data")

        if not os.path.exists(roi_data_dir):
            os.mkdir(roi_data_dir)

        conv_dir2 = os.path.join(roi_data_dir, "convROI_BOLD2")
        tmp = "None"

        if os.path.isdir(conv_dir2):
            tmp = tempfile.mktemp(dir=os.path.dirname(roi_data_dir))
            os.mkdir(tmp)
            shutil.move(conv_dir2, os.path.join(tmp, "convROI_BOLD2"))
            print(
                '\nResample2 brick:\nA "{}" folder already exists, '
                "it will be overwritten by this new "
                "calculation...".format(conv_dir2)
            )
        os.mkdir(conv_dir2)

        if os.path.isdir(tmp):
            shutil.rmtree(tmp)

        # Setting files_to_resample to the resolution of the reference_image
        mask = nib.load(self.reference_image).get_fdata()
        mask_size = mask.shape[:3]

        for roi_file in self.files_to_resample:
            roi_img = nib.load(roi_file)
            roi_data = roi_img.get_fdata()
            resized_roi = resize(roi_data, mask_size)
            # TODO: Should we take info from ROI or from the mask ?
            #       Currently we take from ROI images
            resized_img = nib.Nifti1Image(
                resized_roi, roi_img.affine, roi_img.header
            )

            # Image save
            out_file = os.path.basename(roi_file)
            out_file_no_ext, file_extension = os.path.splitext(out_file)
            out_file = out_file_no_ext + self.suffix + file_extension
            out_file = os.path.join(conv_dir2, out_file)
            nib.save(resized_img, out_file)
            print(
                "\nResample2 brick:{0} saved".format(
                    os.path.basename(out_file)
                )
            )


class RotationMask(ProcessMIA):
    """
    *Compute the rotation mask image*

    Please, see the complete documentation for the
    `RotationMask brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/RotationMask.html>`_

    adapted from:
    https://github.com/nipreps/mriqc/blob/e021008da0a2ef1c48e882baf932139a673349f9/mriqc/interfaces/anatomical.py#L448

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(RotationMask, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = "An existing path file."
        prefix_desc = "Prefix of the output image (a string)."
        suffix_desc = "Suffix of the output image (a string)."

        # Outputs description
        out_file_desc = (
            "Path of the scan after masking "
            "(a pathlike object or string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "suffix",
            traits.String(
                "_rotmasked", output=False, optional=True, desc=suffix_desc
            ),
        )

        self.add_trait(
            "prefix",
            traits.String("", output=False, optional=True, desc=prefix_desc),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(RotationMask, self).list_outputs()

        # Outputs definition
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print(
                    "\nRotationMask brick: The input image format is not "
                    "recognized...!"
                )
                return

            if (
                (not self.suffix)
                or (self.suffix.isspace())
                or self.suffix in [Undefined, "<undefined>"]
            ):
                self.suffix = " "

            if (
                (not self.prefix)
                or (self.prefix.isspace())
                or (self.prefix in [Undefined, "<undefined>"])
            ):
                self.prefix = " "

            path, _ = os.path.split(self.in_file)

            if (
                self.suffix == " "
                and self.prefix == " "
                and path == self.output_directory
            ):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle(
                    "mia_processes - RotationMask brick Warning!"
                )
                msg.setText(
                    "Suffix and prefix input parameters are not "
                    "defined or consist only of one or more white "
                    "spaces.\nThe {0} input parameter will be "
                    "overwritten ...\n Yes or "
                    "Abort?".format(fileName)
                )
                msg.setStandardButtons(QMessageBox.Yes | QMessageBox.Abort)
                retval = msg.exec_()

                if retval != QMessageBox.Abort:
                    print(
                        "\nRotationMask brick warning: the out_file output "
                        "parameter is the same as the in_files input "
                        "parameter (suffix and prefix are not defined):"
                        "\n{0} will be overwrited ...".format(fileName)
                    )

                else:
                    print(
                        "\nRotationMask brick Aborted. Please check your "
                        "input parameters ..."
                    )
                    return

            file = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + fileName
                    + self.suffix.strip()
                    + "."
                    + in_ext
                ),
            )

            if file:
                self.outputs["out_file"] = file

            else:
                print(
                    "- RotationMask brick: There was no output file deducted "
                    "during initialization. Please check the input "
                    "parameters...!"
                )

            # tags inheritance (optional)
            if self.outputs["out_file"]:
                self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(RotationMask, self).run_process_mia()

        file_name = self.in_file

        # Image processing
        try:
            img = nib.load(file_name)

        except (
            nib.filebasedimages.ImageFileError,
            FileNotFoundError,
            TypeError,
        ) as e:
            print(
                "\nRotationMask brick: Error with files to mask, during the "
                "run: ",
                e,
            )
            img = None

        if img is not None:
            data = img.get_fdata()
            mask = data <= 0

            # Pad one pixel to control behavior on borders of binary_opening
            mask = np.pad(
                mask, pad_width=(1,), mode="constant", constant_values=1
            )

            # Remove noise
            struc = sim.generate_binary_structure(3, 2)
            mask = sim.binary_opening(mask, structure=struc).astype(np.uint8)

            # Remove small objects
            label_im, nb_labels = sim.label(mask)
            if nb_labels > 2:
                sizes = sim.sum(mask, label_im, list(range(nb_labels + 1)))
                ordered = list(
                    reversed(sorted(zip(sizes, list(range(nb_labels + 1)))))
                )
                for _, label in ordered[2:]:
                    mask[label_im == label] = 0

            # Un-pad
            mask = mask[1:-1, 1:-1, 1:-1]

            # If mask is small, clean-up
            if mask.sum() < 500:
                mask = np.zeros_like(mask, dtype=np.uint8)

            out_img = img.__class__(mask, img.affine, img.header)
            out_img.header.set_data_dtype(np.uint8)

            # Image save
            _, file_name = os.path.split(file_name)
            file_name_no_ext, file_extension = os.path.splitext(file_name)
            if file_extension == ".gz":
                (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                    file_name_no_ext
                )
                if file_extension_2 == ".nii":
                    file_name_no_ext = file_name_no_ext_2
                    file_extension = ".nii.gz"

            file_out = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + file_name_no_ext
                    + self.suffix.strip()
                    + file_extension
                ),
            )
            nib.save(out_img, file_out)


class Sanitize(ProcessMIA):
    """
    *Sanitize input bold image*

    Please, see the complete documentation for the
    `Sanitize brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/Sanitize.html>`_

    adapted from:
    https://github.com/nipreps/niworkflows/blob/45ab13e1bf6fdbf5e29c90cef44055b0b9cf391b/niworkflows/interfaces/header.py#L394

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Sanitize, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = "An existing path file."
        prefix_desc = "Prefix of the output image (a string)."
        suffix_desc = "Suffix of the output image (a string)."
        n_volumes_to_discard_desc = "number of non steady-state volumes"
        max_32bit_desc = (
            "cast data to float32 if higher precision is " "encountered"
        )
        # Outputs description
        out_file_desc = (
            "Path of the snaitized scan"
            "(a pathlike object or string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "max_32bit",
            traits.Bool(
                False, output=False, optional=True, desc=max_32bit_desc
            ),
        )

        self.add_trait(
            "n_volumes_to_discard",
            traits.Int(
                0, output=False, optional=True, desc=n_volumes_to_discard_desc
            ),
        )

        self.add_trait(
            "prefix",
            traits.String("", output=False, optional=True, desc=prefix_desc),
        )

        self.add_trait(
            "suffix",
            traits.String(
                "_valid", output=False, optional=True, desc=suffix_desc
            ),
        )

        # Outputs traits
        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Sanitize, self).list_outputs()

        # Outputs definition
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print(
                    "\nSanitize brick: The input image format is not "
                    "recognized...!"
                )
                return

            if (
                (not self.suffix)
                or (self.suffix.isspace())
                or self.suffix in [Undefined, "<undefined>"]
            ):
                self.suffix = " "

            if (
                (not self.prefix)
                or (self.prefix.isspace())
                or (self.prefix in [Undefined, "<undefined>"])
            ):
                self.prefix = " "

            path, _ = os.path.split(self.in_file)

            if (
                self.suffix == " "
                and self.prefix == " "
                and path == self.output_directory
            ):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle(
                    "mia_processes - " "Sanitize brick Warning!"
                )
                msg.setText(
                    "Suffix and prefix input parameters are not "
                    "defined or consist only of one or more white "
                    "spaces.\nThe {0} input parameter will be "
                    "overwritten ...\n Yes or "
                    "Abort?".format(fileName)
                )
                msg.setStandardButtons(QMessageBox.Yes | QMessageBox.Abort)
                retval = msg.exec_()

                if retval != QMessageBox.Abort:
                    print(
                        "\nSanitize brick warning: the out_file output "
                        "parameter is the same as the in_files input "
                        "parameter (suffix and prefix are not defined):"
                        "\n{0} will be overwrited ...".format(fileName)
                    )

                else:
                    print(
                        "\nSanitize brick: Aborted. Please check your input "
                        "parameters ..."
                    )
                    return

            file = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + fileName
                    + self.suffix.strip()
                    + "."
                    + in_ext
                ),
            )

            if file:
                self.outputs["out_file"] = file
                # FIXME: In the latest version of mia, indexing of the
                #        database with particular tags defined in the
                #        processes is done only at the end of the
                #        initialisation of the whole pipeline. So we
                #        cannot use the value of these tags in other
                #        processes of the pipeline at the time of
                #        initialisation (see populse_mia #290). Until
                #        better we use a quick and dirty hack with the
                #        set_dbFieldValue() function !

                tag_to_add = dict()
                tag_to_add["name"] = "RepetitionTime"
                tag_to_add["field_type"] = "float"
                tag_to_add["description"] = (
                    "The period of time "
                    "in msec between the "
                    "beginning of a pulse "
                    "sequence and the "
                    "beginning of the "
                    "succeeding pulse "
                    "sequence"
                )
                tag_to_add["visibility"] = True
                tag_to_add["origin"] = "user"
                tag_to_add["unit"] = "ms"
                tag_to_add["default_value"] = None
                tag_to_add["value"] = get_dbFieldValue(
                    self.project, self.in_file, "RepetitionTime"
                )

                if tag_to_add["value"] is not None:
                    set_dbFieldValue(self.project, file, tag_to_add)

                else:
                    print(
                        "Sanitize:\n the 'RepetitionTime' tag could "
                        "not be added to the database for the '{}' "
                        "parameter. This can lead to a subsequent issue "
                        "during initialization!!\n".format(file)
                    )

            else:
                print(
                    "- Sanitize brick: There was no output file deducted "
                    "during initialisation. "
                    "Please check the input parameters...!"
                )

            # tags inheritance (optional)
            if self.outputs["out_file"]:
                self.inheritance_dict[self.outputs["out_file"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Sanitize, self).run_process_mia()

        file_name = self.in_file

        # Image processing
        try:
            img = nib.load(file_name)

        except (
            nib.filebasedimages.ImageFileError,
            FileNotFoundError,
            TypeError,
        ) as e:
            print(
                "\nSanitize brick: Error with files to mask, during the "
                "run: ",
                e,
            )
            img = None

        if img is not None:
            # Retrieve xform codes
            sform_code = int(img.header._structarr["sform_code"])
            qform_code = int(img.header._structarr["qform_code"])

            # Check qform is valid
            valid_qform = False
            try:
                img.get_qform()
                valid_qform = True
            except ValueError:
                pass

            # Matching affines
            matching_affines = valid_qform and np.allclose(
                img.get_qform(), img.get_sform()
            )

            save_file = False

            # Both match, qform valid (implicit with match), codes OK
            # -> do nothing, empty report
            if matching_affines and qform_code > 0 and sform_code > 0:
                save_file = True

            # Row 2:
            elif valid_qform and qform_code > 0:
                img.set_sform(img.get_qform(), qform_code)
                save_file = True
                print(
                    "\nSanitize brick: Note on orientation:"
                    "\nThe sform has been copied from qform."
                )

            # Rows 3-4:
            # Note: if qform is not valid, matching_affines is False
            elif sform_code > 0 and (not matching_affines or qform_code == 0):
                img.set_qform(img.get_sform(), sform_code)
                save_file = True
                print(
                    "\nSanitize brick: Note on orientation:"
                    "\nThe qform has been copied from sform."
                )

                if not valid_qform and qform_code > 0:
                    print(
                        "\nSanitize brick WARNING: Invalid qform "
                        "information.\nThe qform matrix found in the file "
                        "header is invalid. The qform has been copied from "
                        "sform. Checking the original qform information from "
                        "the data produced by the scanner is advised."
                    )

            # Rows 5-6:
            else:
                affine = img.affine
                img.set_sform(affine, nib.nifti1.xform_codes["scanner"])
                img.set_qform(affine, nib.nifti1.xform_codes["scanner"])
                save_file = True
                print(
                    "\nSanitize brick WARNING: Missing orientation "
                    "information.\nOrientation information could not be "
                    "retrieved from the image header. The qform and sform "
                    "matrices have been set to a default LAS-oriented "
                    "affine. Analyses of this dataset MAY BE INVALID!"
                )

            if (
                self.max_32bit and np.dtype(img.get_data_dtype()).itemsize > 4
            ) or self.n_volumes_to_discard:
                # force float32 only if 64 bit dtype is detected
                if (
                    self.max_32bit
                    and np.dtype(img.get_data_dtype()).itemsize > 4
                ):
                    in_data = img.get_fdata(dtype=np.float32)
                else:
                    in_data = img.dataobj

                # fmt: off
                img = nib.Nifti1Image(
                    in_data[:, :, :, self.n_volumes_to_discard:],
                    img.affine,
                    img.header,
                )
                # fmt: on
                save_file = True

            if len(img.header.extensions) != 0:
                img.header.extensions.clear()
                save_file = True

            # Store new file
            if save_file:
                # Image save
                _, file_name = os.path.split(file_name)
                file_name_no_ext, file_extension = os.path.splitext(file_name)
                if file_extension == ".gz":
                    (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                        file_name_no_ext
                    )
                    if file_extension_2 == ".nii":
                        file_name_no_ext = file_name_no_ext_2
                        file_extension = ".nii.gz"

                file_out = os.path.join(
                    self.output_directory,
                    (
                        self.prefix.strip()
                        + file_name_no_ext
                        + self.suffix.strip()
                        + file_extension
                    ),
                )
                nib.save(img, file_out)


class TSNR(ProcessMIA):
    """
    *Computes the time-course SNR for a time series*

    Please, see the complete documentation for the
    `TSNR brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/TSNR.html>`_

    adapted from:
    https://github.com/nipy/nipype/blob/f662acfce8def4717e0c3414618f3a5de5913b31/nipype/algorithms/confounds.py#L899

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(TSNR, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_file_desc = "An existing path file."
        prefix_tsnr_desc = "Prefix of the tsnr image (a string)."
        suffix_tsnr_desc = "Suffix of the tsnr image (a string)."
        prefix_stddev_desc = "Prefix of the tsnr image (a string)."
        suffix_stddev_desc = "Suffix of the tsnr image (a string)."

        # Outputs description
        out_tsnr_file_desc = (
            "Path of the scan after masking "
            "(a pathlike object or string representing a file)."
        )

        out_stddev_file_desc = (
            "Path of the scan after masking "
            "(a pathlike object or string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        self.add_trait(
            "prefix_tsnr",
            traits.String(
                "", output=False, optional=True, desc=prefix_tsnr_desc
            ),
        )

        self.add_trait(
            "suffix_tsnr",
            traits.String(
                "_tsnr", output=False, optional=True, desc=suffix_tsnr_desc
            ),
        )

        self.add_trait(
            "prefix_stddev",
            traits.String(
                "", output=False, optional=True, desc=prefix_stddev_desc
            ),
        )

        self.add_trait(
            "suffix_stddev",
            traits.String(
                "_stddev", output=False, optional=True, desc=suffix_stddev_desc
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_tsnr_file", File(output=True, desc=out_tsnr_file_desc)
        )

        self.add_trait(
            "out_stddev_file",
            File(output=True, optional=True, desc=out_stddev_file_desc),
        )

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(TSNR, self).list_outputs()

        # Outputs definition
        if self.in_file:
            valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

            if not valid_ext:
                print(
                    "\nTSNR brick: The input image format is not "
                    "recognized...!"
                )
                return

            if (
                (not self.suffix_tsnr)
                or (self.suffix_tsnr.isspace())
                or self.suffix_tsnr in [Undefined, "<undefined>"]
            ):
                self.suffix_tsnr = " "

            if (
                (not self.prefix_tsnr)
                or (self.prefix_tsnr.isspace())
                or (self.prefix_tsnr in [Undefined, "<undefined>"])
            ):
                self.prefix_tsnr = " "

            file_tsnr = ""

            path, _ = os.path.split(self.in_file)

            if (
                self.suffix_tsnr == " "
                and self.prefix_tsnr == " "
                and path == self.output_directory
            ):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle("mia_processes - TSNR brick Warning!")
                msg.setText(
                    "suffix_tsnr and prefix_tsnr input parameters are not "
                    "defined or consist only of one or more white "
                    "spaces.\nThe {0} input parameter will be "
                    "overwritten ...\n Yes or "
                    "Abort?".format(fileName)
                )
                msg.setStandardButtons(QMessageBox.Yes | QMessageBox.Abort)
                retval = msg.exec_()

                if retval != QMessageBox.Abort:
                    print(
                        "\nTSNR brick warning: the out_tsnr_file output "
                        "parameter is the same as the in_file input "
                        "parameter (suffix and prefix are not defined):"
                        "\n{0} will be overwrited ...".format(fileName)
                    )

                else:
                    print(
                        "\nTSNR brick: Aborted. Please check your input "
                        "parameters ..."
                    )
                    return

            else:
                file_tsnr = os.path.join(
                    self.output_directory,
                    (
                        self.prefix_tsnr.strip()
                        + fileName
                        + self.suffix_tsnr.strip()
                        + "."
                        + in_ext
                    ),
                )

            if (
                (not self.suffix_stddev)
                or (self.suffix_stddev.isspace())
                or self.suffix_stddev in [Undefined, "<undefined>"]
            ):
                self.suffix_stddev = " "

            if (
                (not self.prefix_stddev)
                or (self.prefix_stddev.isspace())
                or (self.prefix_stddev in [Undefined, "<undefined>"])
            ):
                self.prefix_stddev = " "

            file_stddev = ""

            if (
                self.suffix_stddev == " "
                and self.prefix_stddev == " "
                and path == self.output_directory
            ):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle("mia_processes - TSNR brick Warning!")
                msg.setText(
                    "suffix_stddev and prefix_stddev input parameters are not "
                    "defined or consist only of one or more white "
                    "spaces.\nThe {0} input parameter will be "
                    "overwritten ...\n Yes or "
                    "Abort?".format(fileName)
                )
                msg.setStandardButtons(QMessageBox.Yes | QMessageBox.Abort)
                retval = msg.exec_()

                if retval != QMessageBox.Abort:
                    print(
                        "\nTSNR brick warning: the out_stddev_file output "
                        "parameter is the same as the in_files input "
                        "parameter (suffix and prefix are not defined):"
                        "\n{0} will be overwrited ...".format(fileName)
                    )
                else:
                    print(
                        "\nTSNR brick Aborted. Please check your input "
                        "parameters ..."
                    )
                    return

            else:
                file_stddev = os.path.join(
                    self.output_directory,
                    (
                        self.prefix_stddev.strip()
                        + fileName
                        + self.suffix_stddev.strip()
                        + "."
                        + in_ext
                    ),
                )

            self.outputs["out_tsnr_file"] = file_tsnr
            self.outputs["out_stddev_file"] = file_stddev

            # tags inheritance (optional)
            if self.outputs["out_tsnr_file"]:
                self.inheritance_dict[
                    self.outputs["out_tsnr_file"]
                ] = self.in_file
            if self.outputs["out_stddev_file"]:
                self.inheritance_dict[
                    self.outputs["out_stddev_file"]
                ] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(TSNR, self).run_process_mia()

        file_name = self.in_file

        try:
            img = nib.load(file_name)
        except (
            nib.filebasedimages.ImageFileError,
            FileNotFoundError,
            TypeError,
        ) as e:
            print(
                "\nTSNR brick: Error with files to mask, during the " "run: ",
                e,
            )
            img = None

        if img is not None:
            header = img.header.copy()
            data = img.get_fdata(dtype=np.float32).reshape(
                img.shape[:3] + (-1,)
            )
            data = np.nan_to_num(data)

            if data.dtype.kind == "i":
                header.set_data_dtype(np.float32)
                data = data.astype(np.float32)

            meanimg = np.mean(data, axis=3)
            stddevimg = np.std(data, axis=3)
            tsnr = np.zeros_like(meanimg)
            stddevimg_nonzero = stddevimg > 1.0e-3
            tsnr[stddevimg_nonzero] = (
                meanimg[stddevimg_nonzero] / stddevimg[stddevimg_nonzero]
            )

            # Image save
            _, file_name = os.path.split(file_name)
            file_name_no_ext, file_extension = os.path.splitext(file_name)
            if file_extension == ".gz":
                (file_name_no_ext_2, file_extension_2) = os.path.splitext(
                    file_name_no_ext
                )
                if file_extension_2 == ".nii":
                    file_name_no_ext = file_name_no_ext_2
                    file_extension = ".nii.gz"

            tsnr_file_out = os.path.join(
                self.output_directory,
                (
                    self.prefix_tsnr.strip()
                    + file_name_no_ext
                    + self.suffix_tsnr.strip()
                    + file_extension
                ),
            )
            tsnr_img = nib.Nifti1Image(tsnr, img.affine, header)
            nib.save(tsnr_img, tsnr_file_out)

            stddev_file_out = os.path.join(
                self.output_directory,
                (
                    self.prefix_stddev.strip()
                    + file_name_no_ext
                    + self.suffix_stddev.strip()
                    + file_extension
                ),
            )
            stddev_img = nib.Nifti1Image(stddevimg, img.affine, header)
            nib.save(stddev_img, stddev_file_out)


class TemplateFromTemplateFlow(ProcessMIA):
    """
    *Get template image from templateflow*

    Please, see the complete documentation for the
    `TemplateFromTemplateFlow brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/TemplateFromTemplateFlow.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(TemplateFromTemplateFlow, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_template_desc = "Template Name (a String)."
        resolution_desc = "Resolution of the template  (Int or None)"
        suffix_desc = "BIDS suffix (String or None)"
        atlas_desc = "Name of a particular atlas (String or None)"
        desc_desc = "Description field (String or None)"
        label_desc = "Label field (String or None)"

        # Outputs description
        template_desc = (
            "Path of the template (a pathlike object "
            "or string representing a file)."
        )

        # Inputs traits
        self.add_trait(
            "atlas",
            traits.String("", output=False, optional=True, desc=atlas_desc),
        )

        self.add_trait(
            "desc",
            traits.String("", output=False, optional=True, desc=desc_desc),
        )

        self.add_trait(
            "in_template",
            traits.String(
                default="MNI152NLin2009cAsym",
                output=False,
                optional=True,
                desc=in_template_desc,
            ),
        )

        self.add_trait(
            "label",
            traits.String("", output=False, optional=True, desc=label_desc),
        )

        self.add_trait(
            "resolution",
            traits.Int(2, output=False, optional=True, desc=resolution_desc),
        )

        self.add_trait(
            "suffix",
            traits.String("", output=False, optional=True, desc=suffix_desc),
        )

        # Outputs traits
        self.add_trait("template", File(output=True, desc=template_desc))

        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(TemplateFromTemplateFlow, self).list_outputs()

        if not self.suffix:
            suffix = None
        else:
            suffix = self.suffix

        if not self.atlas:
            atlas = None
        else:
            atlas = self.atlas

        if not self.desc:
            desc = None
        else:
            desc = self.desc

        if not self.label:
            label = None
        else:
            label = self.label

        template_spec = {
            "resolution": self.resolution,
            "suffix": suffix,
            "atlas": atlas,
            "desc": desc,
            "label": label,
        }

        tpl_target_path = get_template(self.in_template, **template_spec)
        if not tpl_target_path:
            print(
                "\nTemplateFromTemplateFlow brick: Could not find template "
                "'{0}' with specs={1}. Please revise your template "
                "argument.".format(self.in_template, template_spec)
            )
            return

        if isinstance(tpl_target_path, list):
            print(
                "\nTemplateFromTemplateFlow brick: The available template "
                "modifiers ({0}) did not select a unique template (got "
                "'{1}'). Please revise your template "
                "argument.".format(self.in_template, template_spec)
            )
            return

        self.outputs["template"] = str(tpl_target_path)

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(TemplateFromTemplateFlow, self).run_process_mia()


class Threshold(ProcessMIA):
    """
    *Makes a binary mask image at a given threshold*

    Please, see the complete documentation for the
    `Threshold brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/other/Threshold.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Threshold, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_files_desc = (
            "A list of items with string elements corresponding "
            "to existing path files."
        )
        threshold_desc = (
            "Value for the applied threshold (a float between 0 " "and 1)."
        )
        GM_filtering_desc = (
            "Filtering (based on the keyword c1) to keep only "
            "grey matter images (a boolean)"
        )
        suffix_desc = "Suffix of the output image (a string)."
        prefix_desc = "Prefix of the output image (a string)."

        # Outputs description
        out_files_desc = (
            "Path of the scan after application of the threshold "
            "(a pathlike object or string representing a file, or "
            "a list of pathlike objects or strings representing a "
            "file)."
        )

        # Inputs traits
        self.add_trait(
            "in_files",
            InputMultiPath(
                traits.Either(File(), traits.List(File())),
                output=False,
                desc=in_files_desc,
            ),
        )

        self.add_trait(
            "GM_filtering",
            traits.Bool(
                default_value=True,
                output=False,
                optional=True,
                desc=GM_filtering_desc,
            ),
        )

        self.add_trait(
            "prefix",
            traits.String("", output=False, optional=True, desc=prefix_desc),
        )

        self.add_trait(
            "suffix",
            traits.String(
                "_002", output=False, optional=True, desc=suffix_desc
            ),
        )

        self.add_trait(
            "threshold",
            traits.Range(
                value=0.3,
                low=0.0,
                high=1.0,
                output=False,
                optional=True,
                desc=threshold_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_files",
            OutputMultiPath(File(), output=True, desc=out_files_desc),
        )

        self.init_default_traits()

    def _namesFilter(self):
        """Filtering of in_files input parameter in order to keep only GM.

        :returns: GM images (containing "c1" in the first 5 characters)
        """
        files = []

        for file_name in self.in_files:
            if isinstance(file_name, (list, TraitListObject)):
                file_name = file_name[0]

            # Take the first 5 characters in the case if the GM was processed
            # ex. normalisation = wc1, normalisation then smooth = swc1 ...
            # This is not the cleaner way ... for example in the case of a
            # file name including `c1` in the first 5 characters for another
            # reason (`c1` in PatientName, etc.).
            # FiXME: what can we do to be cleaner?
            if "c1" in os.path.basename(os.path.normpath(file_name))[:5]:
                files.append(file_name)

        return files

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Threshold, self).list_outputs()

        # Outputs definition
        if self.in_files != Undefined and (
            self.threshold and self.threshold != Undefined
        ):
            if self.GM_filtering is True:
                files_name = self._namesFilter()

            else:
                files_name = self.in_files

            if (
                (not self.suffix)
                or (self.suffix.isspace())
                or self.suffix in [Undefined, "<undefined>"]
            ):
                self.suffix = " "

            if (
                (not self.prefix)
                or (self.prefix.isspace())
                or (self.prefix in [Undefined, "<undefined>"])
            ):
                self.prefix = " "

            files = []
            flag = True  # If False, suf/pref check will not be performed later
            retval = None

            for file_name1 in files_name:
                valid_ext, in_ext, fileName = checkFileExt(file_name1, EXT)

                if not valid_ext:
                    print(
                        "\nThreshold brick: The input image format is not "
                        "recognized...!"
                    )
                    return
                path, _ = os.path.split(file_name1)

                if (
                    self.suffix == " "
                    and self.prefix == " "
                    and path == self.output_directory
                    and flag is True
                ):
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Warning)
                    msg.setWindowTitle(
                        "mia_processes - " "Threshold brick Warning!"
                    )
                    msg.setText(
                        "Suffix and prefix input parameters are not "
                        "defined or consist only of one or more white "
                        "spaces.\nThe {0} input parameter will be "
                        "overwritten ...\n Yes or "
                        "Abort?".format(file_name1)
                    )
                    msg.setStandardButtons(
                        QMessageBox.Yes
                        | QMessageBox.YesToAll
                        | QMessageBox.Abort
                    )
                    retval = msg.exec_()

                    if retval != QMessageBox.Abort:
                        print(
                            "\nThreshold brick warning: the out_files output "
                            "parameter is the same as the in_files input "
                            "parameter (suffix and prefix are not defined):"
                            "\n{0} will be overwrited ...".format(file_name1)
                        )

                        if retval == QMessageBox.YesToAll:
                            flag = False
                            print(
                                "\nThreshold brick: YesToAll selected: end of "
                                "overwrite checks on input images ..."
                            )
                    else:
                        files_name = []
                        print(
                            "\nThreshold brick Aborted. Please check your "
                            "input parameters ..."
                        )
                        return

                files.append(
                    os.path.join(
                        self.output_directory,
                        (
                            self.prefix.strip()
                            + fileName
                            + self.suffix.strip()
                            + "."
                            + in_ext
                        ),
                    )
                )

            if files:
                self.outputs["out_files"] = files
            else:
                print(
                    "- Threshold brick: There was no output file deducted "
                    "during initialisation. Please check the input "
                    "parameters...!"
                )

        # tags inheritance (optional)
        if self.outputs:
            for key, val in self.outputs.items():
                if key == "out_files":
                    for in_val, out_val in zip(files_name, val):
                        _, fileOval = os.path.split(out_val)
                        fileOval_no_ext, _ = os.path.splitext(fileOval)
                        _, fileIval = os.path.split(in_val)
                        fileIval_no_ext, _ = os.path.splitext(fileIval)

                        if (self.prefix) and (
                            fileOval_no_ext.startswith(self.prefix)
                        ):
                            # fmt: off
                            fileOval_no_ext = fileOval_no_ext[
                                len(self.prefix):]
                            # fmt: on

                        if (self.suffix) and (
                            fileOval_no_ext.endswith(self.suffix)
                        ):
                            fileOval_no_ext = fileOval_no_ext[
                                : -len(self.suffix)
                            ]

                        if fileOval_no_ext == fileIval_no_ext:
                            self.inheritance_dict[out_val] = in_val
                            # FIXME: In the latest version of mia, indexing of
                            #        the database with particular tags defined
                            #        in the processes is done only at the end
                            #        of the initialisation of the whole
                            #        pipeline. So we cannot use the value of
                            #        these tags in other processes of the
                            #        pipeline at the time of initialisation
                            #        (see populse_mia #290). Until better we
                            #        use a quick and dirty hack with the
                            #        set_dbFieldValue() function !
                            tag_to_add = dict()
                            tag_to_add["name"] = "PatientName"
                            tag_to_add["field_type"] = "string"
                            tag_to_add["description"] = ""
                            tag_to_add["visibility"] = True
                            tag_to_add["origin"] = "user"
                            tag_to_add["unit"] = None
                            tag_to_add["default_value"] = None
                            tag_to_add["value"] = get_dbFieldValue(
                                self.project, in_val, "PatientName"
                            )

                            if tag_to_add["value"] is not None:
                                set_dbFieldValue(
                                    self.project, out_val, tag_to_add
                                )

                            else:
                                print(
                                    "\nThreshold brick:\nThe 'PatientName'"
                                    " tag could not be added to the "
                                    "database for the '{}' parameter. This "
                                    "can lead to a subsequent issue during "
                                    "initialization!!\n".format(out_val)
                                )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Threshold, self).run_process_mia()

        if self.GM_filtering is True:
            files_name = self._namesFilter()

        else:
            files_name = self.in_files

        for file_name in files_name:
            # Image processing
            img_final = threshold(file_name, self.threshold)

            # Image save
            _, file_name = os.path.split(file_name)
            file_name_no_ext, file_extension = os.path.splitext(file_name)
            out_file = os.path.join(
                self.output_directory,
                (
                    self.prefix.strip()
                    + file_name_no_ext
                    + self.suffix.strip()
                    + file_extension
                ),
            )
            nib.save(img_final, out_file)


def artifact_mask(imdata, airdata, distance, zscore=10.0):
    """Computes a mask of artifacts found in the air region"""

    if not np.issubdtype(airdata.dtype, np.integer):
        airdata[airdata < 0.95] = 0
        airdata[airdata > 0.0] = 1

    bg_img = imdata * airdata
    if np.sum((bg_img > 0).astype(np.uint8)) < 100:
        return np.zeros_like(airdata)

    # Find the background threshold (the most frequently occurring value
    # excluding 0)
    bg_location = np.median(bg_img[bg_img > 0])
    bg_spread = mad(bg_img[bg_img > 0])
    bg_img[bg_img > 0] -= bg_location
    # TODO: This is generating a RuntimeWarning: divide by zero with MRIQC anat
    bg_img[bg_img > 0] /= bg_spread

    # Apply this threshold to the background voxels to identify voxels
    # contributing artifacts.
    qi1_img = np.zeros_like(bg_img)
    qi1_img[bg_img > zscore] = 1
    qi1_img[distance < 0.10] = 0

    # Create a structural element to be used in an opening operation.
    struc = sim.generate_binary_structure(3, 1)
    qi1_img = sim.binary_opening(qi1_img, struc).astype(np.uint8)
    qi1_img[airdata <= 0] = 0

    return qi1_img


def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.
    :param nparray points: an numobservations by numdimensions numpy array
        of observations
    :param float thresh: the modified z-score to use as a threshold.
        Observations with a modified z-score (based on the median absolute
        deviation) greater than this value will be classified as outliers.
    :return: A bolean mask, of size numobservations-length array.
    .. note:: References
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """
    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median) ** 2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    timepoints_to_discard = 0
    for i in range(len(modified_z_score)):
        if modified_z_score[i] <= thresh:
            break
        else:
            timepoints_to_discard += 1

    return timepoints_to_discard


def threshold(file_name, thresh):
    """
    Basic method for image thresholding

    :param file_name: Image to be thresholded
    :param thresh: Threshold value (a float between 0 and 1)
    :returns: Image after thresholding
    """
    img = nib.load(file_name)
    img_data = img.get_fdata()
    _max = img_data.max()
    img_thresh = (img_data > (_max * thresh)).astype(float)
    img_final = nib.Nifti1Image(img_thresh, img.affine, img.header)
    return img_final
