# -*- coding: utf-8 -*-

"""The freesurfer preprocess library of the mia_processes package.

The purpose of this module is to customise the main freesurfer
preprocessing bricks provided by nipype and to correct some things
that do not work directly in populse_mia.

:Contains:
    :Class:
        - Binarize
        - SynthStrip
        - SynthStripMriqc

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

import os

# capsul import
import capsul

# Other import
import nibabel as nb
import numpy as np
import scipy
import torch
from capsul.in_context import freesurfer
from nipype.interfaces.base import File, Undefined
from nitransforms.linear import Affine

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA
from traits.api import Bool, Either, Enum, Float, Int, List, String

# mia_processes import
from mia_processes.utils import checkFileExt

EXT = {"NIFTI_GZ": "nii.gz", "NIFTI": "nii", "MGZ": "mgz"}


class Binarize(ProcessMIA):
    """
    *Binarize a volume (or volume-encoded surface file) using \
FreeSurfer mri_binarize.

    Binarization can be done based on threshold or on matched values.

    Please, see the complete documentation for the `Binarize brick in
    the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/freesurfer/Binarize.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Binarize, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["freesurfer", "nipype"]

        # Mandatory inputs description
        in_file_desc = (
            "Input file (a pathlike object or string " "representing a file)."
        )
        # Optional inputs with default value description
        abs_desc = "take abs of invol first (ie, make unsigned)"
        bin_col_num_desc = "set binarized voxel value to its column number"
        get_count_file_desc = (
            "save number of hits in ascii file" "(hits, ntotvox, pct)"
        )
        invert_desc = "set binval=0, binvalnot=1"
        max_desc = "Maximum voxel threshold(float)."
        min_desc = "Minimum voxel threshold(float)."
        output_type_desc = (
            "Typecodes of the output image formats (one "
            "of NIFTI, MGZ, NIFTI_GZ)."
        )
        out_suffix_desc = "Suffix of the output image (a string)."
        zero_edges_desc = "zero the edge voxels"
        zero_slice_edge_desc = "zero the edge slice voxels"
        # Optional inputs description
        bin_val_desc = "set vox outside range to val (default is 0)"
        bin_val_not_desc = "set vox outside range to val (default is 0)"
        dilate_desc = "niters: dilate binarization in 3D"
        erode_desc = "nerode: erode binarization in 3D " "(after any dilation)"
        erode2d_desc = (
            "nerode2d: erode binarization in 2D " "(after any 3D erosion)"
        )
        frame_no_desc = "use 0-based frame of input (default is 0)"
        match_desc = "Match instead of threshold"

        rmax_desc = "Compute max based on rmax*globalmean."
        rmin_desc = "Compute min based on rmin*globalmean."
        # Outputs description
        count_file_desc = (
            "File that contains number of hits" "(hits, ntotvox, pct)"
        )
        out_file_desc = (
            "The binanized file (a pathlike object or a "
            "string representing a file)."
        )

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits
        self.add_trait(
            "abs",
            Bool(default=False, output=False, optional=True, desc=abs_desc),
        )

        self.add_trait(
            "bin_col_num",
            Bool(
                default=False,
                output=False,
                optional=True,
                desc=bin_col_num_desc,
            ),
        )

        self.add_trait(
            "get_count_file",
            Bool(
                False,
                default=False,
                output=False,
                optional=True,
                desc=get_count_file_desc,
            ),
        )

        self.add_trait(
            "invert",
            Bool(default=False, output=False, optional=True, desc=invert_desc),
        )
        self.add_trait(
            "max",
            Either(
                Float(),
                Undefined,
                default=100.0,
                output=False,
                optional=True,
                desc=max_desc,
            ),
        )
        self.add_trait(
            "min",
            Either(
                Float(),
                Undefined,
                default=0.0,
                output=False,
                optional=True,
                desc=min_desc,
            ),
        )

        self.add_trait(
            "output_type",
            Enum(
                "NIFTI",
                "NIFTI_GZ",
                "MGZ",
                output=False,
                optional=True,
                desc=output_type_desc,
            ),
        )

        self.add_trait(
            "out_suffix",
            String(
                "_thresh", output=False, optional=True, desc=out_suffix_desc
            ),
        )

        self.add_trait(
            "zero_edges",
            Bool(
                default=False,
                output=False,
                optional=True,
                desc=zero_edges_desc,
            ),
        )

        self.add_trait(
            "zero_slice_edge",
            Bool(
                default=False,
                output=False,
                optional=True,
                desc=zero_slice_edge_desc,
            ),
        )

        # Optional inputs traits
        self.add_trait(
            "bin_val",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=bin_val_desc,
            ),
        )

        self.add_trait(
            "bin_val_not",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=bin_val_not_desc,
            ),
        )

        self.add_trait(
            "dilate",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=dilate_desc,
            ),
        )

        self.add_trait(
            "erode",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=erode_desc,
            ),
        )

        self.add_trait(
            "erode2d",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=erode2d_desc,
            ),
        )
        self.add_trait(
            "frame_no",
            Either(
                Undefined,
                Int(),
                default=Undefined,
                output=False,
                optional=True,
                desc=frame_no_desc,
            ),
        )

        self.add_trait(
            "match",
            Either(
                Undefined,
                List(Int()),
                default=Undefined,
                output=False,
                optional=True,
                desc=match_desc,
            ),
        )

        self.add_trait(
            "rmax",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=rmax_desc,
            ),
        )

        self.add_trait(
            "rmin",
            Either(
                Undefined,
                Float(),
                default=Undefined,
                output=False,
                optional=True,
                desc=rmin_desc,
            ),
        )
        # Outputs traits
        self.add_trait(
            "count_file",
            File(output=True, optional=True, desc=count_file_desc),
        )

        self.add_trait("out_file", File(output=True, desc=out_file_desc))

        self.init_default_traits()
        self.init_process("nipype.interfaces.freesurfer.Binarize")

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
        super(Binarize, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if (
            self.min != Undefined
            or self.max != Undefined
            or self.rmin != Undefined
            or self.rmax != Undefined
        ) and self.match != Undefined:
            print(
                '\nInitialisation failed. "match" parameter can not be used '
                'with "min" and/or "max" parameters'
                ' Please, define only "min"/"max" paremeters or "match"'
                "parameters (set the other as Undefined) ...!"
            )
            return

        if self.in_file:
            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is" " not recognized...!")
                    return
                else:
                    output_type = self.output_type
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        os.path.split(self.in_file)[1].replace(
                            "." + in_ext,
                            self.out_suffix + "." + EXT[output_type],
                        ),
                    )

                    self.outputs["count_file"] = os.path.join(
                        self.output_directory,
                        os.path.split(self.in_file)[1].replace(
                            "." + in_ext, self.out_suffix + "_count.txt"
                        ),
                    )

                    self.inheritance_dict[
                        self.outputs["out_file"]
                    ] = self.in_file
            else:
                print("No output_directory was found...!\n")
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Binarize, self).run_process_mia()

        # madatory inputs / outputs
        self.process.in_file = self.in_file
        self.process.binary_file = self.out_file
        self.process.out_type = EXT[self.output_type]

        # optionnal inputs
        self.process.min = self.min
        self.process.max = self.max
        self.process.rmin = self.rmin
        self.process.rmax = self.rmax
        self.process.match = self.match
        if self.get_count_file:
            self.process.count_file = self.count_file
        self.process.bin_val = self.bin_val
        self.process.bin_val_not = self.bin_val_not
        self.process.frame_no = self.frame_no
        self.process.dilate = self.dilate
        self.process.erode = self.erode
        self.process.erode2d = self.erode2d
        self.process.invert = self.invert
        self.process.abs = self.abs
        self.process.bin_col_num = self.bin_col_num
        self.process.zero_edges = self.zero_edges
        self.process.zero_slice_edge = self.zero_slice_edge

        return self.process.run(configuration_dict={})


class SynthStrip(ProcessMIA):
    """
    *Skull stripping using SynthStrip*

    Please, see the complete documentation for the 'SynthStrip brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/freesurfer/SynthStrip.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(SynthStrip, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["freesurfer", "nipype"]

        # Mandatory inputs description
        in_file_desc = "Input image to be brain extracted"
        # Optional inputs with default value description
        border_mm_desc = "Mask border threshold in mm"
        no_csf_desc = "Exclude CSF from brain border"
        output_type_desc = (
            "Typecodes of the output image formats (one "
            "of NIFTI, MGZ, NIFTI_GZ)."
        )
        # Outputs description
        out_file_desc = "Brain-extracted path"
        out_mask_desc = "Brain mask path"

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits
        self.add_trait(
            "border_mm",
            Int(1, output=False, optional=True, desc=border_mm_desc),
        )

        self.add_trait(
            "no_csf",
            Bool(False, output=False, optional=True, desc=no_csf_desc),
        )

        self.add_trait(
            "output_type",
            Enum(
                "NIFTI",
                "NIFTI_GZ",
                "MGZ",
                output=False,
                optional=True,
                desc=output_type_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=True, desc=out_file_desc)
        )

        self.add_trait("out_mask", File(output=True, desc=out_mask_desc))

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
        super(SynthStrip, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is" " not recognized...!")
                    return
                else:
                    output_type = self.output_type
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        os.path.split(self.in_file)[1].replace(
                            "." + in_ext, "_desc-brain." + EXT[output_type]
                        ),
                    )

                    self.outputs["out_mask"] = os.path.join(
                        self.output_directory,
                        os.path.split(self.in_file)[1].replace(
                            "." + in_ext,
                            "_desc-brain_mask." + EXT[output_type],
                        ),
                    )

            else:
                print("No output_directory was found...!\n")
                return

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file
            self.inheritance_dict[self.outputs["out_mask"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(SynthStrip, self).run_process_mia()

        # default input
        fconf = capsul.engine.configurations.get(
            "capsul.engine.module.freesurfer"
        )
        model_path = os.path.join(
            os.path.dirname(fconf["setup"]), "models", "synthstrip.1.pt"
        )

        cmd = [
            "mri_synthstrip",
            "-i",
            self.in_file,
            "-o",
            self.out_file,
            "-m",
            self.out_mask,
            "-b",
            str(self.border_mm),
            "--model",
            model_path,
        ]

        if self.no_csf:
            cmd += ["--no-csf"]

        return freesurfer.freesurfer_call(cmd)


class SynthStripMriqc(ProcessMIA):
    """
    *Skull stripping using SynthStrip*

    STATEMENT OF CHANGES: This class is derived from the work done by
    MRIQC 22.06 and  FreeSurfer 1.0. The original file for this work derives
    from is found at
    https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/synthstrip/cli.py
    and
    https://github.com/nipreps/mriqc/blob/22.0.6/mriqc/synthstrip/model.py
    and
    https://github.com/freesurfer/freesurfer/blob/2995ded957961a7f3704de57eee88eb6cc30d52d/mri_synthstrip/mri_synthstrip

    Please, see the complete documentation for the 'SynthStripMriqc brick
    in the populse.mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/preprocess/freesurfer/SynthStripMriqc.html>`_

    """

    class StripModel(torch.nn.Module):
        """blabla"""

        def __init__(
            self,
            nb_features=16,
            nb_levels=7,
            feat_mult=2,
            max_features=64,
            nb_conv_per_level=2,
            max_pool=2,
            return_mask=False,
        ):
            """blabla"""
            super().__init__()
            # dimensionality
            ndims = 3

            # build feature list automatically
            if isinstance(nb_features, int):
                if nb_levels is None:
                    raise ValueError(
                        "SynthStripMriqc brick: must provide unet nb_levels "
                        "if nb_features is an integer"
                    )
                feats = np.round(
                    nb_features * feat_mult ** np.arange(nb_levels)
                ).astype(int)
                feats = np.clip(feats, 1, max_features)
                nb_features = [
                    np.repeat(feats[:-1], nb_conv_per_level),
                    np.repeat(np.flip(feats), nb_conv_per_level),
                ]

            elif nb_levels is not None:
                raise ValueError(
                    "SynthStripMriqc brick: cannot use nb_levels if "
                    "nb_features is not an integer"
                )

            # extract any surplus (full resolution) decoder convolutions
            enc_nf, dec_nf = nb_features
            nb_dec_convs = len(enc_nf)
            final_convs = dec_nf[nb_dec_convs:]
            dec_nf = dec_nf[:nb_dec_convs]
            self.nb_levels = int(nb_dec_convs / nb_conv_per_level) + 1

            if isinstance(max_pool, int):
                max_pool = [max_pool] * self.nb_levels

            # cache downsampling / upsampling operations
            MaxPooling = getattr(torch.nn, "MaxPool%dd" % ndims)
            self.pooling = [MaxPooling(s) for s in max_pool]
            self.upsampling = [
                torch.nn.Upsample(scale_factor=s, mode="nearest")
                for s in max_pool
            ]

            # configure encoder (down-sampling path)
            prev_nf = 1
            encoder_nfs = [prev_nf]
            self.encoder = torch.nn.ModuleList()

            for level in range(self.nb_levels - 1):
                convs = torch.nn.ModuleList()

                for conv in range(nb_conv_per_level):
                    nf = enc_nf[level * nb_conv_per_level + conv]
                    convs.append(SynthStripMriqc.ConvBlock(ndims, prev_nf, nf))
                    prev_nf = nf

                self.encoder.append(convs)
                encoder_nfs.append(prev_nf)

            # configure decoder (up-sampling path)
            encoder_nfs = np.flip(encoder_nfs)
            self.decoder = torch.nn.ModuleList()

            for level in range(self.nb_levels - 1):
                convs = torch.nn.ModuleList()

                for conv in range(nb_conv_per_level):
                    nf = dec_nf[level * nb_conv_per_level + conv]
                    convs.append(SynthStripMriqc.ConvBlock(ndims, prev_nf, nf))
                    prev_nf = nf

                self.decoder.append(convs)

                if level < (self.nb_levels - 1):
                    prev_nf += encoder_nfs[level]

            # now we take care of any remaining convolutions
            self.remaining = torch.nn.ModuleList()

            for num, nf in enumerate(final_convs):
                self.remaining.append(
                    SynthStripMriqc.ConvBlock(ndims, prev_nf, nf)
                )
                prev_nf = nf

            # final convolutions
            if return_mask:
                self.remaining.append(
                    SynthStripMriqc.ConvBlock(
                        ndims, prev_nf, 2, activation=None
                    )
                )
                self.remaining.append(torch.nn.Softmax(dim=1))

            else:
                self.remaining.append(
                    SynthStripMriqc.ConvBlock(
                        ndims, prev_nf, 1, activation=None
                    )
                )

        def forward(self, x):
            """blabla"""
            # encoder forward pass
            x_history = [x]

            for level, convs in enumerate(self.encoder):
                for conv in convs:
                    x = conv(x)

                x_history.append(x)
                x = self.pooling[level](x)

            # decoder forward pass with upsampling and concatenation
            for level, convs in enumerate(self.decoder):
                for conv in convs:
                    x = conv(x)

                if level < (self.nb_levels - 1):
                    x = self.upsampling[level](x)
                    x = torch.cat([x, x_history.pop()], dim=1)

            # remaining convs at full resolution
            for conv in self.remaining:
                x = conv(x)

            return x

    class ConvBlock(torch.nn.Module):
        """
        *Specific convolutional block followed by leakyrelu for unet*

        """

        def __init__(
            self,
            ndims,
            in_channels,
            out_channels,
            stride=1,
            activation="leaky",
        ):
            """blabla"""
            super().__init__()
            Conv = getattr(torch.nn, "Conv%dd" % ndims)
            self.conv = Conv(in_channels, out_channels, 3, stride, 1)

            if activation == "leaky":
                self.activation = torch.nn.LeakyReLU(0.2)

            elif activation is None:
                self.activation = None

            else:
                raise ValueError(f"Unknown activation: {activation}")

        def forward(self, x):
            """blabla"""
            out = self.conv(x)

            if self.activation is not None:
                out = self.activation(out)

            return out

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(SynthStripMriqc, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["freesurfer"]

        # Mandatory inputs description
        in_file_desc = "Input image to be brain extracted (a path to a file)"
        # Optional inputs with default value description
        model_desc = "Alternative model weights (a path to a file)"
        border_mm_desc = "Mask border threshold in mm (integer, default: 1)"
        # no_csf_desc = 'Exclude CSF from brain border'
        output_type_desc = (
            "Typecodes of the output image formats (one "
            "of NIFTI, MGZ, NIFTI_GZ)."
        )
        gpu_desc = "Use the GPU (bool, default: False)"
        # Outputs description
        out_file_desc = "Brain-extracted path (a path to a file)"
        out_mask_desc = "Brain mask path (a path to a file)"

        # Mandatory inputs traits
        self.add_trait(
            "in_file", File(output=False, optional=False, desc=in_file_desc)
        )

        # Optional inputs with default value traits
        self.add_trait(
            "model", File(output=False, optional=True, desc=model_desc)
        )

        self.add_trait(
            "border_mm",
            Int(1, output=False, optional=True, desc=border_mm_desc),
        )

        self.add_trait(
            "output_type",
            Enum(
                "NIFTI",
                "NIFTI_GZ",
                "MGZ",
                output=False,
                optional=True,
                desc=output_type_desc,
            ),
        )

        self.add_trait(
            "gpu", Bool(False, output=False, optional=True, desc=gpu_desc)
        )

        # Outputs traits
        self.add_trait(
            "out_file", File(output=True, optional=True, desc=out_file_desc)
        )

        self.add_trait("out_mask", File(output=True, desc=out_mask_desc))

        self.init_default_traits()

    def conform(self, input_nii):
        """Resample image as SynthStrip likes it."""
        shape = np.array(input_nii.shape[:3])
        affine = input_nii.affine

        # Get corner voxel centers in index coords
        corner_centers_ijk = (
            np.array(
                [
                    (i, j, k)
                    for k in (0, shape[2] - 1)
                    for j in (0, shape[1] - 1)
                    for i in (0, shape[0] - 1)
                ]
            )
            + 0.5
        )

        # Get corner voxel centers in mm
        corners_xyz = (
            affine
            @ np.hstack(
                (corner_centers_ijk, np.ones((len(corner_centers_ijk), 1)))
            ).T
        )

        # Target affine is 1mm voxels in LIA orientation
        target_affine = np.diag([-1.0, 1.0, -1.0, 1.0])[:, (0, 2, 1, 3)]

        # Target shape
        extent = corners_xyz.min(1)[:3], corners_xyz.max(1)[:3]
        target_shape = ((extent[1] - extent[0]) / 1.0 + 0.999).astype(int)

        # SynthStrip likes dimensions be multiple of 64 (192, 256, or 320)
        target_shape = np.clip(
            np.ceil(np.array(target_shape) / 64).astype(int) * 64, 192, 320
        )

        # Ensure shape ordering is LIA too
        target_shape[2], target_shape[1] = target_shape[1:3]

        # Coordinates of center voxel do not change
        input_c = affine @ np.hstack((0.5 * (shape - 1), 1.0))
        target_c = target_affine @ np.hstack((0.5 * (target_shape - 1), 1.0))

        # Rebase the origin of the new, plumb affine
        target_affine[:3, 3] -= target_c[:3] - input_c[:3]

        nii = Affine(
            reference=nb.Nifti1Image(
                np.zeros(target_shape), target_affine, None
            ),
        ).apply(input_nii)
        return nii

    def resample_like(self, image, target, output_dtype=None, cval=0):
        """Resample the input image to be in the target's grid via identity
        transform."""
        return Affine(reference=target).apply(
            image, output_dtype=output_dtype, cval=cval
        )

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
        super(SynthStripMriqc, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.output_directory:
                valid_ext, in_ext, fileName = checkFileExt(self.in_file, EXT)

                if not valid_ext:
                    print("\nThe input image format is" " not recognized...!")
                    return

                else:
                    output_type = self.output_type
                    self.outputs["out_file"] = os.path.join(
                        self.output_directory,
                        os.path.split(self.in_file)[1].replace(
                            "." + in_ext, "_desc-brain." + EXT[output_type]
                        ),
                    )

                    self.outputs["out_mask"] = os.path.join(
                        self.output_directory,
                        os.path.split(self.in_file)[1].replace(
                            "." + in_ext,
                            "_desc-brain_mask." + EXT[output_type],
                        ),
                    )

            else:
                print("No output_directory was found...!\n")
                return

        if self.outputs:
            self.inheritance_dict[self.outputs["out_file"]] = self.in_file
            self.inheritance_dict[self.outputs["out_mask"]] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(SynthStripMriqc, self).run_process_mia()

        # necessary for speed gains (I think)
        torch.backends.cudnn.benchmark = True
        torch.backends.cudnn.deterministic = True

        # configure GPU device
        if self.gpu:
            os.environ["CUDA_VISIBLE_DEVICES"] = "0"
            device = torch.device("cuda")
            device_name = "GPU"
        else:
            os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
            device = torch.device("cpu")
            device_name = "CPU"

        # configure model
        print(f"SynthStripMriqc brick: Configuring model on the {device_name}")

        with torch.no_grad():
            model = self.StripModel()
            model.to(device)
            model.eval()

        # load model weights
        if self.model in [Undefined, "<undefined>", ""]:
            fconf = capsul.engine.configurations.get(
                "capsul.engine.module.freesurfer"
            )

            if fconf is not None:
                modelfile = os.path.join(
                    os.path.dirname(fconf["setup"]),
                    "models",
                    "synthstrip.1.pt",
                )

            else:
                raise RuntimeError(
                    "SynthStripMriqc: No configuration for "
                    "freeSurfer was found, so the path for "
                    "the model cannot be defined. "
                    "See File > MIA Preferences"
                )

        else:
            modelfile = self.model
            print("SynthStripMriqc brick: Using custom model weights")

        checkpoint = torch.load(modelfile, map_location=device)
        model.load_state_dict(checkpoint["model_state_dict"])

        # load input volume
        print(f"SynthStripMriqc brick: Input image read from: {self.in_file}")

        # normalize intensities
        image = nb.load(self.in_file)
        conformed = self.conform(image)
        in_data = conformed.get_fdata(dtype="float32")
        in_data -= in_data.min()
        in_data = np.clip(in_data / np.percentile(in_data, 99), 0, 1)
        in_data = in_data[np.newaxis, np.newaxis]

        # predict the surface distance transform
        input_tensor = torch.from_numpy(in_data).to(device)

        with torch.no_grad():
            sdt = model(input_tensor).cpu().numpy().squeeze()

        # unconform the sdt and extract mask
        sdt_target = self.resample_like(
            nb.Nifti1Image(sdt, conformed.affine, None),
            image,
            output_dtype="int16",
            cval=100,
        )
        sdt_data = np.asanyarray(sdt_target.dataobj).astype("int16")

        # find largest CC (just do this to be safe for now)
        components = scipy.ndimage.label(sdt_data.squeeze() < self.border_mm)[
            0
        ]
        bincount = np.bincount(components.flatten())[1:]
        mask = components == (np.argmax(bincount) + 1)
        mask = scipy.ndimage.morphology.binary_fill_holes(mask)

        # write the masked output
        if self.out_file not in [Undefined, "<undefined>"]:
            img_data = image.get_fdata()
            bg = np.min([0, img_data.min()])
            img_data[mask == 0] = bg
            nb.Nifti1Image(img_data, image.affine, image.header).to_filename(
                self.out_file
            )
            print(f"Masked image saved to: {self.out_file}")

        # write the brain mask
        if self.out_mask not in [Undefined, "<undefined>"]:
            hdr = image.header.copy()
            hdr.set_data_dtype("uint8")
            nb.Nifti1Image(mask, image.affine, hdr).to_filename(self.out_mask)
            print(
                f"SynthStripMriqc brick: "
                f"Binary brain mask saved to: {self.out_mask}"
            )

        print("If you use SynthStripMriqc in your analysis, please cite:")
        print("----------------------------------------------------")
        print("SynthStrip: Skull-Stripping for Any Brain Image.")
        print("A Hoopes, JS Mora, AV Dalca, B Fischl, M Hoffmann.")
