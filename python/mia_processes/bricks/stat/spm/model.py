# -*- coding: utf-8 -*-

"""The library for the SPM fMRI statistical analysis of the mia_processes
package.

The purpose of this module is to customise the main spm statistical analysis
bricks provided by nipype and to correct some things that do not work directly
in populse_mia.

:Contains:
    :Class:
        - EstimateContrast
        - EstimateModel
        - Level1Design

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

import numpy as np
import scipy.io

# nipype import
from nipype.interfaces.base import (
    File,
    InputMultiPath,
    OutputMultiPath,
    traits,
)
from nipype.interfaces.spm.base import ImageFileSPM
from populse_db.database import FIELD_TYPE_INTEGER, FIELD_TYPE_STRING
from populse_mia.data_manager.database_mia import TAG_ORIGIN_USER

# populse_db and populse_mia import
from populse_mia.data_manager.project import COLLECTION_CURRENT

# populse_mia imports
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# soma-base imports
from soma.qt_gui.qt_backend.Qt import QMessageBox
from traits.api import Undefined

# mia_processes import
from mia_processes.utils import get_dbFieldValue, set_dbFieldValue

# from .stats_pop_ups import _SessionQuery # currently in ec_dev package


class EstimateContrast(ProcessMIA):
    """
    * Estimate contrasts of interest.*

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(EstimateContrast, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["spm", "nipype"]

        # Inputs description
        spm_mat_file_desc = (
            "SPM.mat file (a pathlike object or string " "representing a file)"
        )
        session_type_desc = (
            "Selects the contrast type. One of tcon "
            "(T-contrast), fcon (F-contrast) or tconsess "
            "(T-contrast cond/sess based)"
        )
        contrast_name_desc = "Name of contrasts (a list of string)"
        condition_name_desc = (
            "Conditions information (a list of list of string)"
        )
        contrast_weight_desc = "Contrast weights (list of list of float)"
        session_desc = (
            "Session list (a list of list of float). "
            "Length should be equal to the number of sessions, "
            "with 1. for sessions to include and 0. elsewhere"
        )
        multi_reg_desc = "The regressor files (a list of file)"
        residual_image_desc = "Mean-squared image of the residuals"
        beta_images_desc = "Parameter estimates of the design matrix"
        # out_dir_name_desc = 'Out directory name'

        # TODO : understand use_derivs and group_contrast
        # use_derivs_desc = ("Use derivatives for estimation. "
        #                    "Mutually exclusive with inputs: group_contrast")
        # group_contrast_desc = ("Higher level contrast. "
        #                        "Mutually exclusive with inputs: use_derivs.")

        # Outputs description
        spmT_images_desc = "Stat images from a t-contrast"
        con_images_desc = "Contrast images from a t-contrast"

        # This bricks only works for T contrast
        # ess_images_desc = "Contrast images from an F-contrast"
        # spmF_images_desc = "Stat images from an F-contrast"
        out_spm_mat_file_desc = "Updated SPM mat file."

        # Inputs
        self.add_trait(
            "spm_mat_file",
            File(output=False, copyfile=True, desc=spm_mat_file_desc),
        )

        self.add_trait(
            "session_type",
            traits.Enum(
                "tcon",
                "fcon",
                "tconsess",
                output=False,
                optional=True,
                desc=session_type_desc,
            ),
        )

        self.add_trait(
            "contrast_name",
            traits.List(
                traits.String(),
                output=False,
                optional=True,
                desc=contrast_name_desc,
            ),
        )
        self.contrast_name = ["+"]

        self.add_trait(
            "condition_name",
            traits.List(
                traits.List(traits.String()),
                output=False,
                optional=True,
                desc=condition_name_desc,
            ),
        )
        self.condition_name = [["R1_1"]]

        self.add_trait(
            "contrast_weight",
            traits.List(
                traits.List(traits.Float()),
                output=False,
                optional=True,
                desc=contrast_weight_desc,
            ),
        )
        self.contrast_weight = [[1.0]]

        self.add_trait(
            "session",
            traits.Either(
                Undefined,
                traits.List(
                    traits.Either(Undefined, traits.List(traits.Float()))
                ),
                output=False,
                optional=True,
                desc=session_desc,
            ),
        )
        self.session = Undefined

        self.add_trait(
            "multi_reg",
            traits.Either(
                Undefined,
                traits.List(File()),
                output=False,
                optional=True,
                desc=multi_reg_desc,
            ),
        )

        self.add_trait(
            "beta_images",
            InputMultiPath(
                File(), output=False, copyfile=False, desc=beta_images_desc
            ),
        )
        self.add_trait(
            "residual_image",
            File(output=False, copyfile=False, desc=residual_image_desc),
        )
        # self.add_trait(
        #     "use_derivs",
        #     traits.Bool(output=False, optional=True, xor=["group_contrast"]),
        # )
        # self.add_trait(
        #     "group_contrast",
        #     traits.Bool(output=False, optional=True, xor=["use_derivs"]),
        # )

        # Outputs
        self.add_trait(
            "con_images",
            OutputMultiPath(
                File(), optional=True, output=True, desc=con_images_desc
            ),
        )

        self.add_trait(
            "spmT_images",
            OutputMultiPath(
                File(), optional=True, output=True, desc=spmT_images_desc
            ),
        )

        # self.add_trait("ess_images",
        #                OutputMultiPath(File(),
        #                                optional=True,
        #                                output=True,
        #                                desc=ess_images_desc)
        #                )

        # self.add_trait("spmF_images",
        #                OutputMultiPath(File(),
        #                                optional=True,
        #                                output=True,
        #                                desc=spmF_images_desc)
        #                )

        self.add_trait(
            "out_spm_mat_file",
            File(output=True, copyfile=False, desc=out_spm_mat_file_desc),
        )

        self.init_default_traits()
        self.init_process("nipype.interfaces.spm.EstimateContrast")

    def _get_contrasts(self, session_type):
        """
        Get list of contrasts with each contrast being a list of the form:
        [('contrast_name',
          'stat',
          [condition list],
          [weight list],
          [session list])]
        """

        contrasts = [tuple()]

        if session_type == "tcon":
            stat = "T"
        # elif session_type == "fcon":
        #     stat = "F"
        else:
            stat = None

        if stat is not None:
            if self.session is Undefined:
                contrasts = [
                    (cont_name, stat, condition, cont_weight)
                    for cont_name, condition, cont_weight in zip(
                        self.contrast_name,
                        self.condition_name,
                        self.contrast_weight,
                    )
                ]
            else:
                contrasts = [
                    (cont_name, stat, condition, cont_weight, sess)
                    if sess is not Undefined
                    else (cont_name, stat, condition, cont_weight)
                    for cont_name, condition, cont_weight, sess in zip(
                        self.contrast_name,
                        self.condition_name,
                        self.contrast_weight,
                        self.session,
                    )
                ]

        return contrasts

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
        super(EstimateContrast, self).list_outputs()

        # TODO : add F contrast

        # We cannot directly called Nipype's method to get all the output
        # because we want to avoid to read in the SPM.mat file
        # (potentially not created yet)

        if (self.spm_mat_file) and (
            self.spm_mat_file not in ["<undefined>", Undefined]
        ):
            # The management of self.process.output_directory could be
            # delegated to the
            # populse_mia.user_interface.pipeline_manager.process_mia
            # module. We can't do it at the moment because the
            # sync_process_output_traits() of the
            # capsul/process/nipype_process module raises an exception
            # in nipype if the mandatory parameter are not yet defined!

            spm_mat_dir, spm_mat_file = os.path.split(self.spm_mat_file)

            if self.output_directory:
                # Change output_directory for this process in order to
                # use a specific directory for each analysis

                if os.path.dirname(
                    spm_mat_dir
                ) == self.output_directory and "data" in os.path.basename(
                    spm_mat_dir
                ):
                    # if spm_mat already in a subfolder for a analysis
                    out_directory = spm_mat_dir
                else:
                    # if spm_mat file not in a subfolder(for e.g spm_mat
                    # file in download data)
                    sub_name = get_dbFieldValue(
                        self.project, self.spm_mat_file, "PatientName"
                    )
                    if sub_name is None:
                        print(
                            "Please, fill 'PatientName' tag "
                            "in the database for SPM file"
                        )
                        return self.make_initResult()

                    out_directory = os.path.join(
                        self.output_directory, sub_name + "data"
                    )

                    if not os.path.exists(out_directory):
                        os.mkdir(out_directory)
                self.output_directory = out_directory
            else:
                print("No output_directory was found...!\n")
                return self.make_initResult()

            self.outputs["out_spm_mat_file"] = os.path.join(
                self.output_directory, spm_mat_file
            )

            # Counting the number of spmT and con files to create
            # TODO: understand how to know exactly number of spmT/spmF
            nb_spmT = 0
            if self.multi_reg not in [Undefined, "<undefined>", [], None]:
                for reg_file in self.multi_reg:
                    if os.path.splitext(reg_file)[1] != ".txt":
                        nb_spmT += 1
            else:
                nb_spmT = len(self.contrast_name)

            spmT_files = []
            con_files = []

            for i in range(nb_spmT):
                i += 1
                spmT_files.append(
                    os.path.join(
                        self.output_directory, "spmT_{:04d}.nii".format(i)
                    )
                )
                con_files.append(
                    os.path.join(
                        self.output_directory, "con_{:04d}.nii".format(i)
                    )
                )
                # spmT_files.append(os.path.join(path,
                #                                'spmT_{:04d}.nii'.format(i)))
                # con_files.append(os.path.join(path,
                #                               'con_{:04d}.nii'.format(i)))

            if spmT_files:
                self.outputs["spmT_images"] = spmT_files
                self.outputs["con_images"] = con_files

        if self.outputs:
            self.inheritance_dict[
                self.outputs["out_spm_mat_file"]
            ] = self.spm_mat_file

            if spmT_files:
                # FIXME: Quick and dirty. This brick was written quickly and
                #        will need to be reworked to cover all cases. At that
                #        time, inheritance will also need to be reviewed.
                #        Currently we take the first element of the lists, but
                #        in the general case there can be several elements in
                #        the lists
                self.inheritance_dict[
                    self.outputs["spmT_images"][0]
                ] = self.spm_mat_file
                # FIXME: In the latest version of mia, indexing of the
                #        database with particular tags defined in the
                #        processes is done only at the end of the
                #        initialisation of the whole pipeline. So we
                #        cannot use the value of these tags in other
                #        processes of the pipeline at the time of
                #        initialisation (see populse_mia #290). Unti
                #        better we use a quick and dirty hack with the
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
                    self.project, self.spm_mat_file, "PatientName"
                )

                if tag_to_add["value"] is not None:
                    set_dbFieldValue(self.project, spmT_files[0], tag_to_add)

                else:
                    print(
                        "\nEstimateContrast brick:\nThe 'PatientName' "
                        "tag could not be added to the database for "
                        "the '{}' parameter. This can lead to a "
                        "subsequent issue during "
                        "initialization!!\n".format(spmT_files[0])
                    )

                age = get_dbFieldValue(self.project, self.spm_mat_file, "Age")

                if age is not None:
                    tag_to_add = dict()
                    tag_to_add["name"] = "Age"
                    tag_to_add["field_type"] = "int"
                    tag_to_add["description"] = ""
                    tag_to_add["visibility"] = True
                    tag_to_add["origin"] = "user"
                    tag_to_add["unit"] = None
                    tag_to_add["default_value"] = None
                    tag_to_add["value"] = age
                    set_dbFieldValue(self.project, spmT_files[0], tag_to_add)

                pathology = get_dbFieldValue(
                    self.project, self.spm_mat_file, "Pathology"
                )

                if pathology is not None:
                    tag_to_add = dict()
                    tag_to_add["name"] = "Pathology"
                    tag_to_add["field_type"] = "string"
                    tag_to_add["description"] = ""
                    tag_to_add["visibility"] = True
                    tag_to_add["origin"] = "user"
                    tag_to_add["unit"] = None
                    tag_to_add["default_value"] = None
                    tag_to_add["value"] = pathology
                    set_dbFieldValue(self.project, spmT_files[0], tag_to_add)

                self.inheritance_dict[
                    self.outputs["con_images"][0]
                ] = self.spm_mat_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""

        super(EstimateContrast, self).run_process_mia()

        self.process.spm_mat_file = self.spm_mat_file
        self.process.contrasts = self._get_contrasts(self.session_type)
        self.process.beta_images = self.beta_images
        self.process.residual_image = self.residual_image

        # if self.use_derivs is not None:
        #     self.process.use_derivs = self.use_derivs
        # else:
        #     if self.group_contrast is not None:
        #         self.process.group_contrast = self.group_contrast

        #     else:
        #         self.process.use_derivs = False

        return self.process.run(configuration_dict={})


class EstimateModel(ProcessMIA):
    """
    * Estimation of model parameters using classical
    (ReML - Restricted Maximum Likelihood)*

    Please, see the complete documentation for the EstimateModel brick
    in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/html/documentation/bricks/stats/spm/EstimateModel.html>`

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(EstimateModel, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["spm", "nipype"]

        # Inputs description
        spm_mat_file_desc = (
            "The SPM.mat file that contains the design "
            "specification (a file object)"
        )
        # EstimateModel brick only works for classical method
        # estimation_method_desc = (
        #     "A dictionary of either Classical: 1, "
        #     "Bayesian: 1, or Bayesian2: 1"
        # )
        estimation_method_desc = "A dictionary : Classical: 1"
        write_residuals_desc = "Write individual residual images (a boolean)"
        # flags_desc = (
        #     "Additional arguments (a dictionary with keys which are "
        #     "any value and with values which are any value)"
        # )
        version_desc = "Version of spm (a string)"
        tot_reg_num_desc = (
            "The total number of estimated regression"
            "coefficients (an integer)"
        )

        # Outputs description
        out_spm_mat_file_desc = (
            "The file containing specification of the "
            "design and estimated model parameters (a "
            "pathlike object or string representing an "
            "existing file)"
        )
        beta_images_desc = (
            "Images of estimated regression coefficients (a "
            "list of items which are a pathlike object or "
            "string representing an existing file)"
        )
        mask_image_desc = (
            "The image indicating which voxels were included in "
            "the analysis (a pathlike object or string "
            "representing an existing file)"
        )
        residual_image_desc = (
            "The image of the variance of the error (a "
            "pathlike object or string representing an "
            "existing file)"
        )
        residual_images_desc = (
            "The individual error images (a list of items "
            "which are a pathlike object or string "
            "representing an existing file)"
        )
        RPVimage_desc = (
            "The image of the estimated resolution elements per "
            "voxel (a pathlike object or string representing an "
            "existing file)."
        )
        factor_info_desc = (
            "A list of items which are a dictionary for each "
            "factor with keys which are 'name' or 'levels' "
            "and with values which are a string (name of the "
            "factor) or an integer (number of levels for the "
            "factor)"
        )
        bases_desc = (
            "To define basic functions for modeling hemodynamic "
            "response (a 'none' string or a dictionary with keys "
            "which are 'hrf' or 'fourier' or 'fourier_han' or "
            "'gamma' or 'fir' and with values which are a dictionary "
            "with keys which are 'derivs' or 'length' or 'order' and "
            "with values which are a list or a float or an integer)"
        )
        spmT_images_desc = "Stat images from a t-contrast"
        con_images_desc = "Contrast images from a t-contrast"
        ess_images_desc = "Contrast images from an F-contrast"
        spmF_images_desc = "Stat images from an F-contrast"

        # EstimateModel brick only works for classical method
        # Following outputs are generated by
        # Bayesian 1st-level or Bayesian 2nd level methods
        #
        # labels_desc = "Label file"
        # SDerror_desc = "Images of the standard deviation of the error"
        # ARcoef_desc = (
        #     "Images of the autoregressive (AR) coefficient (a list "
        #     "of items which are a pathlike object or string "
        #     "representing an existing file)."
        # )
        # Cbetas_desc = "Images of the parameter posteriors"
        # SDbetas_desc = ("Images of the standard deviation of "
        #                 "parameter posteriors")

        # Inputs traits
        self.add_trait(
            "spm_mat_file",
            File(output=False, optional=False, desc=spm_mat_file_desc),
        )

        # EstimateModel brick only works for classical method
        self.add_trait(
            "estimation_method",
            traits.Dict(
                # traits.Enum("Classical", "Bayesian", "Bayesian2"),
                traits.Enum("Classical"),
                traits.Enum(1),
                usedefault=True,
                output=False,
                optional=True,
                desc=estimation_method_desc,
            ),
        )
        self.estimation_method = {"Classical": 1}

        self.add_trait(
            "write_residuals",
            traits.Bool(
                usedefault=True,
                output=False,
                optional=True,
                desc=write_residuals_desc,
            ),
        )

        # self.add_trait(
        #     "flags", traits.Dict(output=False,
        #     optional=True, desc=flags_desc)
        # )

        self.add_trait(
            "version",
            traits.String(
                "spm12",
                usedefault=True,
                output=False,
                optional=True,
                desc=version_desc,
            ),
        )

        self.add_trait(
            "tot_reg_num",
            traits.Either(
                (traits.Int(), Undefined),
                output=False,
                optional=True,
                desc=tot_reg_num_desc,
            ),
        )
        self.tot_reg_num = Undefined

        # we need factor info and bases
        # to know if contrasts will be estimated by SPM
        self.add_trait(
            "factor_info",
            traits.List(
                traits.Dict(
                    traits.Enum("name", "levels"),
                    traits.Either(traits.Str, traits.Int),
                ),
                usedefault=True,
                output=False,
                optional=True,
                desc=factor_info_desc,
            ),
        )

        self.add_trait(
            "bases",
            traits.Union(
                traits.Dict(
                    traits.Enum(
                        "hrf", "fourier", "fourier_han", "gamma", "fir"
                    ),
                    traits.Dict(
                        traits.Enum("derivs", "length", "order"),
                        traits.Union(
                            traits.Enum([0, 0], [1, 0], [1, 1]),
                            traits.Int,
                            traits.Float,
                        ),
                    ),
                ),
                traits.Enum(["none"]),
                usedefault=True,
                output=False,
                optional=True,
                desc=bases_desc,
            ),
        )
        self.bases = {"hrf": {"derivs": [0, 0]}}

        # Output traits
        self.add_trait(
            "out_spm_mat_file",
            File(output=True, optional=False, desc=out_spm_mat_file_desc),
        )

        self.add_trait(
            "beta_images",
            OutputMultiPath(
                ImageFileSPM(),
                output=True,
                optional=True,
                desc=beta_images_desc,
            ),
        )

        self.add_trait(
            "mask_image",
            ImageFileSPM(output=True, optional=True, desc=mask_image_desc),
        )

        self.add_trait(
            "residual_image",
            ImageFileSPM(output=True, optional=True, desc=residual_image_desc),
        )

        self.add_trait(
            "residual_images",
            OutputMultiPath(
                ImageFileSPM(),
                output=True,
                optional=True,
                desc=residual_images_desc,
            ),
        )

        self.add_trait(
            "RPVimage",
            ImageFileSPM(output=True, optional=True, desc=RPVimage_desc),
        )

        # If factorial design specified in Level1design
        # SPM also write T and F contrasts
        self.add_trait(
            "con_images",
            OutputMultiPath(
                File(), optional=True, output=True, desc=con_images_desc
            ),
        )

        self.add_trait(
            "spmT_images",
            OutputMultiPath(
                File(), optional=True, output=True, desc=spmT_images_desc
            ),
        )

        self.add_trait(
            "ess_images",
            OutputMultiPath(
                File(), optional=True, output=True, desc=ess_images_desc
            ),
        )

        self.add_trait(
            "spmF_images",
            OutputMultiPath(
                File(), optional=True, output=True, desc=spmF_images_desc
            ),
        )

        # EstimateModel brick only works for classical method
        # Following outputs are generated by
        # Bayesian 1st-level or Bayesian 2nd level methods

        # self.add_trait(
        #     "labels",
        #     ImageFileSPM(output=True, optional=True, desc=labels_desc),
        # )

        # self.add_trait(
        #     "SDerror",
        #     OutputMultiPath(
        #         ImageFileSPM(), output=True, optional=True, desc=SDerror_desc
        #     ),
        # )

        # self.add_trait(
        #     "ARcoef",
        #     OutputMultiPath(
        #         ImageFileSPM(), output=True, optional=True, desc=ARcoef_desc
        #     ),
        # )

        # self.add_trait(
        #     "Cbetas",
        #     OutputMultiPath(
        #         ImageFileSPM(), output=True, optional=True, desc=Cbetas_desc
        #     ),
        # )

        # self.add_trait(
        #     "SDbetas",
        #     OutputMultiPath(
        #         ImageFileSPM(), output=True, optional=True, desc=SDbetas_desc
        #     ),
        # )

        self.init_default_traits()
        self.init_process("nipype.interfaces.spm.EstimateModel")

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dict),
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key 'notInDb' of the
        self.outputs dictionary. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(EstimateModel, self).list_outputs()

        # We cannot directly called Nipype's method to get all the output
        # because we want to avoid to read in the SPM.mat file
        # (potentially not created yet)
        if (
            (self.spm_mat_file)
            and (self.spm_mat_file not in ["<undefined>", Undefined])
            and (self.estimation_method)
            and (self.estimation_method not in ["<undefined>", Undefined])
        ):
            im_form = "nii" if "12" in self.version else "img"
            spm_mat_dir, spm_mat_file = os.path.split(self.spm_mat_file)
            # The management of self.process.output_directory could be
            # delegated to the
            # populse_mia.user_interface.pipeline_manager.process_mia
            # module. We can't do it at the moment because the
            # sync_process_output_traits() of the
            # capsul/process/nipype_process module raises an exception
            # in nipype if the mandatory parameter are not yet defined!

            if self.output_directory:
                # Change output_directory for this process in order to
                # use a specific directory for each analysis

                if os.path.dirname(
                    spm_mat_dir
                ) == self.output_directory and "data" in os.path.basename(
                    spm_mat_dir
                ):
                    # if spm_mat already in a subfolder for a analysis
                    out_directory = spm_mat_dir
                else:
                    # if spm_mat file not in a subfolder(for e.g spm_mat
                    # file in download data)
                    sub_name = get_dbFieldValue(
                        self.project, self.spm_mat_file, "PatientName"
                    )
                    if sub_name is None:
                        print(
                            "Please, fill 'PatientName' tag "
                            "in the database for SPM file"
                        )
                        return self.make_initResult()

                    out_directory = os.path.join(
                        self.output_directory, sub_name + "data"
                    )

                    if not os.path.exists(out_directory):
                        os.mkdir(out_directory)
                self.output_directory = out_directory
            else:
                print("No output_directory was found...!\n")
                return self.make_initResult()

            self.outputs["out_spm_mat_file"] = os.path.join(
                self.output_directory, spm_mat_file
            )
            # Detecting the number of beta files to create
            betas = []

            if self.tot_reg_num in ["<undefined>", Undefined]:
                tot_reg_numb = get_dbFieldValue(
                    self.project, self.spm_mat_file, "Regress num"
                )

                if tot_reg_numb is not None:
                    self.tot_reg_num = tot_reg_numb

                else:
                    print(
                        '\nEstimateModel brick:\n The "tot_reg_num" parameter'
                        " could not be determined automatically because the "
                        '"Regress num" tag for the {} file is not filled '
                        "in the database. Please set this value and launch "
                        "the calculation "
                        "again!\n".format(self.spm_mat_file)
                    )
                    self.outputs = {}
                    return self.make_initResult()

            # Bayesian and Bayesian2 are not yet fully implemented
            # if ("Bayesian" in self.estimation_method.keys()) or (
            #     "Bayesian2" in self.estimation_method.keys()
            # ):
            #     self.outputs["labels"] = os.path.join(
            #         self.output_directory, "labels.{}".format(im_form)
            #     )
            #     """
            #     outputs['SDerror'] = glob(
            #                   os.path.join(path, 'Sess*_SDerror*'))
            #     outputs['ARcoef'] = glob(os.path.join(path, 'Sess*_AR_*'))
            #     """

            if "Classical" in self.estimation_method.keys():
                self.outputs["residual_image"] = os.path.join(
                    self.output_directory, "ResMS.{}".format(im_form)
                )
                self.outputs["RPVimage"] = os.path.join(
                    self.output_directory, "RPV.{}".format(im_form)
                )
                self.outputs["mask_image"] = os.path.join(
                    self.output_directory, "mask.{}".format(im_form)
                )

                if self.write_residuals:
                    nb_dyn = get_dbFieldValue(
                        self.project, self.spm_mat_file, "Dynamic Number"
                    )

                    if nb_dyn:
                        nb_dyn = sum(nb_dyn)
                        ress = []

                        for i in range(nb_dyn):
                            i += 1
                            ress.append("Res_{:04d}.{}".format(i, im_form))

                        if ress:
                            self.outputs["residual_images"] = [
                                os.path.join(self.output_directory, res)
                                for res in ress
                            ]

                    else:
                        print(
                            "\nEstimateModel brick:\nThe number of dynamics "
                            "could not be determined automatically because "
                            'the "Dynamic Number" tag is not filled in the'
                            "database for the {} file.\nAs a result, it is "
                            "not possible to safely create the "
                            '"residual_images" output '
                            "parameter!\n".format(self.spm_mat_file)
                        )

                for i in range(
                    int(
                        0
                        if self.tot_reg_num in ["<undefined>", Undefined, None]
                        else self.tot_reg_num
                    )
                ):
                    i += 1
                    betas.append("beta_{:04d}.{}".format(i, im_form))

                if betas:
                    self.outputs["beta_images"] = [
                        os.path.join(self.output_directory, beta)
                        for beta in betas
                    ]

                else:
                    print(
                        "\nEstimateModel brick:\nNo beta image were found to"
                        " be added to the database ..."
                    )

            if (self.tot_reg_num in ["<undefined>", Undefined, None]) or (
                self.write_residuals and not nb_dyn
            ):
                self.outputs = {}

            # If factorial design used in levelOne  Design,
            # detecting the number of contrasts files to create
            if self.factor_info:
                number_factors = len(self.factor_info)
                number_interactions = int(
                    (number_factors * (number_factors - 1)) / 2
                )
                multiplier = 1
                if "hrf" in self.bases.keys():
                    if "derivs" in self.bases["hrf"].keys():
                        multiplier += sum(self.bases["hrf"]["derivs"])
                # number of F contarst =
                # average effect of condition + main effect for each factor
                # + interaction between each factor
                number_f_contrast = 1 + number_factors + number_interactions
                # number_of_T_contrast =
                # positive effect of condition * derivs multiplier +
                # positive effect of each factor * derivs multiplier +
                # positive interaction between each factor
                number_t_contrast = multiplier * (
                    1 + number_factors + number_interactions
                )

                cons = []
                spmTs = []
                esss = []
                spmFs = []
                for i in range(number_f_contrast):
                    i += 1
                    esss.append("ess_{:04d}.{}".format(i, im_form))
                    spmTs.append("spmF_{:04d}.{}".format(i, im_form))

                    if esss:
                        self.outputs["ess_images"] = [
                            os.path.join(self.output_directory, ess)
                            for ess in esss
                        ]
                    if spmFs:
                        self.outputs["spmF_images"] = [
                            os.path.join(self.output_directory, spmF)
                            for spmF in spmFs
                        ]
                for i in range(number_t_contrast):
                    i += 1
                    i += number_f_contrast
                    cons.append("con_{:04d}.{}".format(i, im_form))
                    spmTs.append("spmT_{:04d}.{}".format(i, im_form))

                    if cons:
                        self.outputs["con_images"] = [
                            os.path.join(self.output_directory, con)
                            for con in cons
                        ]
                    if spmTs:
                        self.outputs["spmT_images"] = [
                            os.path.join(self.output_directory, spmT)
                            for spmT in spmTs
                        ]

        if self.outputs:
            for key, value in self.outputs.items():
                if key in [
                    "out_spm_mat_file",
                    "mask_image",
                    "residual_image",
                    "RPVimage",
                ]:
                    self.inheritance_dict[value] = self.spm_mat_file

                elif key in "residual_images":
                    for fullname in value:
                        self.inheritance_dict[fullname] = self.spm_mat_file

                elif key in [
                    "beta_images",
                    "con_images",
                    "spmT_images",
                    "ess_images",
                    "spmF_images",
                ]:
                    patient_name = get_dbFieldValue(
                        self.project, self.spm_mat_file, "PatientName"
                    )

                    for fullname in value:
                        self.inheritance_dict[fullname] = self.spm_mat_file
                        # FIXME: In the latest version of mia, indexing of the
                        #        database with particular tags defined in the
                        #        processes is done only at the end of the
                        #        initialisation of the whole pipeline. So we
                        #        cannot use the value of these tags in other
                        #        processes of the pipeline at the time of
                        #        initialisation (see populse_mia #290). Unti
                        #        better we use a quick and dirty hack with the
                        #        set_dbFieldValue() function !

                        if patient_name is not None:
                            tag_to_add = dict()
                            tag_to_add["name"] = "PatientName"
                            tag_to_add["field_type"] = "string"
                            tag_to_add["description"] = ""
                            tag_to_add["visibility"] = True
                            tag_to_add["origin"] = "user"
                            tag_to_add["unit"] = None
                            tag_to_add["default_value"] = None
                            tag_to_add["value"] = patient_name
                            set_dbFieldValue(
                                self.project, fullname, tag_to_add
                            )

                        else:
                            print(
                                "\nEstimateModel brick:\nThe 'PatientName' "
                                "tag could not be added to the database for "
                                "the '{}' parameter. This can lead to a "
                                "subsequent issue during "
                                "initialization!!\n".format(fullname)
                            )

                        age = get_dbFieldValue(
                            self.project, self.spm_mat_file, "Age"
                        )

                        if age is not None:
                            tag_to_add = dict()
                            tag_to_add["name"] = "Age"
                            tag_to_add["field_type"] = "int"
                            tag_to_add["description"] = ""
                            tag_to_add["visibility"] = True
                            tag_to_add["origin"] = "user"
                            tag_to_add["unit"] = None
                            tag_to_add["default_value"] = None
                            tag_to_add["value"] = age
                            set_dbFieldValue(
                                self.project, fullname, tag_to_add
                            )

                        pathology = get_dbFieldValue(
                            self.project, self.spm_mat_file, "Pathology"
                        )

                        if pathology is not None:
                            tag_to_add = dict()
                            tag_to_add["name"] = "Pathology"
                            tag_to_add["field_type"] = "string"
                            tag_to_add["description"] = ""
                            tag_to_add["visibility"] = True
                            tag_to_add["origin"] = "user"
                            tag_to_add["unit"] = None
                            tag_to_add["default_value"] = None
                            tag_to_add["value"] = pathology
                            set_dbFieldValue(
                                self.project, fullname, tag_to_add
                            )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(EstimateModel, self).run_process_mia()
        self.process.spm_mat_file = os.path.abspath(self.spm_mat_file)
        self.process.estimation_method = self.estimation_method
        self.process.write_residuals = self.write_residuals

        # self.process.flags = self.flags
        # Removing the image files to avoid a bug
        # I do not observe the bug now, i comment the following lines:
        # TODO: In fact, self.outputs is == {} at this point(see issue #272).
        #       If necessary, we can use the dict4runtine object?
        #
        # for key, value in self.outputs.items():
        #   if key not in ["out_spm_mat_file"]:
        #        if value not in ["<undefined>", Undefined]:

        #            if type(value) in [list, traits.TraitListObject,
        #                               traits.List]:

        #                for element in value:
        #                    if os.path.isfile(element):
        #                        os.remove(element)
        #           else:
        #                if os.path.isfile(value):
        #                    os.remove(value)
        # self.process.run()
        return self.process.run(configuration_dict={})


class Level1Design(ProcessMIA):
    """
    * Specification of the GLM design matrix, fMRI data files and filtering.

    Please, see the complete documentation for the Level1Design brick
    in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/html/documentation/bricks/stats/spm/LevelDesign.html>

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Level1Design, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ["spm", "nipype"]

        # Inputs description
        timing_units_desc = "One of 'scans' or 'secs'"
        interscan_interval_desc = "TR in secs (a float)"
        microtime_resolution_desc = (
            "The number of time-bins per scan used "
            "to build regressors (an integer)"
        )
        microtime_onset_desc = (
            "The reference time-bin at which the "
            "regressors are resampled to coincide with "
            "data acquisition (an integer)"
        )
        sess_scans_desc = (
            "The fMRI scan for each session (a list of items "
            "which are a pathlike object or string representing "
            "an existing file)"
        )
        sess_cond_names_desc = (
            "The name of each condition (list of items "
            "which are a list of items which are a "
            "string)"
        )
        sess_cond_onsets_desc = (
            "The onset times (in seconds or in scans) of "
            "the epochs or events within each condition "
            "(a list of items which are a list of items "
            "which are a list of items which are a "
            "float)"
        )
        sess_cond_durations_desc = (
            "The duration times (in seconds or in "
            "scans) of the epochs within each "
            "condition (a list of items which are a "
            "list of items which are a list of items "
            "which are a float)"
        )
        sess_cond_tmod_desc = (
            "Allows for the characterisation of linear or "
            "nonlinear time effects (a list of items which "
            "are a list of items which are 0 or 1 or 2 or 3 "
            "or 4 or 5 or 6)"
        )
        sess_cond_pmod_names_desc = (
            "The name of the parametric modulation (a "
            "list of items which are a list of items "
            "which are a list of items which are a "
            "string)"
        )
        sess_cond_pmod_values_desc = (
            "The values used for the parametric "
            "modulation, one for each occurence of "
            "the event (a list of items which are a "
            "list of items which are a list of items "
            "which are a list of items which are a "
            "float)"
        )
        sess_cond_pmod_polys_desc = (
            "The polynomial expansion used for the "
            "parametric modulation (a list of items "
            "which are a list of items which are a "
            "list of items which are 1 or 2 or 3 or 4 "
            "or 5 or 6)"
        )
        sess_cond_orth_desc = (
            "Orthogonalise regressors within trial types "
            "(a list of items which are a list of items "
            "which are 0 or 1)"
        )
        sess_multi_desc = (
            "A *.mat file containing details of the multiple "
            "experimental conditions for each session (a list "
            "of items which are a filename)"
        )
        sess_regress_desc = (
            "Additional columns included in the design "
            "matrix, which may model effects that would not "
            "be convolved with the haemodynamic response (a "
            "list of items which are a list of items which "
            "are a dictionary with keys which are 'name' or "
            "'val' and with values which are a string or a "
            "list of float)"
        )
        sess_multi_reg_desc = (
            "The .mat/.txt file(s) containing details of "
            "multiple regressors (a list of items which are "
            "a filename)"
        )
        sess_hpf_desc = "High-pass filter (a list of items which are a float)"
        factor_info_desc = (
            "A list of items which are a dictionary for each "
            "factor with keys which are 'name' or 'levels' "
            "and with values which are a string (name of the "
            "factor) or an integer (number of levels for the "
            "factor)"
        )
        bases_desc = (
            "To define basic functions for modeling hemodynamic "
            "response (a 'none' string or a dictionary with keys "
            "which are 'hrf' or 'fourier' or 'fourier_han' or "
            "'gamma' or 'fir' and with values which are a dictionary "
            "with keys which are 'derivs' or 'length' or 'order' and "
            "with values which are a list or a float or an integer)"
        )
        volterra_expansion_order_desc = (
            "One of 1 or 2 (1: do not model "
            "interactions, 2: model interactions)"
        )
        global_intensity_normalization_desc = (
            "Global intensity normalisation "
            "with scaling or not (one of "
            "'none' or 'scaling')"
        )
        mask_threshold_desc = (
            "Masking threshold, defined as proportion of " "globals (a float)"
        )
        mask_image_desc = (
            "Image for explicitly masking the analysis (a "
            "pathlike object or string representing a file)"
        )
        model_serial_correlations_desc = (
            "one of AR(1), or FAST or none "
            "(AR(1): autoregressive model, "
            "FAST: available in SPM12, "
            "none: serial correlation is "
            "ignored)"
        )

        # Outputs description
        spm_mat_file_desc = (
            "SPM.mat file (a pathlike object or string " "representing a file"
        )

        # Inputs traits
        self.add_trait(
            "timing_units",
            traits.Enum(
                "scans",
                "secs",
                usedefault=True,
                output=False,
                optional=True,
                desc=timing_units_desc,
            ),
        )

        # Plan to retrieve this parameter automatically from the database?
        self.add_trait(
            "interscan_interval",
            traits.Either(
                traits.Float(),
                Undefined,
                usedefault=True,
                output=False,
                optional=True,
                desc=interscan_interval_desc,
            ),
        )
        self.interscan_interval = Undefined

        self.add_trait(
            "microtime_resolution",
            traits.Int(
                16,
                usedefault=True,
                output=False,
                optional=True,
                desc=microtime_resolution_desc,
            ),
        )

        # In study without slice-timing correction, as cevastoc32,
        # it should be 8?
        self.add_trait(
            "microtime_onset",
            traits.Int(
                8,
                usedefault=True,
                output=False,
                optional=True,
                desc=microtime_onset_desc,
            ),
        )
        # TODO: mictotime onset (fmri_spec.timing.fmri_t0) = 1 in Amigo

        self.add_trait(
            "sess_scans",
            InputMultiPath(
                File(exists=True),
                usedefault=True,
                output=False,
                desc=sess_scans_desc,
            ),
        )

        self.add_trait(
            "sess_cond_names",
            traits.List(
                traits.Either(traits.List(traits.String()), None),
                value=[[]],
                output=False,
                optional=True,
                desc=sess_cond_names_desc,
            ),
        )

        self.add_trait(
            "sess_cond_onsets",
            traits.List(
                traits.Either(traits.List(traits.List(traits.Float())), None),
                value=[[[]]],
                output=False,
                optional=True,
                desc=sess_cond_onsets_desc,
            ),
        )

        self.add_trait(
            "sess_cond_durations",
            traits.List(
                traits.Either(traits.List(traits.List(traits.Float())), None),
                value=[[[]]],
                output=False,
                optional=True,
                desc=sess_cond_durations_desc,
            ),
        )

        self.add_trait(
            "sess_cond_tmod",
            traits.List(
                traits.Either(
                    traits.List(traits.Enum(0, 1, 2, 3, 4, 5, 6)), None
                ),
                value=[[0]],
                output=False,
                optional=True,
                desc=sess_cond_tmod_desc,
            ),
        )

        self.add_trait(
            "sess_cond_pmod_names",
            traits.List(
                traits.Either(
                    traits.List(
                        traits.Either(traits.List(traits.String()), None)
                    ),
                    None,
                ),
                value=[[[]]],
                output=False,
                optional=True,
                desc=sess_cond_pmod_names_desc,
            ),
        )

        self.add_trait(
            "sess_cond_pmod_values",
            traits.List(
                traits.Either(
                    traits.List(
                        traits.Either(
                            traits.List(traits.List(traits.Float())), None
                        )
                    ),
                    None,
                ),
                value=[[[[]]]],
                output=False,
                optional=True,
                desc=sess_cond_pmod_values_desc,
            ),
        )

        self.add_trait(
            "sess_cond_pmod_polys",
            traits.List(
                traits.Either(
                    traits.List(
                        traits.Either(
                            traits.List(traits.Enum(1, 2, 3, 4, 5, 6)), None
                        )
                    ),
                    None,
                ),
                value=[[[]]],
                output=False,
                optional=True,
                desc=sess_cond_pmod_polys_desc,
            ),
        )

        self.add_trait(
            "sess_cond_orth",
            traits.List(
                traits.Either(traits.List(traits.Enum(0, 1)), None),
                value=[[]],
                output=False,
                optional=True,
                desc=sess_cond_orth_desc,
            ),
        )

        self.add_trait(
            "sess_multi",
            traits.List(
                traits.File(),
                value=[],
                output=False,
                optional=True,
                desc=sess_multi_desc,
            ),
        )

        self.add_trait(
            "sess_regress",
            traits.List(
                traits.List(
                    traits.Dict(
                        traits.Enum("name", "val"),
                        traits.Union(traits.Str, traits.List(traits.Float())),
                    )
                ),
                value=[[]],
                output=False,
                optional=True,
                desc=sess_regress_desc,
            ),
        )

        self.add_trait(
            "sess_multi_reg",
            InputMultiPath(
                traits.Either(traits.List(traits.File()), None),
                value=[[]],
                output=False,
                optional=True,
                desc=sess_multi_reg_desc,
            ),
        )

        self.add_trait(
            "sess_hpf",
            traits.List(
                traits.Float(),
                value=[427.2],
                usedefault=True,
                output=False,
                optional=True,
                desc=sess_hpf_desc,
            ),
        )
        # TODO: 427.2 corresponds to the value used in Amigo
        #       Duration * TR * 3.56 = 427.2 if Duration == 40 and
        #       TR == 3s; why 3.56 ?
        #       I was expecting rather to:
        #       (time between first block start
        #        - second block start) * TR *2 = 80 *3 *2 = 480
        #       Can we code an automatic recovery procedure for the
        #       hpf parameter ?
        #       (in this case, the user would have to declare additional
        #       tags in the database, like the block duration !)

        self.add_trait(
            "factor_info",
            traits.List(
                traits.Dict(
                    traits.Enum("name", "levels"),
                    traits.Either(traits.Str, traits.Int),
                ),
                usedefault=True,
                output=False,
                optional=True,
                desc=factor_info_desc,
            ),
        )

        # traits.Union available from traits 6.1
        self.add_trait(
            "bases",
            traits.Union(
                traits.Dict(
                    traits.Enum(
                        "hrf", "fourier", "fourier_han", "gamma", "fir"
                    ),
                    traits.Dict(
                        traits.Enum("derivs", "length", "order"),
                        traits.Union(
                            traits.Enum([0, 0], [1, 0], [1, 1]),
                            traits.Int,
                            traits.Float,
                        ),
                    ),
                ),
                traits.Enum(["none"]),
                usedefault=True,
                output=False,
                optional=True,
                desc=bases_desc,
            ),
        )
        self.bases = {"hrf": {"derivs": [0, 0]}}

        self.add_trait(
            "volterra_expansion_order",
            traits.Enum(
                1,
                2,
                usedefault=True,
                output=False,
                optional=True,
                desc=volterra_expansion_order_desc,
            ),
        )

        self.add_trait(
            "global_intensity_normalization",
            traits.Enum(
                "none",
                "scaling",
                usedefault=True,
                output=False,
                optional=True,
                desc=global_intensity_normalization_desc,
            ),
        )

        self.add_trait(
            "mask_threshold",
            traits.Float(
                0.8,
                usedefault=True,
                output=False,
                optional=True,
                desc=mask_threshold_desc,
            ),
        )

        self.add_trait(
            "mask_image",
            ImageFileSPM(output=False, optional=True, desc=mask_image_desc),
        )

        self.add_trait(
            "model_serial_correlations",
            traits.Enum(
                "AR(1)",
                "FAST",
                "none",
                usedefault=True,
                output=False,
                optional=True,
                desc=model_serial_correlations_desc,
            ),
        )

        # Output traits
        self.add_trait(
            "spm_mat_file", File(output=True, desc=spm_mat_file_desc)
        )

        # Special parameter used as a messenger for the run_process_mia method
        self.add_trait(
            "dict4runtime",
            traits.Dict(output=False, optional=True, userlevel=1),
        )

        self.init_default_traits()
        self.init_process("nipype.interfaces.spm.Level1Design")

    def _get_conditions(self, idx_session, idx_cond):
        """Generate the condition dictionary.

        One condition dictionary is generated for each condition.
        The condition dictionary has the following keys:
        - name (a string)
        - onset (a list of float)
        - duration (a list of float)
        - tmod (an integer between 0 and 6)
        - pmod a dictionary
        - orth (0 or 1)

        :param idx_session: session index.
        :param idx_cond: condition index
        :returns: conditions dictionary
        """
        cond = dict()
        cond["name"] = self.sess_cond_names[idx_session][idx_cond]
        onset = []

        for i, value in enumerate(
            self.sess_cond_onsets[idx_session][idx_cond]
        ):
            onset.insert(i, np.float64(value))

        cond["onset"] = onset

        if (
            len(self.sess_cond_durations[idx_session][idx_cond])
            == len(self.sess_cond_onsets[idx_session][idx_cond])
        ) or (len(self.sess_cond_durations[idx_session][idx_cond]) == 1):
            duration = []

            for i, value in enumerate(
                self.sess_cond_durations[idx_session][idx_cond]
            ):
                duration.insert(i, np.float64(value))
            cond["duration"] = duration

        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle("Level1Design Error!")
            msg.setText(
                "Warning: The number of values in the "
                "sess_cond_durations parameter does not correspond to "
                "the number of values in the sess_cond_onsets "
                "parameter, for the session '{0}' and the "
                "condition '{1}'! Pease, check your "
                "settings ...".format(idx_session, cond["name"])
            )
            msg.setStandardButtons(QMessageBox.Close)
            msg.buttonClicked.connect(msg.close)
            msg.exec()
            return False

        try:
            cond["tmod"] = self.sess_cond_tmod[idx_session][idx_cond]

        except IndexError:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle("Level1Design Error!")
            msg.setText(
                "Warning: It seems that not all sess_cond_tmod "
                "parameter values are defined for all conditions. "
                "Each condition must have a value for the "
                "sess_cond_tmod parameter. Please, check the "
                "sess_cond_tmod parameter and try again to initialize "
                "the pipeline ...!"
            )
            msg.setStandardButtons(QMessageBox.Close)
            msg.buttonClicked.connect(msg.close)
            msg.exec()
            return False

        if (self.sess_cond_pmod_names) != [[[]]]:
            if self.sess_cond_pmod_names[idx_session][idx_cond]:
                pmods = []

                for i, val in enumerate(
                    self.sess_cond_pmod_names[idx_session][idx_cond]
                ):
                    pmod = dict()
                    pmod["name"] = val
                    param = []

                    for j, val in enumerate(
                        self.sess_cond_pmod_values[idx_session][idx_cond][i]
                    ):
                        param.insert(j, np.float64(val))

                    pmod["param"] = param
                    pmod["poly"] = self.sess_cond_pmod_polys[idx_session][
                        idx_cond
                    ][i]
                    pmods.append(pmod)

                cond["pmod"] = pmods

        try:
            cond["orth"] = self.sess_cond_orth[idx_session][idx_cond]

        except IndexError:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle("Level1Design Error!")
            msg.setText(
                "Warning: It seems that not all sess_cond_orth "
                "parameter values are defined for all conditions. "
                "Each condition must have a value for the "
                "sess_cond_orth parameter. Please, check the "
                "sess_cond_orth parameter and try again to initialize "
                "the pipeline ...!"
            )
            msg.setStandardButtons(QMessageBox.Close)
            msg.buttonClicked.connect(msg.close)
            msg.exec()
            return False

        return cond

    def _get_session_info(self, idx_session):
        """Generate the session_info parameter for each session.

        The 'session_info' dictonary is equivalent to the 'sess'
        structure array spm parameter, for each session. It therefore includes
        the 'scans', 'cond', 'multi', 'regress', 'multi_reg' and 'hpf' keys,
        like for the 'sess' parameter of spm.

        :param idx_session: session index.
        :returns: session_info, beta_sess parameters and condition number
        """
        session_info = dict()  # The information for this session
        beta_sess = 0  # The regressors number for this session
        cond_nb = 0  # The condition number for this session

        # scans
        session_info["scans"] = self.sess_scans[idx_session]

        # When HFR used, get model derivatives to multiply number of beta
        multiplier = 1
        if "hrf" in self.bases.keys():
            if "derivs" in self.bases["hrf"].keys():
                multiplier += sum(self.bases["hrf"]["derivs"])

        # cond
        if self.sess_cond_names[idx_session]:
            cond = []

            # idx_cond: condition index
            for idx_cond in range(len(self.sess_cond_names[idx_session])):
                condition = self._get_conditions(idx_session, idx_cond)
                cond.append(condition)
                beta_sess += 1 * multiplier
                cond_nb += 1

                if "pmod" in condition.keys():
                    for i in condition["pmod"]:
                        beta_sess += i["poly"]

                if "tmod" in condition.keys():
                    beta_sess += condition["tmod"]

            session_info["cond"] = cond

        # multi
        if self.sess_multi and self.sess_multi[idx_session]:
            session_info["multi"] = {self.sess_multi[idx_session]}
            mat = scipy.io.loadmat(
                self.sess_multi[idx_session], squeeze_me=True
            )

            if "names" in mat:
                beta_sess += mat["names"].shape[0]
                cond_nb += mat["names"].shape[0]

            if "pmod" in mat:
                for i in range(len(mat["pmod"])):
                    if isinstance(mat["pmod"][i][2], int):
                        beta_sess += mat["pmod"][i][2]

                    elif isinstance(mat["pmod"][i][2], np.ndarray):
                        for j in range(len(mat["pmod"][i][2])):
                            beta_sess += mat["pmod"][i][2][j]

            if "tmod" in mat:
                for i in range(len(mat["tmod"])):
                    beta_sess += mat["tmod"][i]

        # regress
        if (
            (self.sess_regress not in ["<undefined>", Undefined])
            and (not all([not elem for elem in self.sess_regress]))
            and (self.sess_regress[idx_session])
        ):
            for regressor in self.sess_regress[idx_session]:
                val = []

                for j, value in enumerate(regressor["val"]):
                    val.insert(j, np.float64(value))

                regressor["val"] = val
                beta_sess += 1

            session_info["regress"] = self.sess_regress[idx_session]

        # multi_reg
        if (
            (self.sess_multi_reg not in ["<undefined>", Undefined])
            and (not all([not elem for elem in self.sess_multi_reg]))
            and (self.sess_multi_reg[idx_session])
        ):
            session_info["multi_reg"] = []

            for reg_file in self.sess_multi_reg[idx_session]:
                if os.path.basename(reg_file)[0:3] == "rp_":
                    beta_sess += 6

                if os.path.splitext(reg_file)[1] == ".mat":
                    mat = scipy.io.loadmat(reg_file)
                    beta_sess += mat["R"].shape[1]

                reg_to_add = [{reg_file}]
                session_info["multi_reg"].append(reg_to_add)

        # hpf
        if self.sess_hpf[idx_session]:
            session_info["hpf"] = self.sess_hpf[idx_session]

        beta_sess += 1
        return session_info, beta_sess, cond_nb

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dict),
        The optional self.inheritance_dict can have two structures. On the one
        hand it can be a dictionary whose keys are the documents to inherit
        metadata and the values the documents used for the inheritance of
        these metadata. On the other hand, it can be a dictionary whose keys
        are the documents to inherit metadata and the values are dictionary
        with two keys; parent (for usual inheritance) or own_tags (to add a
        new tag or modify an existing one). In order not to include an output
        in the database, the name of the plug related to this output must be
        an element of the list corresponding to the value of the optional key
        "notInDb" of the self.outputs dictionary. To work properly this method
        must return self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Level1Design, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if (
            (self.sess_scans)
            and (self.sess_scans not in ["<undefined>", Undefined])
            and (self.sess_scans[0] not in ["<undefined>", Undefined])
        ):
            # The management of self.process.output_directory could be
            # delegated to the
            # populse_mia.user_interface.pipeline_manager.process_mia
            # module. We can't do it at the moment because the
            # sync_process_output_traits() of the
            # capsul/process/nipype_process module raises an exception
            # in nipype if the mandatory parameter are not yet defined!

            dir_name = ""
            subjects_names = []
            for idx_session in range(len(self.sess_scans)):
                sub_name = get_dbFieldValue(
                    self.project, self.sess_scans[idx_session], "PatientName"
                )
                if sub_name is None:
                    print("Please, fill 'PatientName' tag " "in the database")
                    return self.make_initResult()
                if sub_name not in subjects_names:
                    subjects_names.append(sub_name)
                    dir_name += sub_name + "_"

            if self.output_directory:
                # Create a directory for this analysis
                out_directory = os.path.join(
                    self.output_directory, dir_name + "data"
                )
                if not os.path.exists(out_directory):
                    os.mkdir(out_directory)
                self.process.output_directory = out_directory
                self.dict4runtime["out_directory"] = out_directory
            else:
                print("No output_directory was found...!\n")
                return self.make_initResult()

            self.outputs["spm_mat_file"] = os.path.join(
                out_directory, "SPM.mat"
            )

            sessions = []  # The total session_info parameter for nipype
            beta = 0  # The total number of regressors in the GLM design matrix
            cond_totNb = 0  # The total number of condition in cond and multi

            # idx_session: session index
            for idx_session in range(len(self.sess_scans)):
                # session_info: All parameters by session (Scans, Conditions,
                # Multiple conditions, Regressors, Multiple regressors,
                # High-pass filter)
                # beta_sess: Number of regressors by session
                session_info, beta_sess, cond_nb = self._get_session_info(
                    idx_session
                )
                sessions.append(session_info)
                beta += beta_sess
                cond_totNb += cond_nb

            if self.volterra_expansion_order == 2:
                beta += int((cond_totNb + 1) * cond_totNb / 2)

            # self.sessions = sessions  # The session_info for nipype
            self.dict4runtime["sessions"] = sessions

            # Some tests to check that the definition of the parameters for
            # the self.sessions went well:
            # - If self.sessions have no key in ['cond', 'multi', 'regress',
            #   'multi_reg'], or a condition is False, we can think there is
            #   an initialisation issue ...
            check = True  # Flag to stop the check processes
            init_res = []  # Initialisation result; True: Ok, False: Fail

            for i in sessions:
                for j in ["cond", "multi", "regress", "multi_reg"]:
                    if j in i.keys():
                        init_res.append(True)

                if i.get("cond") is not None and i.get("cond") is False:
                    init_res.append(False)

                if (True not in init_res or False in init_res) and (check):
                    check = False
                    self.outputs = {}
                    print(
                        "\nThere seems to be a problem in the definition of "
                        "session parameters. The initialisation failed, "
                        "please, check your settings  ..."
                    )

            # - If sess_multi_reg is plugged and sessions have no multi_reg
            #   key, we can think there is an initialisation issue ...

            if self.sess_multi_reg != [[]] and check:
                init_res = []
                [
                    init_res.append(True)
                    for i in sessions
                    if "multi_reg" in i.keys()
                ]

                if not any(init_res):
                    self.outputs = {}
                    print(
                        "\nThe sess_multi_reg plug is linked to a node. "
                        "However, no parameters for multi_reg have been "
                        "defined in the nipype session_info during "
                        "initialisation. This leads to an initialisation "
                        "failure. Please, unplug the sess_multi_reg plug if "
                        "not needed or check your settings ..."
                    )

            # - If sess_multi is plugged and sessions have no multi key,
            #   we can think there is an initialisation issue ...

            # if is_plugged['sess_multi'] and check:
            if self.sess_multi != [] and check:
                init_res = []
                [
                    init_res.append(True)
                    for i in sessions
                    if "multi" in i.keys()
                ]

                if not any(init_res):
                    self.outputs = {}
                    print(
                        "\nThe sess_multi plug is linked to a node. However, "
                        "no parameters for sess_multi have been defined in "
                        "the nipype session_info during initialisation. This "
                        "leads to an initialisation failure. Please, unplug "
                        "the sess_multi plug if not needed or check your "
                        "settings ..."
                    )

        if self.outputs:
            self.inheritance_dict[self.outputs["spm_mat_file"]] = dict()
            # FIXME: Currently, spm_mat_file will only inherit the first scan
            #        if there are several scans in self.sess_scans. This
            #        requires some thought on how to operate in a more general
            #        framework
            self.inheritance_dict[self.outputs["spm_mat_file"]][
                "parent"
            ] = self.sess_scans[0]
            self.inheritance_dict[self.outputs["spm_mat_file"]][
                "own_tags"
            ] = []
            # Add tag for number of regressors
            tag_to_add = dict()
            tag_to_add["name"] = "Regress num"
            tag_to_add["field_type"] = FIELD_TYPE_INTEGER
            tag_to_add["description"] = "Total number of regressors"
            tag_to_add["visibility"] = True
            tag_to_add["origin"] = TAG_ORIGIN_USER
            tag_to_add["unit"] = None
            tag_to_add["default_value"] = None
            tag_to_add["value"] = beta
            self.inheritance_dict[self.outputs["spm_mat_file"]][
                "own_tags"
            ].append(tag_to_add)
            # FIXME: In the latest version of mia, indexing of the
            #        database with particular tags defined in the
            #        processes is done only at the end of the
            #        initialisation of the whole pipeline. So we
            #        cannot use the value of these tags in other
            #        processes of the pipeline at the time of
            #        initialisation (see populse_mia #290). Unti
            #        better we use a quick and dirty hack with the
            #        set_dbFieldValue() function !
            set_dbFieldValue(
                self.project, self.outputs["spm_mat_file"], tag_to_add
            )

            # if self.factor_info:
            #     # Add tag for factorial design
            #     tag_to_add = dict()
            #     tag_to_add["name"] = "Factorial design info"
            #     tag_to_add["field_type"] = FIELD_TYPE_LIST_JSON
            #     tag_to_add["description"] = ("Level1Design factorial "
            #                                  "design info")
            #     tag_to_add["visibility"] = True
            #     tag_to_add["origin"] = TAG_ORIGIN_USER
            #     tag_to_add["unit"] = None
            #     tag_to_add["default_value"] = None
            #     tag_to_add["value"] = self.factor_info
            #     self.inheritance_dict[self.outputs["spm_mat_file"]][
            #         "own_tags"
            #     ].append(tag_to_add)
            #     # FIXME: In the latest version of mia, indexing of the
            #     #        database with particular tags defined in the
            #     #        processes is done only at the end of the
            #     #        initialisation of the whole pipeline. So we
            #     #        cannot use the value of these tags in other
            #     #        processes of the pipeline at the time of
            #     #        initialisation (see populse_mia #290). Unti
            #     #        better we use a quick and dirty hack with the
            #     #        set_dbFieldValue() function !
            #     set_dbFieldValue(
            #         self.project, self.outputs["spm_mat_file"], tag_to_add
            #     )

            dyn_num = 0

            for scan in self.sess_scans:
                patient_name = get_dbFieldValue(
                    self.project, scan, "PatientName"
                )

                if patient_name is not None:
                    tag_to_add = dict()
                    tag_to_add["name"] = "PatientName"
                    tag_to_add["field_type"] = FIELD_TYPE_STRING
                    tag_to_add["description"] = ""
                    tag_to_add["visibility"] = True
                    tag_to_add["origin"] = TAG_ORIGIN_USER
                    tag_to_add["unit"] = None
                    tag_to_add["default_value"] = None
                    tag_to_add["value"] = patient_name
                    self.inheritance_dict[self.outputs["spm_mat_file"]][
                        "own_tags"
                    ].append(tag_to_add)
                    # FIXME: In the latest version of mia, indexing of the
                    #        database with particular tags defined in the
                    #        processes is done only at the end of the
                    #        initialisation of the whole pipeline. So we
                    #        cannot use the value of these tags in other
                    #        processes of the pipeline at the time of
                    #        initialisation (see populse_mia #290). Unti
                    #        better we use a quick and dirty hack with the
                    #        set_dbFieldValue() function !
                    set_dbFieldValue(
                        self.project, self.outputs["spm_mat_file"], tag_to_add
                    )

                else:
                    print(
                        "\nLevel1Design:\n The PatientName tag is not filled "
                        "in the database for the {} file ...\nThis may cause "
                        "issues in the further operation of the "
                        "pipeline...\n".format(scan)
                    )

                age = get_dbFieldValue(self.project, scan, "Age")

                if age is not None:
                    tag_to_add = dict()
                    tag_to_add["name"] = "Age"
                    tag_to_add["field_type"] = "int"
                    tag_to_add["description"] = ""
                    tag_to_add["visibility"] = True
                    tag_to_add["origin"] = "user"
                    tag_to_add["unit"] = None
                    tag_to_add["default_value"] = None
                    tag_to_add["value"] = age
                    set_dbFieldValue(
                        self.project, self.outputs["spm_mat_file"], tag_to_add
                    )

                pathology = get_dbFieldValue(self.project, scan, "Pathology")

                if pathology is not None:
                    tag_to_add = dict()
                    tag_to_add["name"] = "Pathology"
                    tag_to_add["field_type"] = "string"
                    tag_to_add["description"] = ""
                    tag_to_add["visibility"] = True
                    tag_to_add["origin"] = "user"
                    tag_to_add["unit"] = None
                    tag_to_add["default_value"] = None
                    tag_to_add["value"] = pathology
                    set_dbFieldValue(
                        self.project, self.outputs["spm_mat_file"], tag_to_add
                    )

                dimensions = get_dbFieldValue(
                    self.project,
                    scan,
                    "Dataset dimensions (Count, X,Y,Z,T...)",
                )

                if (
                    (dimensions is not None)
                    and (isinstance(dimensions, list))
                    and (len(dimensions) == 5)
                ):
                    dyn_num += dimensions[4]
                    tag_to_add = dict()
                    tag_to_add["name"] = "Dynamic Number"
                    tag_to_add["field_type"] = FIELD_TYPE_INTEGER
                    tag_to_add["description"] = (
                        "Total number of dynamics in " "the functionals"
                    )
                    tag_to_add["visibility"] = True
                    tag_to_add["origin"] = TAG_ORIGIN_USER
                    tag_to_add["unit"] = None
                    tag_to_add["default_value"] = None
                    tag_to_add["value"] = dyn_num
                    self.inheritance_dict[self.outputs["spm_mat_file"]][
                        "own_tags"
                    ].append(tag_to_add)
                    # FIXME: In the latest version of mia, indexing of the
                    #        database with particular tags defined in the
                    #        processes is done only at the end of the
                    #        initialisation of the whole pipeline. So we
                    #        cannot use the value of these tags in other
                    #        processes of the pipeline at the time of
                    #        initialisation (see populse_mia #290). Unti
                    #        better we use a quick and dirty hack with the
                    #        set_dbFieldValue() function !
                    set_dbFieldValue(
                        self.project, self.outputs["spm_mat_file"], tag_to_add
                    )

                else:
                    print(
                        '\nLevel1Design:\nThe "dynamics number" tag is not '
                        "filled in the database for the {} file ...\nThis "
                        "may cause issues in the further operation of the "
                        "pipeline...\n".format(scan)
                    )

            if (self.interscan_interval is Undefined) and (
                "RepetitionTime"
                in self.project.session.get_fields_names(COLLECTION_CURRENT)
            ):
                # FIXME: Currently, spm_mat_file will only inherit the first
                #        scan if there are several scans in self.sess_scans.
                #        This requires some thought on how to operate in a
                #        more general framework
                rep_time = get_dbFieldValue(
                    self.project, self.sess_scans[0], "RepetitionTime"
                )
                if rep_time is not None:
                    self.interscan_interval = rep_time[0] / 1000

                else:
                    self.outputs = {}
                    print(
                        "\nLevel1Design:\n The interscan_interval parameter "
                        "could not be determined automatically because the "
                        '"repetition time" tag for the {} file is not filled '
                        "in the database. Please set this value and launch "
                        "the calculation again!\n".format(self.sess_scans[0])
                    )

            elif self.interscan_interval is Undefined:
                self.outputs = {}
                print(
                    "\nLevel1Design:\nThe interscan_interval parameter "
                    "(repetition time in seconds) could not be determined "
                    "automatically. Please set this value and launch the "
                    "calculation again!\n"
                )

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Level1Design, self).run_process_mia()
        # Removing the spm_mat_file to avoid a bug (nipy/nipype Issues #2612)
        if self.output_directory:
            out_file = os.path.join(
                self.dict4runtime["out_directory"], "SPM.mat"
            )

            if os.path.isfile(out_file):
                os.remove(out_file)
        else:
            print("No output_directory was found...!\n")
        self.process.spm_mat_dir = self.dict4runtime["out_directory"]
        self.process.timing_units = self.timing_units
        self.process.interscan_interval = self.interscan_interval
        self.process.microtime_resolution = self.microtime_resolution
        self.process.microtime_onset = self.microtime_onset
        self.process.session_info = self.dict4runtime["sessions"]
        self.process.factor_info = self.factor_info
        self.process.bases = self.bases
        self.process.volterra_expansion_order = self.volterra_expansion_order
        (
            self.process.global_intensity_normalization
        ) = self.global_intensity_normalization
        self.process.mask_threshold = self.mask_threshold

        # Only one mask can be defined in spm. If more than one is given, only
        # the first one will be used for all sessions ...
        if self.mask_image != Undefined:
            if type(self.mask_image) in [
                list,
                traits.TraitListObject,
                traits.List,
            ]:
                self.process.mask_image = os.path.abspath(self.mask_image[0])

            else:
                self.process.mask_image = os.path.abspath(self.mask_image)

        self.process.model_serial_correlations = self.model_serial_correlations

        return self.process.run(configuration_dict={})
