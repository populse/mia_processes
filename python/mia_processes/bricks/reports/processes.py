"""The other preprocess library of the mia_processe package.

The purpose of this module is to provide bricks generally necessary for the
pre-processing steps, which are not found in nipype.

:Contains:
    :Class:
        - QC_AnatReport

    :Function:
        -
"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# nibabel import
import nibabel as nib

# nipype import
from nipype.interfaces.base import (OutputMultiPath, InputMultiPath, File,
                                    traits, TraitListObject, Undefined,
                                    DictStrStr, Str)

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# soma-base imports
from soma.qt_gui.qt_backend.Qt import QMessageBox

# Other import
import os


class QC_AnatReport(ProcessMIA):
    """
        * Computes the artifact mask using the method described in [Mortamet2009] *

        Please, see the complete documentation for the `ArtifactMask' brick in the populse.mia_processes website
        https://populse.github.io/mia_processes/documentation/bricks/preprocess/other/ArtifactMask.html

        """

    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(QC_AnatReport, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description

        # Outputs description

        # Inputs traits

        # Outputs traits

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
        super(QC_AnatReport, self).list_outputs()

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(QC_AnatReport, self).run_process_mia()