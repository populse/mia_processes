# -*- coding: utf-8 -*- #

"""The toobox to personalise the run in MIA.

Basically, this module provides the necessary tools for
a custom execution of the bricks in populse_mia.

:Contains:
    :Class:
        - Process_Mia

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# Populse_MIA imports
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# Soma-base imports
from soma.controller.trait_utils import relax_exists_constraint

# Other imports
import os
from traits.api import Undefined
import traits.api as traits

from capsul.api import Process

# class Process_Mia(Process): Ne marche pas avoir !
class Process_Mia(ProcessMIA):
    """Class overriding the ProcessMIA class, in order to personalise 
       the run in MIA.

        Methods:
            - _run_processes: call the run_process_mia method in the 
               Process_Mia subclass
            - list_outputs: override the outputs of the process
            - make_initResult: make the final dictionnary for outputs,
                               inheritance and requirement from the
                               initialisation of a brick
            - relax_nipype_exists_constraints: relax the exists constraint of
                                               the process.inputs traits
            - requirements: capsul Process.requirements() implementation using
                            MIA's Process_Mia.requirement attribute
            - run_process_mia: implements specific runs for Process_Mia
                               subclasses
    
    """
    
    def __init__(self, *args, **kwargs):
        super(Process_Mia, self).__init__(*args, **kwargs)
        self.change_dir = False
        self.requirement = None
        self.outputs = {}
        self.inheritance_dict = {}

    def _run_process(self):
        """ call the run_process_mia method in the Process_Mia subclass"""
        self.run_process_mia()

    def list_outputs(self):
        """Override the outputs of the process."""
        self.relax_nipype_exists_constraints()

    def make_initResult(self):
        """Make the initResult_dict from initialisation."""        
        if ((self.requirement is None) or
            (not self.inheritance_dict) or
            (not self.outputs)):
            print('\nDuring the {0} process initialisation, some possible '
                   'problems were detected:'.format(self))
             
            if self.requirement is None:
                 print('- requirement attribute was not found ...')

            if not self.inheritance_dict:
                print('- inheritance_dict attribute was not found ...')

            if not self.outputs:
                print('- outputs attribute was not found ...')

            print()

        return {'requirement': self.requirement, 'outputs': self.outputs,
                'inheritance_dict': self.inheritance_dict}

    def relax_nipype_exists_constraints(self):
        """Relax the exists constraint of the process.inputs traits"""
        if hasattr(self, 'process') and hasattr(self.process, 'inputs'):
            ni_inputs = self.process.inputs
            for name, trait in ni_inputs.traits().items():
                relax_exists_constraint(trait)

    def requirements(self):
        """Capsul Process.requirements() implementation using MIA's
        Process_Mia.requirement attribute
        """
        if self.requirement:
            return {req: 'any' for req in self.requirement}
        return {}

    def run_process_mia(self):
        """
        Implements specific runs for Process_Mia subclasses
        """
        pass
 
