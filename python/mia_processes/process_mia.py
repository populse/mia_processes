##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

import os
from traits.api import Undefined

# Populse_MIA imports
from populse_mia.pipeline_manager.process_mia import ProcessMIA


class Process_Mia(ProcessMIA):
    """Class overriding the ProcessMIA class, in order to personalize the run in MIA

        Methods:
            - _run_processes: overrides the ProcessMIA's _run_process that is called to run the process
            - run_process_mia: (need to be overridden)
            - switch_to_scripts_dir:  Changes the current working directory to the scripts directory
            - switch_to_cur_work_dir: Changes the scripts directory to the current working directory
    """
   
    def __init__(self):
        super(Process_Mia, self).__init__()
        self.change_dir = False

    def _run_process(self):
        """
        Method overriding the ProcessMIA's _run_process that is called to run the process
        """
        if self.change_dir:
            self.switch_to_scripts_dir()

        self.run_process_mia()

        if self.change_dir:
            self.switch_to_cur_work_dir()

    def run_process_mia(self):
        if self.change_dir:
            self.manage_matlab_launch_parameters()

    def switch_to_scripts_dir(self):
        """Method that changes the current working directory to the scripts directory.

           Returns
           -------
           cwd: The current working directory"""

        try:
            self.cwd = os.getcwd()

        except OSError:
            self.cwd = None

        if not hasattr(self, 'output_directory') or self.output_directory is None or self.output_directory is Undefined:
            raise ValueError('output_directory is not set but is mandatory to run a Process_Mia')

        print('\nChanging from {0} directory to {1} directory ...\n'.format(self.cwd, self.output_directory))
        os.chdir(self.output_directory)
        
    def switch_to_cur_work_dir(self):
        """Method that changes the scripts directory to the current working directory."""
        try:
            cwd1 = os.getcwd()

        except OSError:
            cwd1 = None

        try:
            os.chdir(self.cwd)
            print('Changing from {0} directory to {1} directory ...\n'.format(cwd1, self.cwd))

        except Exception as e:
            print('{0}: {1}'.format(e.__class__, e))
