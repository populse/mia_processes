# -*- coding: utf-8 -*- #

"""The toobox to run a brick using Matlab.

Basically, this module provides the necessary tools for
the launch of the bricks using Matlab.

:Contains:
    :Class:
        - ProcessMatlab

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
from populse_mia.software_properties import Config

# Other imports
import subprocess


class ProcessMatlab(ProcessMIA):
    """Class overriding the ProcessMIA class, in order to use Matlab.

    .. Methods:
        - add_exit: Add an exit to the Matlab script
        - add_path: Add a Matlab path to the Matlab script
        - change_directory: Change the working directory in the Matlab script
        - display_parameter: Display a given parameter in the Matlab script
        - run: Run the Matlab script
        - set_global_variable: Add a global variable to the Matlab script
        - set_variable: "Assign a value to a variable in the Matlab script
    """

    def __init__(self):
        super(ProcessMatlab, self).__init__()
        self.use_project = True
        self.matlab_script = ""

    def add_exit(self):
        """Adds an exit to the Matlab script."""
        self.matlab_script += 'exit'

    def add_path(self, path):
        """Adds a Matlab path to the Matlab script.

        :param path: Matlab path
        """
        self.matlab_script += 'addpath("{0}");'.format(path)

    def change_directory(self, directory):
        """Changes the working directory in the Matlab script.

        :param directory: directory
        """
        self.matlab_script += 'cd("{0}");'.format(directory)

    def display_parameter(self, parameter_name):
        """Displays a given parameter in the Matlab script.

        :param parameter_name: name of the parameter
        """
        self.matlab_script += 'disp({0});'.format(parameter_name)

    def run(self):
        """Runs the Matlab script."""
        config = Config()
        subprocess.run([config.get_matlab_path(), '-nodisplay', '-r',
                        self.matlab_script])

    def set_global_variable(self, variable_name):
        """Adds a global variable to the Matlab script.

        :param variable_name: name of the variable
        """
        self.matlab_script += 'global {0};'.format(variable_name)

    def set_variable(self, variable_name, value):
        """Assign a value to a variable in the Matlab script.

        :param variable_name: name of the variable
        :param value: value of the variable
        """
        if type(value) is str:
            self.matlab_script += '{0} = "{1}";'.format(variable_name, value)
        else:
            self.matlab_script += '{0} = {1};'.format(variable_name, value)
