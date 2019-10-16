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

import os
from traits.api import Undefined
import subprocess

# Populse_MIA imports
from mia_processes.process_mia import Process_Mia
from populse_mia.software_properties import Config


class ProcessMatlab(Process_Mia):
    """Class overriding the Process_MIA class, in order to use Matlab.

    Attributes:
        - project: current project in the software

    Methods:

    """

    
    def __init__(self):
        super(ProcessMatlab, self).__init__()

        self.use_project = True
        self.matlab_script = ""

    def set_variable(self, variable_name, value):
        """Adds "variable_name = value;" to the Matlab script.

        :param variable_name: name of the variable
        :param value: value of the variable
        """
        if type(value) is str:
            self.matlab_script += '{0} = "{1}";'.format(variable_name, value)
        else:
            self.matlab_script += '{0} = {1};'.format(variable_name, value)

    def set_global_variable(self, variable_name):
        """
        Adds a global variable to the Matlab script
        :param variable_name: name of the variable
        :return:
        """
        self.matlab_script += 'global {0};'.format(variable_name)

    def change_directory(self, directory):
        """Changes the working directory.

        :param directory: directory
        """
        self.matlab_script += 'cd("{0}");'.format(directory)

    def display_parameter(self, parameter_name):
        """Displays a given parameter.

        :param parameter_name: name of the parameter
        """
        self.matlab_script += 'disp({0});'.format(parameter_name)

    def add_path(self, path):
        """Adds a Matlab path.

        :param path: path

        """
        
        self.matlab_script += 'addpath("{0}");'.format(path)

    def add_exit(self):
        """Adds an exit to the Matlab script."""
        self.matlab_script += 'exit'

    def run(self):
        """Runs the Matlab script."""
        config = Config()
        subprocess.run([config.get_matlab_path(), '-nodisplay', '-r',
                        self.matlab_script])

