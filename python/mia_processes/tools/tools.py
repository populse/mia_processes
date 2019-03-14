##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

import os

# Trait import
from nipype.interfaces.base import traits

# MIA import
#from populse_mia.project.filter import Filter
from mia_processes.process_mia import Process_Mia

                        ### Bricks/classes in this file: ###
                        # - Files_To_List                  #
                        ####################################
                        ### No methods in this files     ###
                        ####################################

class Files_To_List(Process_Mia):
    """
-  Files_To_List (mia_processes.tools.tools.Files_To_List)
*** From 2 file names, generating a list containing all theses file names ***
    '/home/ArthurBlair/data/raw_data/Anat.nii ' + '/home/ArthurBlair/data/raw_data/Func.nii' -> Files_To_List ->
    ['/home/ArthurBlair/data/raw_data/Anat.nii', '/home/ArthurBlair/data/raw_data/Func.nii']
    * Input parameters:
        * file1: a string corresponding to an existing path file (traits.File)
            <ex. /home/ArthurBlair/data/raw_data/Anat.nii>
        * file2: a string corresponding to an existing path file (traits.File)
            <ex. /home/ArthurBlair/data/raw_data/Func.nii>
    * Output parameters:
        * file_list: a list (traits.List)
            <ex. ['/home/ArthurBlair/data/raw_data/Anat.nii',
                 '/home/ArthurBlair/data/Func.nii']>
    """

    def __init__(self):
        super(Files_To_List, self).__init__()

        # Inputs description
        file1_desc = 'A string corresponding to an existing path file.'
        file2_desc = 'A string corresponding to an existing path file.'

        # Outputs description
        file_list_desc = 'A list of items which are an existing file name.'
        
        # Inputs traits
        self.add_trait("file1",
                       traits.File(output=False,
                                   desc=file1_desc))
        
        self.add_trait("file2",
                       traits.File(output=False,
                                   desc=file2_desc))

        # Outputs traits
        self.add_trait("file_list",
                       traits.List(output=True,
                                   desc=file_list_desc))

    def list_outputs(self):
        super(Files_To_List, self).list_outputs()
        out_list = [self.file1, self.file2]
        return {'file_list': out_list}, {}

    def run_process_mia(self):
        out_list = [self.file1, self.file2]
        self.file_list = out_list
