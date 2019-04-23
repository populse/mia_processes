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
from populse_mia.project.filter import Filter
from mia_processes.process_mia import Process_Mia

                        ### Bricks/classes in this file: ###
                        # - Files_To_List                  #
                        # - Input_Filter                   #
                        ####################################
                        # No function in this files        #
                        ####################################

class Files_To_List(Process_Mia):
    """ Files_To_List (mia_processes.tools.tools.Files_To_List)
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

class Input_Filter(Process_Mia):
    """
-  Input_Filter (mia_processes.tools.tools.Input_Filter)
*** Process (brick) used to filter the content of the Data Browser tab, to provide an input for another brick. ***
    [Preselected data (or not) from the Data Browser] -> ['/home/ArthurBlair/data/raw_data/Anat.nii']
    * Input parameters:
        * input: a list corresponding to the preselected data (or not) from the Data Browser (traits.List)
            <ex. ['/home/ArthurBlair/data/raw_data/Anat.nii', '/home/ArthurBlair/data/Func.nii']>
    * Output parameters:
        * output: a list with the result of the filter applied (traits.List)
            <ex. ['/home/ArthurBlair/data/raw_data/Anat.nii']>
    * To run properly this process (node, brick) needs, that on the one hand, in the DataBrowser, the
      "send documents to the Pipeline Manager" button was clicked previously (with a selection of data
      made, or not) and on the other hand, at the input_filter brick level, that a right click then the
      selection of the option "Export to database_scans" has been made. Finally, a right click on the
      input_filter brick will allow to filter the input data by selecting "open filter".

    """
    
    def __init__(self):
        super(Input_Filter, self).__init__()

        self.filter = Filter(None, [''], [''], [['FileName']], [], ['CONTAINS'], "")
        self.add_trait("input", traits.List(traits.File, output=False))
        self.add_trait("output", traits.List(traits.File, output=True))

    def list_outputs(self):
        super(Input_Filter, self).list_outputs()
        # TODO: MAYBE WE DON'T NEED THAT, IT SHOULD BE DONE IN open_filter OF PipelineEditor
        
        if self.input:
            self.scans_list = self.input
        else:  # It seems that currently self.project.database have no
               # get_documents_names() method. Currently, if self.input
               # is nothing ([]) we observe this message in the terminal:
               # "No output list method for the process input_filter1" without
               # crash but the process stop after the next line ...
            self.scans_list = self.project.database.get_documents_names()
        
        filt = self.filter
        output = filt.generate_filter(self.project, self.scans_list, self.project.session.get_shown_tags())

        for idx, element in enumerate(output):
            full_path = os.path.abspath(os.path.join(self.project.folder, element))
            output[idx] = full_path

        return {'output': output}

    def run_process_mia(self):
        return
