# -*- coding: utf-8 -*- #

"""The tools library of the mia_processes package.

Basically, this module is dedicated to the low-level processes
needed to run other higher-level bricks.

:Contains:
    :Class:
        - Files_To_List
        - Input_Filter
        - List_Duplicate

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# mia_processes import
from mia_processes.process_mia import Process_Mia

# nipype import
from nipype.interfaces.base import traits

# populse_mia import
from populse_mia.data_manager.filter import Filter
from populse_mia.data_manager.project import COLLECTION_CURRENT

# Other import
import os


class Files_To_List(Process_Mia):
    """ Files_To_List (from mia_processes.tools.tools)
*** From 2 file names, generating a list containing all theses file names ***
    '/home/ArthurBlair/data/raw_data/Anat.nii '
                        + 
    '/home/ArthurBlair/data/raw_data/Func.nii'
                -> Files_To_List ->
    ['/home/ArthurBlair/data/raw_data/Anat.nii',
     '/home/ArthurBlair/data/raw_data/Func.nii']
    If only file1 is specified, returns a list containing file1 only.
    * Input parameters:
        * file1: a string corresponding to an existing path file (traits.File)
            <ex. /home/ArthurBlair/data/raw_data/Anat.nii>
        * file2: an optional string corresponding to an existing path
                 file (traits.File)
            <ex. /home/ArthurBlair/data/raw_data/Func.nii>
    * Output parameters:
        * file_list: a list (traits.List)
            <ex. ['/home/ArthurBlair/data/raw_data/Anat.nii',
                  '/home/ArthurBlair/data/Func.nii']>

    """

    
    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Files_To_List, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = [] # no need of third party software!

        # Inputs description
        file1_desc = 'A string corresponding to an existing path file.'
        file2_desc = ('An optional string corresponding '
                     'to an existing path file.')

        # Outputs description
        file_list_desc = 'A list of items which are an existing file name.'
        
        # Inputs traits
        self.add_trait("file1",
                       traits.File(output=False,
                                   desc=file1_desc))
        
        self.add_trait("file2",
                       traits.File(output=False,
                                   desc=file2_desc, optional=True))

        # Outputs traits
        self.add_trait("file_list",
                       traits.List(output=True,
                                   desc=file_list_desc))

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return 
        self.make_initResult() object.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Files_To_List, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.file1 and not self.file1 in ["<undefined>", traits.Undefined]:
        
            if ((not self.file2) or
                (self.file2 in ["<undefined>", traits.Undefined])):
                self.outputs['file_list'] = [self.file1]

            else:
                self.outputs['file_list'] = [self.file1, self.file2]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        if not self.file2 or self.file2 in ["<undefined>", traits.Undefined]:
            out_list = [self.file1]
            
        else:
            out_list = [self.file1, self.file2]
            
        self.file_list = out_list

        
class Input_Filter(Process_Mia):
    """ Input_Filter (from mia_processes.tools.tools)
*** Process (brick) used to filter the content of the Data Browser tab 
    or the output data of another brick. ***
    [Preselected data (or not) from the Data Browser
                                          or the output data from another brick]
                         ->
    ['/home/ArthurBlair/data/raw_data/Anat.nii']
    * Input parameters:
        * input: a list corresponding to the data from the Data Browser
                 or the output data from another brick (traits.List)
            <ex. ['/home/ArthurBlair/data/raw_data/Anat.nii',
                  '/home/ArthurBlair/data/Func.nii']>
    * Output parameters:
        * output: a list with the result of the filter applied (traits.List)
            <ex. ['/home/ArthurBlair/data/raw_data/Anat.nii']>
    * To run properly, this process (node, brick) needs that, on the one hand,
      in the DataBrowser, the "send documents to the Pipeline Manager" button
      was clicked previously (with a selection of data made, or not), and on
      the other hand, at the input_filter brick level, that a right click then
      the selection of the option "Export to database_scans" has been made.
      Finally, a right click on the input_filter brick will allow to filter the
      input data by selecting "open filter".

    """

    
    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Input_Filter, self).__init__()
        
        # Third party softwares required for the execution of the brick
        self.requirement = [] # no need of third party software!

        # Inputs description
        input_desc = 'Input plug description.'

        # Outputs description
        output_desc = 'Output plug description.'

        # Inputs traits
        self.add_trait("input",
                       traits.List(traits.File,
                                   output=False,
                                   desc=input_desc))

        # Outputs traits
        self.add_trait("output",
                       traits.List(traits.File,
                                   output=True,
                                   desc=output_desc))
        
        # Instantiation of the filter object
        self.filter = Filter(None, [''], [''], [['FileName']],
                             [], ['CONTAINS'], "")

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return 
        self.make_initResult() object.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Input_Filter, self).list_outputs()
        
        # TODO: MAYBE WE DON'T NEED THAT, IT SHOULD BE DONE IN
        # open_filter OF PipelineEditor?

        # Outputs definition and tags inheritance (optional)
        # Getting the input:
        # The current way for obtaining input data and filtering them is very
        # permissive. Indeed, it is possible to obtain the data from the output
        # of a previous brick or from the database (using Export to
        # database_scans) and of course to perform a filtering from it.
        # But also, if no data is sent in, we get all the data from the database
        # (self.project.session.get_documents_names(COLLECTION_CURRENT), below),
        # it can maybe lead to side effects (however it also allows for example
        # not to connect the input to something and to have the whole database
        # by default... Is it really desirable????
        if self.input: 
            self.scans_list = self.input
            
        else: 
            self.scans_list = self.project.session.get_documents_names(
                                                             COLLECTION_CURRENT)

        # The data to filter are always a relative path
        if (self.scans_list) and (os.path.isabs(self.scans_list[0])):
            self.scans_list = [os.path.relpath(i, self.project.folder)
                                                       for i in self.scans_list]

        # Apply the filter to the input
        self.outputs['output'] = self.filter.generate_filter(self.project,
                                                             self.scans_list,
                                          self.project.session.get_shown_tags())

        # The output data are always an absolute path
        for idx, element in enumerate(self.outputs['output']):
            full_path = os.path.abspath(os.path.join(self.project.folder,
                                                     element))
            self.outputs['output'][idx] = full_path
            
        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class List_Duplicate(Process_Mia):
    """List_Duplicate (from mia_processes.tools.tools)
*** From a file name, generating a list containing this file name and the
    file name itself ***
    '/home/ArthurBlair/Anat.nii'
       -> List_Duplicate ->
    ['/home/ArthurBlair/Anat.nii'] + '/home/ArthurBlair/Anat.nii'
    * Input parameters:
        * file_name: a string corresponding to an existing
                     path file (traits.File)
            <ex. /home/ArthurBlair/data/Func.nii>
    * Output parameters:
        * out_file: a string corresponding to an existing
                    path file (traits.File)
            <ex. /home/ArthurBlair/data/Func.nii>
        * out_list: a list with one string element corresponding to an
                    existing path file (traits.List)
            <ex. ['/home/ArthurBlair/data/Func.nii']>

    """

    
    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(List_Duplicate, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = [] # no need of third party software!

        # Inputs description
        file_name_desc = 'A string corresponding to an existing path file.'

        # Outputs description
        out_file_desc = 'A string corresponding to an existing path file.'
        out_list_desc = ('A list with one string element corresponding '
                         ' to an existing path file.')
        
        # Inputs traits
        self.add_trait("file_name",
                       traits.File(output=False,
                                   desc=file_name_desc))

        # Outputs traits
        self.add_trait("out_file",
                       traits.File(output=True,
                                   desc=out_file_desc))
        
        self.add_trait("out_list",
                       traits.List(output=True,
                                   desc=out_list_desc))

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return 
        self.make_initResult() object.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(List_Duplicate, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.file_name:
            self.outputs["out_list"] = [self.file_name]
            self.outputs["out_file"] =  self.file_name

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        if self.file_name:
            self.out_list = [self.file_name]
            self.out_file = self.file_name
