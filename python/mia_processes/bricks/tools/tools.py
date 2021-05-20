# -*- coding: utf-8 -*- #

"""The toolbox library of the mia_processes package.

Basically, this module is dedicated to the low-level processes
needed to run other higher-level bricks.

:Contains:
    :Class:
        - Files_To_List
        - Filter_Files_List
        - Input_Filter
        - List_Duplicate
        - List_To_File

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# nipype import
from nipype.interfaces.base import traits

# populse_mia imports
from populse_mia.data_manager.filter import Filter
from populse_mia.data_manager.project import COLLECTION_CURRENT
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# Other imports
import os


class Files_To_List(ProcessMIA):
    """
    *From 2 file names, generating a list containing all theses file names*

    Please, see the complete documention for the `Files_To_List in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/tools/Files_To_List.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.
        
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
                                   optional=True,
                                   desc=file2_desc))

        # Outputs traits
        self.add_trait("file_list",
                       traits.List(output=True,
                                   desc=file_list_desc))

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
        super(Files_To_List, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.file1 and self.file1 not in ["<undefined>", traits.Undefined]:
        
            if ((not self.file2) or
                (self.file2 in ["<undefined>", traits.Undefined])):
                self.outputs['file_list'] = [self.file1]

            else:
                self.outputs['file_list'] = [self.file1, self.file2]
                
            if self.outputs:
                self.outputs["notInDb"] = ["file_list"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class Filter_Files_List(ProcessMIA):
    """
    *Selects one or more (slicing) element(s) from a list*

    Please, see the complete documention for the `Filter_Files_List in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/tools/Filter_Files_List.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Filter_Files_List, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = [] # no need of third party software!

        # Inputs description
        in_list_desc = 'The list of elements to be filtered.'
        index_filter_desc = 'A list of 0 to 2 indexes for filtering.'

        # Outputs description
        filtered_list_desc = 'The corresponding filtering result (a list).'
        
        # Inputs traits
        self.add_trait("in_list",
                        traits.List(traits.File,
                                    output=False,
                                    desc=in_list_desc))
        
        self.add_trait("index_filter",
                       traits.List(value=[1],
                                   trait=traits.Range(low=1, high=None),
                                   minlen=0,
                                   maxlen=2,
                                   output=False,
                                   optional=True,
                                   desc=index_filter_desc))

        # Outputs traits
        self.add_trait("filtered_list",
                       traits.List(output=True,
                                   desc=filtered_list_desc))

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
        super(Filter_Files_List, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if (self.in_list and
               self.in_list not in ["<undefined>", traits.Undefined]):

            if not self.index_filter:
                self.outputs['filtered_list'] = [self.in_list[0]]

            if len(self.index_filter) == 1:

                if self.index_filter[0] <= len(self.in_list):
                    self.outputs['filtered_list'] = [self.in_list[
                                                        self.index_filter[0]-1]]

                else:
                    print('\nThe initialisation of the Filter_Files_List brick '
                          'failed because the index_filter parameter is '
                          'greater than the length of the in_list '
                          'parameter ...\n')

            if len(self.index_filter) == 2:

                if self.index_filter[0] < self.index_filter[1]:

                    if self.index_filter[0] <= len(self.in_list):

                        if self.index_filter[1] <= len(self.in_list):
                            self.outputs['filtered_list'] = self.in_list[
                                    self.index_filter[0]-1:self.index_filter[1]]

                        else:
                            print('\nThe initialisation of the '
                                  'Filter_Files_List brick failed because the '
                                  'second value of the index_filter parameter '
                                  'is greater than the length of the in_list '
                                  'parameter ...\n')

                    else:
                        print('\nThe initialisation of the Filter_Files_List '
                              'brick failed because the first value of the '
                              'index_filter parameter is greater than the '
                              'length of the in_list parameter ...\n')

                else:
                     print('\nThe initialisation of the Filter_Files_List brick'
                           ' failed because the first value of the index_filter'
                           ' parameter is greater than the second ...\n')

            if self.outputs:
                self.outputs["notInDb"] = ["filtered_list"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class Input_Filter(ProcessMIA):
    """
    *To filter the content of the Data Browser tab or the output data of another brick*

    Please, see the complete documention for the `Input_Filter in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/tools/Input_Filter.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Input_Filter, self).__init__()
        
        # Third party softwares required for the execution of the brick
        self.requirement = [] # no need of third party software!

        # Inputs description
        input_desc = ('A list corresponding to the data from the Data Browser '
                      'or the output data from another brick.')

        # Outputs description
        output_desc = 'A list with the result of the filter applied.'

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
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key 'notInDb' of the
        self.outputs dictionary. To work properly this method must return 
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
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
        # database_scans, or not) and of course to perform a filtering from it.
        # But also, if no data is sent in, we get all the data from the database
        # (self.project.session.get_documents_names(COLLECTION_CURRENT), below).
        # It can maybe lead to side effects (however it also allows for example
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
        self.outputs["notInDb"] = ["output"]

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


class List_Duplicate(ProcessMIA):
    """
    *From a file name, generating a list containing this file name and the file name itself*

    Please, see the complete documention for the `List_Duplicate in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/tools/List_Duplicate.html>`_

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
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key 'notInDb' of the
        self.outputs dictionary. To work properly this method must return 
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(List_Duplicate, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.file_name:
            self.outputs["out_list"] = [self.file_name]
            self.outputs["out_file"] =  self.file_name
            self.outputs["notInDb"] = ["out_list", "out_file"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class List_To_File(ProcessMIA):
    """
    *From several filenames, selects and generates a file.*

    Please, see the complete documention for the `List_To_File in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/tools/List_To_File.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(List_To_File, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = [] # no need of third party software!

        # Inputs description
        file_list_desc = 'The list of elements to be filtered.'
        index_filter_desc = 'A list of 0 to 1 indexes for filtering.'

        # Outputs description
        file_desc = 'The corresponding filtering result (a file).'
        
        # Inputs traits
        self.add_trait("file_list",
                       traits.List(traits.File(),
                                   output=False,
                                   desc=file_list_desc))
        
        self.add_trait("index_filter",
                       traits.List(value=[1],
                                   trait=traits.Range(low=1, high=None),
                                   minlen=0,
                                   maxlen=1,
                                   output=False,
                                   optional=True,
                                   desc=index_filter_desc))

        # Outputs traits
        self.add_trait("file",
                       traits.File(output=True,
                                   desc=file_desc))

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
        super(List_To_File, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if (self.file_list and
               self.file_list not in ["<undefined>", traits.Undefined]):

            if not self.index_filter:
                self.outputs['file'] = [self.file_list[0]]

            if len(self.index_filter) == 1:

                if self.index_filter[0] <= len(self.file_list):
                    self.outputs['file'] = self.file_list[
                                                         self.index_filter[0]-1]

                else:
                    print('\nThe initialisation of the List_To_File brick '
                          'failed because the index_filter parameter is '
                          'greater than the length of file_list '
                          'parameter ...\n')
        
        if self.outputs:
            self.outputs["notInDb"] = ["file"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return
