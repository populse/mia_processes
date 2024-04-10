# -*- coding: utf-8 -*-

"""The toolbox library of the mia_processes package.

Basically, this module is dedicated to the low-level processes
needed to run other higher-level bricks.

:Contains:
    :Class:
        - Concat_to_list
        - Concat_to_list_of_list
        - Delete_data
        - Files_To_List
        - Filter_Files_List
        - Find_In_List
        - Get_Conditions_From_csv
        - Get_Patient_Name
        - Import_Data
        - Input_Filter
        - List_Duplicate
        - List_To_File
        - Make_A_List
        - Make_CVR_reg_physio

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################
# fmt: off
# fmt: on
import os
import re
import shutil
import tempfile
import traceback

# Other imports
from ast import literal_eval

import numpy as np
import pandas as pd

# nipype import
from nipype.interfaces.base import (
    File,
    InputMultiPath,
    OutputMultiPath,
    traits,
)
from nipype.interfaces.spm.base import ImageFileSPM

# populse_mia imports
from populse_mia.data_manager.data_history_inspect import (
    get_data_history_processes,
    get_filenames_in_value,
)
from populse_mia.data_manager.filter import Filter
from populse_mia.data_manager.project import BRICK_OUTPUTS, COLLECTION_CURRENT
from populse_mia.software_properties import Config
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA
from scipy import io, signal
from scipy.interpolate import interp1d
from scipy.special import gammaln
from traits.api import Either, Undefined

# mia_processes imports
from mia_processes.utils import checkFileExt, get_dbFieldValue


class Concat_to_list(ProcessMIA):
    """
    *Make an output list corresponding to the concatenation of list1 and list2*

    Ex. ['a', 'b', 'c'] and ['d', 'e'] gives
        ['a', 'b', 'c', 'd', 'e']

    Please, see the complete documentation for the `Concat_to_list brick in the
    mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/Concat_to_list.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Concat_to_list, self).__init__()

        # Inputs description
        list1_desc = "A list"
        list2_desc = "A list"

        # Outputs description
        out_list_desc = (
            "A list corresponding to the concatenation of list1 and list2"
        )

        # Inputs traits
        self.add_trait(
            "list1",
            traits.List(output=False, optional=False, desc=list1_desc),
        )
        self.list1 = traits.Undefined

        self.add_trait(
            "list2",
            traits.List(output=False, optional=False, desc=list2_desc),
        )
        self.list2 = traits.Undefined

        # Outputs traits
        self.add_trait(
            "out_list", traits.List(output=True, desc=out_list_desc)
        )
        self.out_list = traits.Undefined

        self.init_default_traits()

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
        super(Concat_to_list, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.list1 not in [[], traits.Undefined] and self.list2 not in [
            [],
            traits.Undefined,
        ]:
            self.outputs["out_list"] = self.list1 + self.list2

        if self.outputs:
            self.outputs["notInDb"] = ["out_list"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class Concat_to_list_of_list(ProcessMIA):
    """
    *Iteration of input list1 with each element of input list2*

    Ex. ['a', 'b', 'c'] and ['1', '2'] gives
        [['a', '1'], ['a', '2'],
         ['b', '1'], ['b', '2'],
         ['c', '1'], ['c', '2']

    Please, see the complete documentation for the `Concat_to_list_of_list
    brick in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/Concat_to_list_of_list.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Concat_to_list_of_list, self).__init__()

        LIST1 = [
            "ACA",
            "ACM",
            "ACP",
            "PICA",
            "ROI-CING",
            "ROI-FRON",
            "ROI-INSULA",
            "ROI-OCC",
            "ROI-PAR",
            "ROI-STR",
            "ROI-TEMP",
            "ROI-THA",
            "SCA",
        ]
        LIST2 = ["L", "R"]

        # Inputs description
        list1_desc = "A list of string"
        list2_desc = "A list of string"

        # Outputs description
        listOflist_desc = "A list of list"

        # Inputs traits
        self.add_trait(
            "list1",
            traits.List(LIST1, output=False, optional=True, desc=list1_desc),
        )

        self.add_trait(
            "list2",
            traits.List(LIST2, output=False, optional=True, desc=list2_desc),
        )

        # Outputs traits
        self.add_trait(
            "listOflist", traits.List(output=True, desc=listOflist_desc)
        )
        self.listOflist = traits.Undefined

        self.init_default_traits()

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
        super(Concat_to_list_of_list, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.list1 not in [[], traits.Undefined] and self.list2 not in [
            [],
            traits.Undefined,
        ]:
            self.outputs["listOflist"] = [
                [i, j] for i in self.list1 for j in self.list2
            ]

            if self.outputs:
                self.outputs["notInDb"] = ["listOflist"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class Delete_data(ProcessMIA):
    """
    *Delete from database data*

    Please, see the complete documentation for the `Delete_data brick in the
    mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/Delete_data.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Delete_data, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []  # no need of third party software!

        # Inputs description
        in_file_desc = (
            "Output file from a brick or a pipeline. "
            "The outputs from the history of this file "
            "will be removed."
        )
        to_keep_filters_desc = (
            "A list of regex. Files that match those "
            "regex will be kept and the others files "
            "will be deleted. "
            "Mutually exclusive with to_remove_filters"
        )
        to_remove_filters_desc = (
            "A list of regex.  Files that match those "
            "regex will be deleted and the others files "
            "will be kept. "
            "Mutually exclusive with to_keep_filters"
        )

        # Outputs description
        files_removed_desc = "List of files removed from database"

        # Inputs traits
        self.add_trait("in_file", traits.File(output=False, desc=in_file_desc))
        self.in_file = traits.Undefined

        self.add_trait(
            "to_keep_filters",
            traits.List(output=False, desc=to_keep_filters_desc),
        )
        # default for mriqc pipeline
        self.to_keep_filters = [
            "(.)*pdf",
            "(.)*_qc.json",
            "(.)*desc-carpet(.)*",
        ]

        self.add_trait(
            "to_remove_filters",
            traits.List(output=False, desc=to_remove_filters_desc),
        )
        # Outputs traits
        self.add_trait(
            "files_removed", traits.List(output=True, desc=files_removed_desc)
        )
        self.files_removed = traits.Undefined

        self.init_default_traits()

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

        super(Delete_data, self).list_outputs()

        if self.to_keep_filters and self.to_remove_filters:
            print(
                "\nDelete_data brick: Initialisation failed. to_keep_filters "
                "and to_remove_filter parameters are mutually exclusive...!"
            )
            return

        # Get history of the file
        # to get all the outputfile of the pipeline/brick
        file_position = (
            self.in_file.find(self.project.getName())
            + len(self.project.getName())
            + 1
        )
        file_database_name = self.in_file[file_position:]

        procs, _ = get_data_history_processes(file_database_name, self.project)
        outputs_filename = set()
        if procs:
            for proc in procs.values():
                if proc.used:
                    for value in proc.brick[BRICK_OUTPUTS].values():
                        filenames = get_filenames_in_value(
                            value, self.project, allow_temp=False
                        )
                        outputs_filename.update(filenames)

        outputs_filename = list(outputs_filename)
        for file in outputs_filename:
            # remove file to keep using to_keep_filter
            if self.to_keep_filters:
                for filter in self.to_keep_filters:
                    if re.match(r"" + filter, file):
                        outputs_filename.remove(file)
            # get file to remove using to_remove_filter
            if self.to_remove_filters:
                for filter in self.to_remove_filters:
                    if not re.match(r"" + filter, file):
                        outputs_filename.remove(file)

        self.files_removed = outputs_filename
        for doc in outputs_filename:
            # Remove document from database
            self.project.session.remove_document(COLLECTION_CURRENT, doc)
            full_doc_path = os.path.join(self.project.folder, doc)
            # Remove from project
            if os.path.isfile(full_doc_path):
                os.remove(full_doc_path)

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class Files_To_List(ProcessMIA):
    """
    *From 3 file names, generating a list containing all these file names*

    Please, see the complete documentation for the `Files_To_List brick in
    the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/Files_To_List.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Files_To_List, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []  # no need of third party software!

        # Inputs description
        file1_desc = "A string corresponding to an existing path file."
        file2_desc = (
            "An optional string corresponding to an existing path file."
        )
        file3_desc = (
            "An optional string corresponding to an existing path file."
        )

        # Outputs description
        file_list_desc = "A list of items which are path files."

        # Inputs traits
        self.add_trait("file1", traits.File(output=False, desc=file1_desc))
        self.file1 = traits.Undefined

        self.add_trait(
            "file2", traits.File(output=False, optional=True, desc=file2_desc)
        )
        self.file2 = traits.Undefined

        self.add_trait(
            "file3", traits.File(output=False, optional=True, desc=file3_desc)
        )
        self.file3 = traits.Undefined

        # Outputs traits
        self.add_trait(
            "file_list", traits.List(output=True, desc=file_list_desc)
        )
        self.file_list = traits.Undefined

        self.init_default_traits()

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
        if self.file1 not in ["<undefined>", traits.Undefined]:
            if self.file2 in [
                "<undefined>",
                traits.Undefined,
            ] and self.file3 in ["<undefined>", traits.Undefined]:
                self.outputs["file_list"] = [self.file1]

            elif self.file3 in [
                "<undefined>",
                traits.Undefined,
            ] and self.file2 not in ["<undefined>", traits.Undefined]:
                self.outputs["file_list"] = [self.file1, self.file2]

            elif self.file2 in [
                "<undefined>",
                traits.Undefined,
            ] and self.file3 not in ["<undefined>", traits.Undefined]:
                self.outputs["file_list"] = [self.file1, self.file3]

            else:
                self.outputs["file_list"] = [
                    self.file1,
                    self.file2,
                    self.file3,
                ]

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

    Please, see the complete documentation for the `Filter_Files_List brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/Filter_Files_List.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Filter_Files_List, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []  # no need of third party software!

        # Inputs description
        in_list_desc = "The list of elements to be filtered."
        index_filter_desc = "A list of 0 to 2 indexes for filtering."

        # Outputs description
        filtered_list_desc = "The corresponding filtering result (a list)."

        # Inputs traits
        self.add_trait(
            "in_list",
            traits.List(traits.File, output=False, desc=in_list_desc),
        )
        self.in_list = traits.Undefined

        self.add_trait(
            "index_filter",
            traits.List(
                value=[1],
                trait=traits.Range(low=1, high=None),
                minlen=0,
                maxlen=2,
                output=False,
                optional=True,
                desc=index_filter_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "filtered_list", traits.List(output=True, desc=filtered_list_desc)
        )
        self.filtered_list = traits.Undefined

        self.init_default_traits()

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
        if self.in_list and self.in_list not in [
            "<undefined>",
            traits.Undefined,
        ]:
            if not self.index_filter:
                self.outputs["filtered_list"] = [self.in_list[0]]

            if len(self.index_filter) == 1:
                if self.index_filter[0] <= len(self.in_list):
                    self.outputs["filtered_list"] = [
                        self.in_list[self.index_filter[0] - 1]
                    ]

                else:
                    print(
                        "\nFilter_Files_List brick: Initialisation failed "
                        "because the index_filter parameter is greater than "
                        "the length of the in_list parameter ...\n"
                    )

            if len(self.index_filter) == 2:
                if self.index_filter[0] < self.index_filter[1]:
                    if self.index_filter[0] <= len(self.in_list):
                        if self.index_filter[1] <= len(self.in_list):
                            # fmt: off
                            self.outputs["filtered_list"] = self.in_list[
                                self.index_filter[0] - 1:self.index_filter[1]
                            ]
                            # fmt: on
                        else:
                            print(
                                "\nFilter_Files_List brick: Initialisation "
                                "failed because the second value of the "
                                "index_filter parameter is greater than the "
                                "length of the in_list parameter ...\n"
                            )

                    else:
                        print(
                            "\nFilter_Files_List brick: Initialisation failed"
                            " because the first value of the index_filter "
                            "parameter is greater than the length of the "
                            "in_list parameter ...\n"
                        )

                else:
                    print(
                        "\nFilter_Files_List brick: Initialisation failed "
                        "because the first value of the index_filter parameter"
                        " is greater than the second ...\n"
                    )

            if self.outputs:
                self.outputs["notInDb"] = ["filtered_list"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class Find_In_List(ProcessMIA):
    """
    *From a list of files, select the 1rst element that contains a pattern*

    Please, see the complete documentation for the `Find_In_List brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/Find_In_List.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Find_In_List, self).__init__()

        # Inputs description
        in_list_desc = (
            "List of files to parse. A list of items which are "
            "an existing, uncompressed file "
            "(valid extensions: [.img, .nii, .hdr])."
        )
        pattern_desc = "Pattern to look for (a string)."

        # Outputs description
        out_file_desc = (
            "The first file found in the in_list parameter with "
            "the pattern in its name."
        )

        # Input traits
        self.add_trait(
            "in_list",
            InputMultiPath(
                traits.Either(ImageFileSPM(), traits.Undefined),
                copyfile=False,
                output=False,
                optional=False,
                desc=in_list_desc,
            ),
        )
        self.add_trait(
            "pattern",
            traits.String(
                "0001", output=False, optional=True, desc=pattern_desc
            ),
        )

        # Output traits
        self.add_trait(
            "out_file", ImageFileSPM(output=True, desc=out_file_desc)
        )

        self.init_default_traits()

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
        super(Find_In_List, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_list and self.in_list not in [
            "<undefined>",
            traits.Undefined,
        ]:
            for fil in self.in_list:
                if self.pattern in fil:
                    self.outputs["out_file"] = fil
                    break

            if self.outputs:
                self.outputs["notInDb"] = ["out_file"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class Get_Conditions_From_csv(ProcessMIA):
    """
    *Get conditions information (conditions names, onsets and durations)
    for Level1Design brick using csv files.*

    Please, see the complete documentation for
    the `Get_Conditions_From_csv brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/Get_Conditions_From_csv.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Get_Conditions_From_csv, self).__init__()

        # Inputs description
        csv_files_desc = (
            ".csv files contening the onset, one for each session "
            "(existing .csv files)"
        )
        design_type_desc = (
            "Type of design for each session (List of string among "
            "bloc or event-related)."
        )

        # Outputs description
        cond_names_desc = ""
        cond_onsets_desc = ""
        cond_durations_desc = ""

        # Input traits
        self.add_trait(
            "csv_files",
            InputMultiPath(
                traits.File(),
                output=False,
                optional=False,
                desc=csv_files_desc,
            ),
        )

        self.add_trait(
            "design_type",
            traits.List(
                traits.Enum("bloc", "event-related"),
                output=False,
                optional=True,
                desc=design_type_desc,
            ),
        )

        # Output traits
        self.add_trait(
            "cond_names",
            traits.List(
                traits.List(traits.String()),
                value=[[]],
                output=True,
                optional=True,
                desc=cond_names_desc,
            ),
        )

        self.add_trait(
            "cond_onsets",
            traits.List(
                traits.List(traits.List(traits.Float())),
                value=[[[]]],
                output=True,
                optional=True,
                desc=cond_onsets_desc,
            ),
        )

        self.add_trait(
            "cond_durations",
            traits.List(
                traits.List(traits.List(traits.Float())),
                value=[[[]]],
                output=True,
                optional=True,
                desc=cond_durations_desc,
            ),
        )

        self.init_default_traits()

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
        super(Get_Conditions_From_csv, self).list_outputs()

        if len(self.csv_files) != len(self.design_type):
            print(
                "\nGet_Conditions_From_csv brick: Initialization failed... "
                "Please precise design type for each csv file"
            )
            return self.make_initResult()

        all_cond_names = []
        all_cond_onsets = []
        all_cond_durations = []

        for i in range(len(self.csv_files)):
            cond_names = []
            cond_onsets = []
            cond_durations = []

            csv_file = self.csv_files[i]
            design = self.design_type[i]
            # Check extension
            valid_bool, in_ext, file_name = checkFileExt(
                csv_file, {"csv": "csv"}
            )
            if not valid_bool:
                print(
                    "\nGet_Conditions_From_csv brick: Initialization "
                    "failed... One of the file is not a .csv file ...!"
                )
                return self.make_initResult()

            # Get infos into csv
            df = pd.read_csv(csv_file)
            col_names = list(df.columns)
            if design == "bloc":
                cond_names = []
                for i in range(len(col_names)):
                    if "duration" not in col_names[i]:
                        # Check that we have duration for each condition
                        try:
                            df.loc[:, col_names[i] + " duration"]
                            cond_names.append(col_names[i])
                        except Exception:
                            print(
                                "\nGet_Conditions_From_csv brick: "
                                "Initialization failed... "
                                "For block design, duration should be "
                                "in the csv file ...!"
                            )
                            return self.make_initResult()
            elif design == "event-related":
                cond_names = col_names

            for name in cond_names:
                if design == "event-related":
                    cond_onsets.append(df.loc[:, name].tolist())
                    cond_durations.append([0])
                elif design == "bloc":
                    cond_onsets.append(df.loc[:, name].tolist())
                    cond_durations.append(
                        df.loc[:, name + " duration"].tolist()
                    )

            all_cond_names.append(cond_names)
            all_cond_onsets.append(cond_onsets)
            all_cond_durations.append(cond_durations)

        # Outputs definition and tags inheritance (optional)
        if all_cond_names:
            self.outputs["cond_names"] = all_cond_names
            self.outputs["cond_onsets"] = all_cond_onsets
            self.outputs["cond_durations"] = all_cond_durations

        if self.outputs:
            self.outputs["notInDb"] = [
                "cond_names",
                "cond_onsets",
                "cond_durations",
            ]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        # super(Get_Conditions_From_csv, self).run_process_mia()
        return


class Get_Patient_Name(ProcessMIA):
    """
    *Get patient name from a file*

    Please, see the complete documentation for the `Get_Patient_Name brick
    in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/Get_Patient_Name.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Get_Patient_Name, self).__init__()

        # Inputs description
        in_file_desc = "In file"

        # Outputs description
        patient_name_desc = "Patient name."

        # Input traits
        self.add_trait(
            "in_file",
            File(
                output=False,
                copyfile=True,
                optional=False,
                desc=in_file_desc,
            ),
        )

        # Output traits
        self.add_trait(
            "patient_name", traits.String(output=True, desc=patient_name_desc)
        )

        self.init_default_traits()

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
        super(Get_Patient_Name, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            sub_name = get_dbFieldValue(
                self.project, self.in_file, "PatientName"
            )
            if sub_name is None:
                print(
                    "(Get_Patient_Name brick: Please, fill 'PatientName' tag "
                    "in the database for in file"
                )
                return self.make_initResult()
            self.outputs["patient_name"] = sub_name

            if self.outputs:
                self.outputs["notInDb"] = ["patient_name"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class Import_Data(ProcessMIA):
    """
    *Import reference data into the current pipeline*

    - This brick was originally written to select regions of interest.
      The rois_list parameter is used to filter the data to be imported from
      a library defined by lib_dir. If roi_list is a list, each element of
      it will be a filename filter applied for retrieval. If roi_list is a
      list of lists, the filters will result from concatenating the elements
      of each internal list (e.g. [["foo", "1"], ["faa", "2"]] gives two
      filters, "foo_1" and "faa_2".
    - If lib_dir is not set, use the miaresources/ROIs/ file.
    - The file_in_db file is used only to retrieve the value of the
      associated PatientName tag.
    - The reference data is imported into the
      output_directory/PatientName_data/ROI_data/raw_data directory.
    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Import_Data, self).__init__()

        # Inputs description
        rois_list_desc = (
            "A list or a list of lists of strings, defining the "
            "data to import."
        )
        lib_dir_desc = (
            "The path to a data library (if not defined, it uses "
            "the default resources path, i.e. miaresources/ROIs/)."
        )
        file_in_db_desc = (
            "A file in database, only used to catch "
            "the PatientName tag value."
        )
        starts_with_desc = (
            "If True applies the file filter only to the "
            "beginning of the names, otherwise to the whole "
            "names (a boolean)."
        )
        # Outputs description
        rois_files_desc = "The list of resulting available files."

        # Inputs traits
        self.add_trait(
            "rois_list",
            traits.Either(
                traits.List(traits.String()),
                traits.List(
                    traits.List(
                        traits.String(),
                        minlen=2,
                        maxlen=2,
                    )
                ),
                output=False,
                optional=False,
                desc=rois_list_desc,
            ),
        ),

        self.add_trait(
            "lib_dir",
            traits.Directory(output=False, optional=True, desc=lib_dir_desc),
        )
        self.lib_dir = traits.Undefined

        self.add_trait(
            "file_in_db",
            File(output=False, optional=False, desc=file_in_db_desc),
        )

        self.add_trait(
            "starts_with",
            traits.Bool(
                default_value=True,
                output=False,
                optional=True,
                desc=starts_with_desc,
            ),
        )

        # Outputs traits
        self.add_trait(
            "rois_files",
            OutputMultiPath(File(), output=True, desc=rois_files_desc),
        )

        # Special parameter used as a messenger for the run_process_mia method
        self.add_trait(
            "dict4runtime",
            traits.Dict(output=False, optional=True, userlevel=1),
        )

        self.init_default_traits()

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
        super(Import_Data, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if (
            self.rois_list not in (traits.Undefined, [], [[]])
            and self.file_in_db != traits.Undefined
        ):
            patient_name = get_dbFieldValue(
                self.project, self.file_in_db, "PatientName"
            )

            if patient_name is None:
                print(
                    "\nImport_Data brick:\n The 'PatientName' tag is not "
                    "filled in the database for the {} file ...\n The "
                    "calculation is aborted...".format(self.file_in_db)
                )
                return self.make_initResult()

            self.dict4runtime["patient_name"] = patient_name
            roi_raw_data_dir = os.path.join(
                self.output_directory,
                patient_name + "_data",
                "ROI_data",
                "raw_data",
            )

            if self.lib_dir in (traits.Undefined, ""):
                config = Config()
                self.lib_dir = os.path.join(
                    config.get_resources_path(), "ROIs"
                )

            # list the content of self.lib_dir
            elts = os.listdir(self.lib_dir)
            # filtering only the files
            all_ref_files = [
                f
                for f in elts
                if os.path.isfile(os.path.join(self.lib_dir, f))
            ]

            if isinstance(self.rois_list[0], list):
                filters = [f[0] + "_" + f[1] for f in self.rois_list]

            else:
                filters = self.rois_list

            list_out = []

            for ref_file in all_ref_files:
                for _filter in filters:
                    if self.starts_with is True and ref_file.startswith(
                        _filter
                    ):
                        keep = True

                    elif self.starts_with is False and _filter in ref_file:
                        keep = True

                    else:
                        keep = False

                    if keep is True:
                        list_out.append(
                            os.path.join(roi_raw_data_dir, ref_file)
                        )
                        self.dict4runtime[
                            os.path.join(self.lib_dir, ref_file)
                        ] = os.path.join(roi_raw_data_dir, ref_file)

            self.outputs["rois_files"] = sorted(list_out)

            if self.outputs:
                self.outputs["notInDb"] = ["rois_files"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        # No need the next line (we don't use self.process and SPM)
        # super(Import_Data, self).run_process_mia()

        patient_name = self.dict4runtime["patient_name"]

        pat_name_dir = os.path.join(
            self.output_directory, patient_name + "_data"
        )

        if not os.path.exists(pat_name_dir):
            os.mkdir(pat_name_dir)

        roi_data_dir = os.path.join(pat_name_dir, "ROI_data")

        if not os.path.exists(roi_data_dir):
            os.mkdir(roi_data_dir)

        roi_raw_data_dir = os.path.join(roi_data_dir, "raw_data")

        tmp = "None"

        if os.path.isdir(roi_raw_data_dir):
            tmp = tempfile.mktemp(dir=os.path.dirname(roi_data_dir))
            os.mkdir(tmp)
            shutil.move(roi_raw_data_dir, os.path.join(tmp, "raw_data"))
            print(
                '\nImportData brick:\nA "{}" folder already exists, '
                "it will be overwritten by this new "
                "calculation...".format(roi_raw_data_dir)
            )
        os.mkdir(roi_raw_data_dir)

        if os.path.isdir(tmp):
            shutil.rmtree(tmp)

        # Copying the selected ROIs from the selected resources folder
        # to roi_raw_data_dir
        for k in self.dict4runtime:
            if k != "patient_name":
                shutil.copyfile(k, self.dict4runtime[k])


class Input_Filter(ProcessMIA):
    """
    *To filter the Data Browser tab or the output data of another brick*

    Please, see the complete documentation for the
    `Input_Filter in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/Input_Filter.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Input_Filter, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []  # no need of third party software!

        # Inputs description
        input_desc = (
            "A list corresponding to the data from the Data Browser "
            "or the output data from another brick."
        )

        # Outputs description
        output_desc = "A list with the result of the filter applied."

        # Inputs traits
        self.add_trait(
            "input", traits.List(traits.File, output=False, desc=input_desc)
        )
        self.input = traits.Undefined

        # Outputs traits
        self.add_trait(
            "output", traits.List(traits.File, output=True, desc=output_desc)
        )
        self.output = traits.Undefined

        # Instantiation of the filter object
        self.filter = Filter(
            None, [""], [""], [["FileName"]], [], ["CONTAINS"], ""
        )

        self.init_default_traits()

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
        #       open_filter of PipelineEditor?

        # Outputs definition and tags inheritance (optional)

        # Getting the input:
        # The current way for obtaining input data and filtering them is very
        # permissive. Indeed, it is possible to obtain the data from the output
        # of a previous brick or from the database (using Export to
        # database_scans, or not) and of course to perform a filtering from it.
        # But also, if no data is sent in, we get all the data from the
        # database (self.project.session.get_documents_names(
        # COLLECTION_CURRENT), below). It can maybe lead to side effects
        # (however it also allows for example not to connect the input to
        # something and to have the whole database by default...
        # Is it really desirable????
        if self.input:
            self.scans_list = self.input

        else:
            self.scans_list = self.project.session.get_documents_names(
                COLLECTION_CURRENT
            )

        # The data to filter are always a relative path
        if (self.scans_list) and (os.path.isabs(self.scans_list[0])):
            self.scans_list = [
                os.path.relpath(i, self.project.folder)
                for i in self.scans_list
            ]

        # Apply the filter to the input
        self.outputs["output"] = self.filter.generate_filter(
            self.project,
            self.scans_list,
            self.project.session.get_shown_tags(),
        )
        self.outputs["notInDb"] = ["output"]

        # The output data are always an absolute path
        for idx, element in enumerate(self.outputs["output"]):
            full_path = os.path.abspath(
                os.path.join(self.project.folder, element)
            )
            self.outputs["output"][idx] = full_path

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class List_Duplicate(ProcessMIA):
    """
    *From a file name, generate a list containing that file name and the
    file name itself*

    Please, see the complete documentation for the
    `List_Duplicate in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/List_Duplicate.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(List_Duplicate, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []  # no need of third party software!

        # Inputs description
        file_name_desc = "A string corresponding to an existing path file."

        # Outputs description
        out_file_desc = "A string corresponding to an existing path file."
        out_list_desc = (
            "A list with one string element corresponding "
            " to an existing path file."
        )

        # Inputs traits
        self.add_trait("file_name", File(output=False, desc=file_name_desc))

        # Outputs traits
        self.add_trait(
            "out_file", traits.File(output=True, desc=out_file_desc)
        )
        self.out_file = traits.Undefined

        self.add_trait(
            "out_list", traits.List(output=True, desc=out_list_desc)
        )
        self.out_list = traits.Undefined

        self.init_default_traits()

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
            self.outputs["out_file"] = self.file_name
            self.outputs["notInDb"] = ["out_list", "out_file"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class List_To_File(ProcessMIA):
    """
    *From several filenames, selects and generates a file*

    Please, see the complete documentation for the
    `List_To_File in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/List_To_File.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(List_To_File, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []  # no need of third party software!

        # Inputs description
        file_list_desc = "The list of elements to be filtered."
        index_filter_desc = "A list of 0 to 1 indexes for filtering."

        # Outputs description
        file_desc = "The corresponding filtering result (a file)."

        # Inputs traits
        self.add_trait(
            "file_list",
            traits.List(traits.File(), output=False, desc=file_list_desc),
        )
        self.file_list = traits.Undefined

        self.add_trait(
            "index_filter",
            traits.List(
                value=[1],
                trait=traits.Range(low=1, high=None),
                minlen=0,
                maxlen=1,
                output=False,
                optional=True,
                desc=index_filter_desc,
            ),
        )

        # Outputs traits
        self.add_trait("file", traits.File(output=True, desc=file_desc))
        self.file = traits.Undefined

        self.init_default_traits()

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
        if self.file_list and self.file_list not in [
            "<undefined>",
            traits.Undefined,
        ]:
            if not self.index_filter:
                self.outputs["file"] = [self.file_list[0]]

            if len(self.index_filter) == 1:
                if self.index_filter[0] <= len(self.file_list):
                    self.outputs["file"] = self.file_list[
                        self.index_filter[0] - 1
                    ]

                else:
                    print(
                        "\nList_To_File brick brick: Initialisation failed "
                        "because the index_filter parameter is greater than "
                        "the length of file_list parameter ...\n"
                    )

        if self.outputs:
            self.outputs["notInDb"] = ["file"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class Make_A_List(ProcessMIA):
    """
    *From 2 objects, generating a list containing all these objects*

    Please, see the complete documentation for the
    `Make_A_List in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/Make_A_List.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Make_A_List, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []  # no need of third party software!

        # Inputs description
        obj1_desc = "An object ..."
        obj2_desc = "An optional object ..."

        # Outputs description
        obj_list_desc = "A list of objects."

        # Inputs traits
        self.add_trait("obj1", traits.Any(output=False, desc=obj1_desc))
        self.obj1 = traits.Undefined

        self.add_trait(
            "obj2", traits.Any(output=False, optional=True, desc=obj2_desc)
        )
        self.obj2 = traits.Undefined

        # Outputs traits
        self.add_trait(
            "obj_list", traits.List(output=True, desc=obj_list_desc)
        )
        self.obj_list = traits.Undefined

        self.init_default_traits()

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
        super(Make_A_List, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.outputs:
            self.outputs = {}

        if self.obj1 and self.obj1 not in ["<undefined>", traits.Undefined]:

            try:
                objt1 = literal_eval(self.obj1)

            except (SyntaxError, ValueError):
                #  TODO: is it enough to do just that?
                objt1 = self.obj1

            if (not self.obj2) or (
                self.obj2 in ["<undefined>", traits.Undefined]
            ):
                self.outputs["obj_list"] = [objt1]

            else:

                try:
                    objt2 = literal_eval(self.obj2)

                except (SyntaxError, ValueError):
                    #  TODO: is it enough to do just that?
                    objt2 = self.obj2

                self.outputs["obj_list"] = [objt1, objt2]

            self.outputs["notInDb"] = ["obj_list"]

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        return


class Make_CVR_reg_physio(ProcessMIA):
    """
    *Generate the physiological regressor for CVR*

    Please, see the complete documentation for the
    `Make_CVR_reg_physio in the mia_processes website
    <https://populse.github.io/mia_processes/html/documentation/bricks/tools/Make_CVR_reg_physio.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Make_CVR_reg_physio, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []  # no need of third party software!

        # Inputs description
        func_file_desc = "The fMRI scan under hypercapnia (a file)"
        trigger_data_desc = (
            "The trigger data (a file with extension "
            "in ['.txt', '.csv', '.log'])"
        )
        physio_data_desc = (
            "The physiological data (a file with extension "
            "in ['.txt', '.csv'])"
        )

        # Outputs description
        cvr_reg_desc = "The physiological regressor for CVR (a file)"

        # Inputs traits
        self.add_trait(
            "func_file",
            Either(
                ImageFileSPM(),
                Undefined,
                copyfile=False,
                output=False,
                optional=False,
                desc=func_file_desc,
            ),
        )
        self.func_file = traits.Undefined

        self.add_trait(
            "trigger_data",
            File(output=False, optional=True, desc=trigger_data_desc),
        )
        self.trigger_data = traits.Undefined

        self.add_trait(
            "physio_data",
            File(output=False, optional=True, desc=physio_data_desc),
        )
        self.physio_data = traits.Undefined

        # Outputs traits
        self.add_trait("cvr_reg", File(output=True, desc=cvr_reg_desc))
        self.cvr_reg = traits.Undefined

        self.init_default_traits()

    def gaussfir(self, bt, nt=3, of=2):
        """Generate a Gaussian FIR filter

        with a specified bandwidth-time product (bt),
        number of taps (nt), and oversampling factor (of)
        """

        def convert_to_t(of, nt):
            """Convert to t in which to compute the filter coefficients."""
            # Filter Length
            filt_len = 2 * of * nt + 1
            t = np.linspace(-nt, nt, filt_len)
            return t

        t = convert_to_t(of, nt)
        alpha = np.sqrt(np.log(2) / 2) / (bt)
        h = (np.sqrt(np.pi) / alpha) * np.exp(-((t * np.pi / alpha) ** 2))
        # Normalize coefficients
        h = h / np.sum(h)
        return h

    def gfb_conv(self, a, b, shape="full"):
        """Return a subsection of the convolution, as specified by the
        shape parameter.

        :param a: a vector (a ndarray numpy object)
        :param b: a vector (a ndarray numpy object)
        :param shape: sub-section of the convolution (a string in
                      ['full', 'same', 'valid'])
        """
        if shape not in ["full", "same", "valid"]:
            raise ValueError(
                "Make_CVR_reg_physio brick: Invalid value for 'shape'. "
                "It must be 'full', 'same', or 'valid'."
            )

        na = np.prod(a.shape)
        nb = np.prod(b.shape)
        # Check if a is a column vector
        iscolumn = np.size(np.shape(a)) == 1
        # Reshape a and b into column vectors
        a = a.flatten()
        b = b.flatten()

        # Determine the size of the output convolution array
        if shape == "full":
            k_max = max(na + nb - 1, 0)
            ia = np.arange(1, k_max + 1)

        elif shape == "same":
            k_max = na
            ia = np.arange(1, k_max + 1) + np.floor(nb / 2)

        elif shape == "valid":
            k_max = max(na - max(nb - 1, 0), 0)
            ia = np.arange(1, k_max + 1) + max(nb - 1, 0)

        ia = np.repeat(ia.reshape(-1, 1), nb, axis=1)
        j = np.repeat(np.arange(1, nb + 1).reshape(1, -1), k_max, axis=0)
        ia = ia + 1 - j
        valid = np.logical_and(ia > 0, ia <= na)
        ia[~valid] = na + 1
        a = np.append(a, 0)  # Add zero to end of a
        c = np.sum(a[ia - 1] * b, axis=1)

        return c if iscolumn else c.reshape(1, -1)

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
        super(Make_CVR_reg_physio, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        # if self.outputs:
        #     self.outputs = {}

        if self.func_file not in ["<undefined>", traits.Undefined]:
            patient_name = get_dbFieldValue(
                self.project, self.func_file, "PatientName"
            )

        else:
            print(
                "Make_CVR_reg_physio brick: The mandatory 'func_file' input "
                "parameter is not defined ..."
            )
            return self.make_initResult()

        if patient_name is None:
            print(
                "Make_CVR_reg_physio brick: Please, fill 'PatientName' "
                "tag in the database for {} ...".format(self.func_file)
            )
            return self.make_initResult()

        if self.output_directory:
            # Create a data directory for this patient
            out_directory = os.path.join(
                self.output_directory, patient_name + "_data"
            )

            if not os.path.exists(out_directory):
                os.mkdir(out_directory)

        else:
            print(
                "Make_CVR_reg_physio brick: No output_directory was "
                "found...!\n"
            )
            return self.make_initResult()

        fname_reg = os.path.join(out_directory, "CVR_physio_reg.mat")
        all_tags_to_add = []
        # Add the reg state (individual or standard), in database
        tag_to_add = dict()
        tag_to_add["name"] = "Regressor state"
        tag_to_add["field_type"] = "string"
        tag_to_add["description"] = (
            "Indicates the state of the physiological regressor for CVR "
            "(Individual or Standard)"
        )
        tag_to_add["visibility"] = True
        tag_to_add["origin"] = "user"
        tag_to_add["unit"] = None
        tag_to_add["default_value"] = None

        try:
            tag_to_add["value"] = "Individual"
            # TODO: The following attributes are currently hard-coded, but
            #       we could add them as input of the brick if necessary?
            # the number of last dyns to be deleted
            end_dyn_sup = 0
            # delay between respiration and EtCO2 variation (s)
            delay_for_etco2 = 4.8
            tr = (
                get_dbFieldValue(
                    self.project, self.func_file, "RepetitionTime"
                )[0]
                / 1000
            )
            # TODO: glover_hrf(tr) don't give good result yet. While waiting
            #       to make a good hrf, we use the hrf calculated as in spm
            # from nilearn.glm.first_level.hemodynamic_models import glover_hrf
            # hrf = glover_hrf(tr)
            hrf = self.spm_hrf(tr, [6, 16, 2, 1, np.inf, 0, 32])
            nb_dyn = (
                get_dbFieldValue(
                    self.project,
                    self.func_file,
                    "Dataset dimensions (Count, X,Y,Z,T...)",
                )[4]
                - end_dyn_sup
            )
            # ---- Start of making phys_trig_data ----
            # record some data prior to start and after end of scan,
            # typically the duration of one hrf
            time_margin = 30  # seconds
            status_sample_freq = 37.127
            # ---- Start of making trigdata ----
            starttime = np.nan
            endtime = np.inf
            n_pulses = 0
            # Read data from different types of files
            file_extension = self.trigger_data.lower().split(".")[-1]

            if file_extension == "txt":
                # TODO: Not yet tested

                with open(self.trigger_data, "r") as fil:
                    lines = fil.readlines()

                    for line in lines:

                        if np.isnan(starttime):
                            colonidxs = [
                                i for i, c in enumerate(line) if c == ":"
                            ]

                            if len(colonidxs) == 2:
                                # fmt: off
                                this_time_s = (
                                    3600
                                    * int(
                                        line[colonidxs[0] - 2:colonidxs[0]]
                                    )
                                    + 60
                                    * int(
                                        line[
                                            colonidxs[0] + 1:colonidxs[1] - 1
                                        ]
                                    )
                                    + int(
                                        line[
                                            colonidxs[1] + 1:colonidxs[1] + 2
                                        ]
                                    )
                                )
                                # fmt: on
                                starttime = this_time_s

                        if line.startswith("E11"):
                            this_data = line.split(",")
                            dt_ms = float(this_data[5])
                            dt_s = dt_ms / 1000
                            next_line = lines[lines.index(line) + 1]
                            trig_times_ms = [
                                int(x) for x in next_line.split(",")
                            ]
                            trig_times_s = np.array(trig_times_ms) / 1000
                            endtime = starttime + len(trig_times_s) * dt_s
                            break

                tr = np.median(np.diff(trig_times_s))

            elif file_extension == "log":
                # tested: OK
                with open(self.trigger_data, "r") as fil:
                    lines = fil.readlines()[5:]  # Skip header lines
                    trig_times = []

                    for line in lines:

                        if np.isnan(starttime):
                            colonidxs = [
                                i for i, c in enumerate(line) if c == ":"
                            ]

                            if len(colonidxs) == 2:
                                # fmt: off
                                this_time_s = (
                                    3600
                                    * int(
                                        line[colonidxs[0] - 2:colonidxs[0]]
                                    )
                                    + 60
                                    * int(
                                        line[colonidxs[0] + 1:colonidxs[1]]
                                    )
                                    + float(
                                        line[
                                            colonidxs[1] + 1:colonidxs[1] + 5
                                        ]
                                    )
                                )
                                # fmt: on
                                starttime_abs = this_time_s
                                tab_pos = line.split("\t")
                                starttime_rel = float(tab_pos[3]) / 10000
                                starttime = starttime_abs - starttime_rel

                        if "Pulse" in line:
                            tab_pos = line.split("\t")
                            trig_times.append(float(tab_pos[3]))

                trig_times_s = np.array(trig_times) / 10000
                starttime = starttime + trig_times_s[0]
                dt_s = np.diff(trig_times_s)
                tr = np.median(dt_s)
                endtime = starttime + trig_times_s[-1] - trig_times_s[0] + tr
                trig_times_s = trig_times_s - trig_times_s[0]
                n_pulses = int((endtime - starttime) / tr)
                # look for missing triggers
                dn = np.round(np.diff(trig_times_s) / tr) - 1
                dn_pos = np.where(dn != 0)[0]

                if dn_pos.size != 0:

                    for pos in dn_pos[::-1]:
                        # fmt: off
                        trig_times_s = np.concatenate(
                            [
                                trig_times_s[:pos + 1],
                                [tr + trig_times_s[pos]],
                                trig_times_s[pos + 1:],
                            ]
                        )
                        # fmt: on

            elif file_extension == "csv":
                # TODO: Not yet tested
                # Extract start time from file name
                short_fname = os.path.splitext(
                    os.path.basename(self.trigger_data)
                )[0]
                starttime = int(short_fname.split("_")[1])

                with open(self.trigger_data, "r") as fil:
                    lines = fil.readlines()[2:]  # Skip header lines
                    trig_times = []

                    for line in lines:
                        data = line.split(";")
                        trig_times.append(float(data[0]))

                trig_times_s = np.array(trig_times)
                starttime = starttime + trig_times_s[0]
                tr = np.median(np.diff(trig_times_s))
                endtime = starttime + trig_times_s[-1] + tr
                trig_times_s = trig_times_s - trig_times_s[0]
                n_pulses = len(trig_times_s)

            else:
                raise FileNotFoundError(
                    "Make_CVR_reg_physio brick: The Only "
                    ".csv, .log and .txt extensions are "
                    "allowed for the trigger_data "
                    "parameter. The individual regressor "
                    "cannot be generated... "
                )
            # ---- End of making trigdata ----
            starttime = float(starttime)
            endtime = float(endtime)

            if self.physio_data.lower().endswith(".csv"):
                # tested: OK
                print("Data from Magdata software detected ...")
                times = []
                data_s = []
                paramnames = None

                with open(self.physio_data, "r") as fid:

                    for line in fid:
                        times.append(line[:9])

                        if line[9:].endswith("\n"):
                            data_s.append(line[9:-2])

                        else:
                            data_s.append(line[9:])

                nlines = len(times)
                phys_trig_data = {}
                time_matrix = np.zeros(nlines)
                datapoint = 0
                have_format = False

                for line in range(nlines):
                    this_times = times[line]

                    try:
                        # convert to seconds since midnight
                        this_time_s = sum(
                            [
                                time * factor
                                for time, factor in zip(
                                    map(int, this_times.split(":")),
                                    [3600, 60, 1],
                                )
                            ]
                        )

                    except ValueError:

                        if have_format is False:
                            formatline = this_times + data_s[line]
                            paramnames = formatline.split(",")
                            # We're not interested in the 'Time' field
                            # at the beginning
                            paramnames = paramnames[1:]
                            # Replace spaces with underscores and
                            # '%' with 'perc'
                            paramnames = [
                                name.strip()
                                .replace(" ", "_")
                                .replace("%", "perc")
                                for name in paramnames
                            ]

                            # If 'IBP2_wm2' is not present,
                            # insert the missing field right after 'IBP2_wm1'
                            if "IBP2_wm2" not in paramnames:
                                position = paramnames.index("IBP2_wm1")
                                paramnames.insert(position + 1, "IBP2_wm2")

                            # Fix typos in magdata-provided field names
                            paramnames = [
                                name.replace("Satus", "Status")
                                for name in paramnames
                            ]
                            n_params = len(paramnames)
                            # Mark data as unread
                            data_matrix = np.full((nlines, n_params), np.nan)
                            have_format = True

                    else:

                        if starttime == 0:
                            starttime = this_time_s
                            this_time_s = 0

                        else:

                            if this_time_s < starttime - time_margin:
                                continue  # discard this line

                            elif this_time_s > endtime + time_margin:
                                # don't expect to find any interesting
                                # data after this point
                                break

                            this_time_s = this_time_s - starttime

                        this_data = [
                            int(val) if val.strip() else np.nan
                            for val in data_s[line].split(",")
                        ]
                        this_nvalues = len(this_data) - 1
                        datapoint += 1
                        data_matrix[datapoint - 1, :this_nvalues] = this_data[
                            1:
                        ]
                        time_matrix[datapoint - 1] = this_time_s

                n_datapoints = datapoint
                data_matrix = data_matrix[:n_datapoints, :]
                time_matrix = time_matrix[:n_datapoints]

                # Correcting time data
                status_idx = np.where(~np.isnan(data_matrix[:, -1]))[0]
                time_stat = time_matrix[status_idx]
                n_timepoints = len(time_stat)

                # Smoothing time data
                n_interp = min(round(25 * status_sample_freq), n_timepoints)
                coeff_beg = np.linalg.lstsq(
                    np.vstack(
                        [(np.arange(n_interp) + 1), np.ones(n_interp)]
                    ).T,
                    time_stat[:n_interp],
                    rcond=None,
                )[0]
                coeff_end = np.linalg.lstsq(
                    np.vstack(
                        [
                            (
                                np.arange(
                                    n_timepoints - n_interp, n_timepoints
                                )
                                + 1
                            ),
                            np.ones(n_interp),
                        ]
                    ).T,
                    time_stat[-n_interp:],
                    rcond=None,
                )[0]
                h = self.gaussfir(0.3, 2, 2 * round(status_sample_freq))
                n_extrap_beg = int(np.ceil((len(h) - 1) / 2))
                n_extrap_end = len(h) - 1 - n_extrap_beg
                time_extrap_beg = np.dot(
                    np.hstack(
                        (
                            np.arange(-n_extrap_beg + 1, 1).reshape(-1, 1),
                            np.ones((n_extrap_beg, 1)),
                        )
                    ),
                    coeff_beg,
                )
                time_extrap_end = np.dot(
                    np.hstack(
                        (
                            np.arange(
                                n_timepoints + 1,
                                n_timepoints + n_extrap_end + 1,
                            ).reshape(-1, 1),
                            np.ones((n_extrap_end, 1)),
                        )
                    ),
                    coeff_end,
                )
                time_extrap = np.concatenate(
                    (time_extrap_beg, time_stat, time_extrap_end)
                )
                time_stat = self.gfb_conv(time_extrap, h, "valid")
                # We can also use np.convolve()!
                # time_stat = np.convolve(
                #     time_extrap.astype(int), h, mode="valid"
                # )

                # Check for points where time stands still
                time_stat_d = np.concatenate(([0], np.diff(time_stat)))
                # time_stat_d = np.diff(time_stat, prepend=0)
                stall_points = time_stat_d < 0.75 * (1 / status_sample_freq)
                stall_end = (
                    np.where(stall_points)[0][-1]
                    if np.any(stall_points)
                    else 0
                )
                regress_valid = np.vstack(
                    (
                        (np.arange(stall_end + 1, n_timepoints + 1)),
                        np.ones((n_timepoints - stall_end,)),
                    )
                )
                time_stat[stall_end:] = regress_valid.T.dot(
                    np.linalg.lstsq(
                        regress_valid.T, time_stat[stall_end:], rcond=None
                    )[0]
                )

                while stall_end > 0:
                    stall_start = (
                        np.where(
                            time_stat_d[:stall_end] > 1 / status_sample_freq
                        )[0][-1]
                        if np.any(
                            time_stat_d[:stall_end] > 1 / status_sample_freq
                        )
                        else 0
                    )
                    sample_interv_ok = (
                        time_stat_d[: stall_start + 1]
                        <= 1 / status_sample_freq
                    )
                    avg_sample_interv_ok = (
                        time_stat[stall_end] - time_stat[: stall_start + 1]
                    ) / (
                        stall_end - np.arange(stall_start + 1)
                    ) >= 1 / status_sample_freq
                    drift_start = (
                        np.where(sample_interv_ok & avg_sample_interv_ok)[0][
                            -1
                        ]
                        if np.any(sample_interv_ok & avg_sample_interv_ok)
                        else 0
                    )
                    old_stall_end = stall_end
                    time_stat_d = np.diff(time_stat[: drift_start + 1])
                    stall_points = time_stat_d < 0.75 * (
                        1 / status_sample_freq
                    )
                    stall_end = (
                        np.where(stall_points)[0][-1]
                        if np.any(stall_points)
                        else 0
                    )
                    regress_valid = np.vstack(
                        [
                            (np.arange(stall_end, drift_start + 1) + 1),
                            np.ones(drift_start - stall_end + 1),
                        ]
                    ).T
                    # fmt: off
                    time_stat[stall_end:drift_start + 1] = np.dot(
                        regress_valid,
                        np.linalg.lstsq(
                            regress_valid,
                            time_stat[stall_end:drift_start + 1],
                            rcond=None,
                        )[0],
                    )
                    time_stat[drift_start:old_stall_end + 1] = np.linspace(
                        time_stat[drift_start],
                        time_stat[old_stall_end],
                        old_stall_end - drift_start + 1,
                    )
                    # fmt: on

                if paramnames is None:
                    raise ValueError(
                        "Make_CVR_reg_physio brick: the "
                        "names of the physiological "
                        "parameters could not be found "
                        "in the physiological parameters "
                        "file!\nThe individual regressor "
                        "cannot be generated... "
                    )

                # Transfer data to the phys_trig_data dictionary
                for field in range(len(paramnames)):
                    this_fieldname = paramnames[field]
                    this_field_datapoints = ~np.isnan(data_matrix[:, field])
                    this_field_datapoints[: min(status_idx)] = False
                    # fmt: off
                    this_field_datapoints[max(status_idx) + 1:] = False
                    # fmt: on

                    if np.any(this_field_datapoints):
                        indx = np.where(this_field_datapoints)[0]
                        time = np.interp(indx, status_idx, time_stat)
                        data = data_matrix[this_field_datapoints, field]
                        phys_trig_data[this_fieldname] = {
                            "indx": indx,
                            "time": time,
                            "data": data,
                        }

                # Correct other bugs in magdata software
                problem_timepoints = (
                    phys_trig_data["Waveform_Status3"]["data"] != 0
                )

                if (
                    phys_trig_data["Pleth_wm2"]["data"].shape[0]
                    == n_timepoints
                ):
                    problem_timepoints = np.logical_or(
                        problem_timepoints,
                        phys_trig_data["Pleth_wm2"]["data"] == 103,
                    )

                if np.any(problem_timepoints):
                    valid_timepoints = np.where(~problem_timepoints)[0]

                    if (
                        phys_trig_data["Resp_wm"]["data"].shape[0]
                        == phys_trig_data["Capno_wm1"]["data"].shape[0]
                    ):
                        phys_trig_data["Capno_wm1"]["data"][
                            problem_timepoints
                        ] = phys_trig_data["Resp_wm"]["data"][
                            problem_timepoints
                        ]
                        problem_problem_timepoints = np.logical_and(
                            problem_timepoints,
                            phys_trig_data["Resp_wm"]["data"] == 0,
                        )
                        phys_trig_data["Capno_wm1"]["data"][
                            problem_problem_timepoints
                        ] = phys_trig_data["Waveform_Status1"]["data"][
                            problem_problem_timepoints
                        ]
                        phys_trig_data["Resp_wm"]["data"][
                            problem_timepoints
                        ] = phys_trig_data["Pleth_wm2"]["data"][
                            problem_timepoints
                        ]
                        phys_trig_data["Pleth_wm2"]["data"][
                            problem_timepoints
                        ] = phys_trig_data["Pleth_wm1"]["data"][
                            problem_timepoints
                        ]
                        problem_timepoints = np.where(problem_timepoints)[0]
                        phys_trig_data["Waveform_Status1"]["data"][
                            problem_timepoints
                        ] = np.maximum(
                            phys_trig_data["Waveform_Status1"]["data"][
                                np.maximum(problem_timepoints - 1, 0)
                            ],
                            phys_trig_data["Waveform_Status1"]["data"][
                                np.minimum(
                                    problem_timepoints + 1,
                                    phys_trig_data["Waveform_Status1"][
                                        "data"
                                    ].shape[0],
                                )
                            ],
                        )
                        phys_trig_data["Waveform_Status4"]["data"][
                            problem_timepoints
                        ] = 128 * (
                            (
                                phys_trig_data["Waveform_Status4"]["data"][
                                    np.maximum(problem_timepoints - 1, 0)
                                ]
                                > 64
                            ).any()
                            or (
                                phys_trig_data["Waveform_Status4"]["data"][
                                    np.minimum(
                                        problem_timepoints + 1,
                                        phys_trig_data["Waveform_Status4"][
                                            "data"
                                        ].shape[0]
                                        - 1,
                                    )
                                ]
                                > 64
                            ).any()
                        )

                    else:
                        phys_trig_data["Capno_wm1"]["data"][
                            problem_timepoints
                        ] = interp1d(
                            phys_trig_data["Capno_wm1"]["time"][
                                valid_timepoints
                            ],
                            phys_trig_data["Capno_wm1"]["data"][
                                valid_timepoints
                            ],
                            fill_value="extrapolate",
                        )(
                            phys_trig_data["Capno_wm1"]["time"][
                                problem_timepoints
                            ]
                        )

                        if (
                            data_matrix["Pleth_wm2"]["data"].shape[0]
                            == n_timepoints
                        ):
                            phys_trig_data["Pleth_wm2"]["data"][
                                problem_timepoints
                            ] = interp1d(
                                phys_trig_data["Pleth_wm2"]["time"][
                                    valid_timepoints
                                ],
                                phys_trig_data["Pleth_wm2"]["data"][
                                    valid_timepoints
                                ],
                                fill_value="extrapolate",
                            )(
                                phys_trig_data["Pleth_wm2"]["time"][
                                    problem_timepoints
                                ]
                            )
                        problem_timepoints = np.where(problem_timepoints)[0]
                        phys_trig_data["Waveform_Status1"]["data"][
                            problem_timepoints
                        ] = np.maximum(
                            phys_trig_data["Waveform_Status1"]["data"][
                                np.maximum(problem_timepoints - 1, 0)
                            ],
                            phys_trig_data["Waveform_Status1"]["data"][
                                np.minimum(
                                    problem_timepoints + 1,
                                    phys_trig_data["Waveform_Status1"][
                                        "data"
                                    ].shape[0],
                                )
                            ],
                        )
                        phys_trig_data["Waveform_Status4"]["data"][
                            problem_timepoints
                        ] = 128 * (
                            (
                                phys_trig_data["Waveform_Status4"]["data"][
                                    np.maximum(problem_timepoints - 1, 0)
                                ]
                                > 64
                            ).any()
                            or (
                                phys_trig_data["Waveform_Status4"]["data"][
                                    np.minimum(
                                        problem_timepoints + 1,
                                        phys_trig_data["Waveform_Status4"][
                                            "data"
                                        ].shape[0]
                                        - 1,
                                    )
                                ]
                                > 64
                            ).any()
                        )

                if n_pulses > 0:
                    phys_trig_data["trigger"] = {
                        "time": trig_times_s,
                        "data": np.ones_like(trig_times_s),
                    }

                phys_trig_data["starttime"] = starttime

            # Read CoolTerm capture log
            elif self.physio_data.lower().endswith(".txt"):
                # TODO: Not yet tested
                print("Data from CoolTerm application detected ...")
                # Process data and construct phys_trig_data dictionary
                # ...
                with open(self.physio_data, "r") as fid:
                    lines = fid.readlines()

                data = []

                for line in lines:
                    line_data = line.strip().split(",")
                    data.append(line_data)

                (
                    comp_date,
                    comp_time,
                    monitor_date,
                    monitor_time,
                    HR,
                    SPO2t,
                    _,
                    _,
                    _,
                    _,
                    _,
                    _,
                    _,
                    _,
                    _,
                    SpO2,
                    BR,
                    ETCO2,
                    ICO2,
                    _,
                    _,
                    _,
                    _,
                    _,
                    _,
                ) = zip(*data)

                # Convert sample start and end times to seconds since midnight
                monitortime = np.array(
                    [list(map(int, t.split(":"))) for t in monitor_time]
                )
                monitortime = (
                    3600 * monitortime[:, 0]
                    + 60 * monitortime[:, 1]
                    + monitortime[:, 2]
                )
                sampletime = np.array(
                    [list(map(int, t.split(":"))) for t in comp_time]
                )
                sampletime = (
                    3600 * sampletime[:, 0]
                    + 60 * sampletime[:, 1]
                    + sampletime[:, 2]
                )
                keep_idx = np.where(
                    np.logical_and(
                        sampletime >= starttime - time_margin,
                        sampletime <= endtime + time_margin,
                    )
                )[0]
                sampletime = sampletime[keep_idx] - starttime
                HR = np.array(HR)[keep_idx].astype(float)
                BR = np.array(BR)[keep_idx].astype(float)
                ICO2 = np.array(ICO2)[keep_idx].astype(float)
                ETCO2 = np.array(ETCO2)[keep_idx].astype(float)
                SpO2 = np.array(SpO2)[keep_idx].astype(float)

                # Transfer data to structure
                phys_trig_data = {
                    "HR": {
                        "indx": np.where(~np.isnan(HR))[0],
                        "time": sampletime[~np.isnan(HR)],
                        "data": HR[~np.isnan(HR)],
                    },
                    "Resp_Rate": {
                        "indx": np.where(~np.isnan(BR))[0],
                        "time": sampletime[~np.isnan(BR)],
                        "data": BR[~np.isnan(BR)],
                    },
                    "ICO2": {
                        "indx": np.where(~np.isnan(ICO2))[0],
                        "time": sampletime[~np.isnan(ICO2)],
                        "data": ICO2[~np.isnan(ICO2)],
                    },
                    "ETCO2": {
                        "indx": np.where(~np.isnan(ETCO2))[0],
                        "time": sampletime[~np.isnan(ETCO2)],
                        "data": ETCO2[~np.isnan(ETCO2)],
                    },
                    "SpO2_perc": {
                        "indx": np.where(~np.isnan(SpO2))[0],
                        "time": sampletime[~np.isnan(SpO2)],
                        "data": SpO2[~np.isnan(SpO2)],
                    },
                }

                if n_pulses > 0:  # We have read a trigger file
                    phys_trig_data["trigger"] = {
                        "time": trig_times_s,
                        "data": np.ones_like(trig_times_s),
                    }

                # In units of seconds since midnight
                phys_trig_data["starttime"] = starttime

            else:
                raise FileNotFoundError(
                    "Make_CVR_reg_physio brick: The Only "
                    ".csv and .txt extensions are "
                    "allowed for the physio_data "
                    "parameter. The individual regressor "
                    "cannot be generated... "
                )
            # ---- End of making phys_trig_data ----

            # performed for all odd frames
            triggers = phys_trig_data["trigger"]["time"]
            trig_TR = np.mean(np.diff(triggers))
            etco2data = phys_trig_data["ETCO2"]["data"]
            etco2time = phys_trig_data["ETCO2"]["time"] - delay_for_etco2
            outliers_etco2 = np.logical_or(etco2data > 100, etco2data < 10)
            valid_etco2 = ~outliers_etco2

            if not np.all(np.diff(etco2time[valid_etco2]) > 0):
                print(
                    "\nMake_CVR_reg_physio brick: The time relative to the "
                    "EtCO2 recording does not seem to be strictly monotonic "
                    "and increasing !!!\n"
                    "This makes no physical sense and the calculation of "
                    "the hypercapnic individual regressor will certainly "
                    "fail.\n"
                )

                print(
                    "\nMake_CVR_reg_physio brick: Trying to fix the "
                    "issue!!!\nPlease check the result carefully, this "
                    "is an automatic process ....\n\n"
                )

                x = etco2time[valid_etco2]
                x2, index = np.unique(x, return_index=True)
                etco2data[outliers_etco2[index]] = interp1d(
                    etco2time[outliers_etco2], x2, kind="cubic"
                )(etco2data[index])

            else:
                etco2data[outliers_etco2] = interp1d(
                    etco2time[valid_etco2],
                    etco2data[valid_etco2],
                    kind="cubic",
                )(etco2time[outliers_etco2])

            for iter in range(1, 3):
                # clean the etco2 data:
                # second derivative to get local outliers
                etco2data_pp = np.diff(etco2data.astype(np.int16), n=2)
                # detect outliers (very crude)
                outliers = np.abs(etco2data_pp) > 15
                # a single outlier creates three points of high curvature.
                # retain central point only:
                # TODO: Here, we can also use self.gfb_conv(). Test it out!
                outliers = (
                    np.convolve(
                        outliers.astype(int), np.array([1, 1, 1]), mode="same"
                    )
                    == 3
                )
                # correct shift due to derivative above
                outliers = np.where(outliers)[0] + 1
                etco2data = np.delete(etco2data, outliers)
                etco2time = np.delete(etco2time, outliers)

            # sort data and combine data from identical timepoints
            etco2time, sort_idx = np.sort(etco2time), np.argsort(etco2time)
            etco2data = etco2data[sort_idx]
            step_points = np.where(np.diff(etco2time) == 0)[0]
            step_points = np.union1d(step_points, step_points + 1)
            step_times = np.unique(etco2time[step_points])

            for k in range(len(step_times) - 1, -1, -1):
                this_points = np.where(etco2time == step_times[k])[0]
                etco2data[this_points[0]] = np.mean(etco2data[this_points])
                etco2data = np.delete(etco2data, this_points[1:])
                etco2time = np.delete(etco2time, this_points[1:])

            # add artificial triggers before and after scan (while we have
            # data) to help with convolution of regressors with hrf
            num_init_trigs = int(
                np.floor((min(triggers) - min(etco2time)) / trig_TR)
            )
            num_post_trigs = int(
                np.floor((max(etco2time) - max(triggers)) / trig_TR)
            )
            triggers = np.concatenate(
                (
                    trig_TR * (-np.arange(num_init_trigs, 0, -1))
                    + min(triggers),
                    triggers,
                    max(triggers) + trig_TR * np.arange(1, num_post_trigs + 1),
                )
            )
            etco2data = interp1d(etco2time, etco2data, kind="cubic")(triggers)
            # value separating 'hypercapnia' from 'normocapnia'
            median_etco2 = np.median(etco2data)
            # average of 'normocapnia'
            baseline_etco2 = np.mean(etco2data[etco2data < median_etco2])
            # shift baseline hc to zero, approximately
            etco2data_shift = etco2data - baseline_etco2
            etco2data_bold = signal.convolve(etco2data_shift, hrf, mode="full")
            # fmt: off
            io.savemat(
                fname_reg,
                {
                    "R": etco2data_bold[
                        num_init_trigs:num_init_trigs + nb_dyn
                    ].reshape((-1, 1))
                },
            )
            # fmt: on
            print(
                "\nMake_CVR_reg_physio brick: "
                "individual regressor generated!"
            )

        except Exception as e:
            print(
                "\nMake_CVR_reg_physio brick: The generation of the "
                "individual regressor has failed!\nTraceback:"
            )
            print("".join(traceback.format_tb(e.__traceback__)), end="")
            print("{0}: {1}\n".format(e.__class__.__name__, e))
            print("\n...Using the standard regressor...\n")
            tag_to_add["value"] = "Standard"
            config = Config()
            shutil.copy(
                os.path.join(
                    config.get_resources_path(),
                    "reference_population_data",
                    "regressor_physio_EtCO2_standard.mat",
                ),
                fname_reg,
            )
            print(
                "\nMake_CVR_reg_physio brick: " "standard regressor generated!"
            )

        all_tags_to_add.append(tag_to_add)
        # TODO: Can we simply add new tags without inheritance (I'm not sure
        #       this possibility was considered when we coded the
        #       tags_inheritance function)? As in this case, it might be
        #       preferable not to inherit from the functional but simply to
        #       add new tags.
        self.tags_inheritance(
            self.func_file, fname_reg, own_tags=all_tags_to_add
        )
        self.outputs["cvr_reg"] = fname_reg
        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def spm_hrf(self, rt, p=[6, 16, 1, 1, 6, 0, 32], t=16):
        """
        Hemodynamic response function.

        Parameters:
        - rt: scan repeat time
        - p: parameters of the response function (two Gamma functions)
        (Default: [6, 16, 1, 1, 6, 0, 32])
        - t: microtime resolution (Default: 16)

        Returns:
        - hrf: hemodynamic response function
        - p: parameters of the response function

        Note: Adapted from SPM12 matlab code (spm_hrf.m)
        """

        def gam_dis_pdf(ga, sh, sc):
            """
            Probability Density Function (PDF) of Gamma distribution.

            Parameters:
            - ga: Gamma-variate (Gamma has range [0,Inf))
            - sh: Shape parameter (h > 0)
            - sc: Scale parameter (l > 0)

            Returns:
            - f: PDF of Gamma-distribution with shape & scale parameters

            Note: Adapted from SPM12 matlab code (spm_Gpdf.m).
                  Warning: very specific to the spm_hrf case,
                  no argument validity test is performed.
            """
            # Initialise result to zeros
            f = np.zeros(len(ga))
            # Compute
            f[1:] = np.exp(
                (sh - 1) * np.log(ga[1:])
                + sh * np.log(sc)
                - sc * ga[1:]
                - gammaln(sh)
            )

            return f

        # from scipy.stats import gamma
        # Modelled hemodynamic response function - {mixture of Gammas}
        dt = rt / t
        u = np.arange(0, np.ceil(p[6] / dt) + 1) - p[5] / dt
        # hrf = (gamma.pdf(u, p[0] / p[2], scale=dt / p[2]) -
        #        gamma.pdf(u, p[1] / p[3], scale=dt / p[3]) / p[4])
        # We can also use the self.gam_dis_pdf() function:
        hrf = (
            gam_dis_pdf(u, p[0] / p[2], dt / p[2])
            - gam_dis_pdf(u, p[1] / p[3], dt / p[3]) / p[4]
        )
        hrf = hrf[(np.arange(0, np.floor(p[6] // rt) + 1) * t).astype(int)]
        hrf /= np.sum(hrf)

        return hrf

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        # The regressor is created at initialisation time rather than at
        # runtime because some bricks that use this regressor need to consult
        # it at initialisation time (e.g. Level1Design)!
        return
