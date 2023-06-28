# -*- coding: utf-8 -*-
"""
General bricks needed to build pipelines or to operate other
bricks (input data manipulations, etc.).

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from .tools import (  # noqa: F401
    Concat_to_list,
    Concat_to_list_of_list,
    Delete_data,
    Files_To_List,
    Filter_Files_List,
    Find_In_List,
    Get_Conditions_From_csv,
    Import_Data,
    Input_Filter,
    List_Duplicate,
    List_To_File,
    Make_A_List,
)
