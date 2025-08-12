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
    Deconv_from_aif,
    Delete_data,
    Files_To_List,
    Filter_Files_List,
    Find_In_List,
    Get_Conditions_From_BIDS_tsv,
    Get_Conditions_From_csv,
    Get_Eprime_Info_GE2REC,
    Get_Patient_Name,
    Get_Regressors_From_csv,
    Import_Data,
    Input_Filter,
    List_Duplicate,
    List_Of_List_To_List,
    List_To_File,
    Make_A_List,
    Make_AIF,
    Make_CVR_reg_physio,
)
