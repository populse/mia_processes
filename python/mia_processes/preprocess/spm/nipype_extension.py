# -*- coding: utf-8 -*- #

"""This module is just a workaround solution to make the nipype NewSegment
brick available more quickly.

The purpose of this module is to override the special __init__ method of
the NewSegment class in order to remove the SPMCommand().version method
that is too long to execute when creating a pipeline.

:Contains:
    :Class:
        - NewSegmentMia

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# nipype import
from nipype.interfaces.spm.preprocess import (NewSegment,
                                              NewSegmentInputSpec,
                                              NewSegmentOutputSpec)
from nipype.interfaces.spm.base import SPMCommand


class NewSegmentMia(NewSegment):
    """ Derived class to avoid to instance SPMCommand().version that takes
    around 6 seconds. 

    """


    input_spec = NewSegmentInputSpec
    output_spec = NewSegmentOutputSpec

    def __init__(self, **inputs):
        self._jobtype = 'spatial'
        self._jobname = 'preproc'
        SPMCommand.__init__(self, **inputs)
