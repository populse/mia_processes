# -*- coding: utf-8 -*- #

"""tests sandbox"""

# nipype import
from nipype.interfaces.base import traits

# populse_mia imports
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

from nipype.interfaces.afni import GCOR

import os
import os.path as op
import re
import numpy as np

from nipype.utils.filemanip import load_json, save_json, split_filename
from nipype.interfaces.base import (
    CommandLineInputSpec,
    CommandLine,
    Directory,
    TraitedSpec,
    traits,
    isdefined,
    File,
    InputMultiObject,
    InputMultiPath,
    Undefined,
    Str,
)
from nipype.external.due import BibTeX
from nipype.interfaces.afni.base import (
    AFNICommandBase,
    AFNICommand,
    AFNICommandInputSpec,
    AFNICommandOutputSpec,
    AFNIPythonCommandInputSpec,
    AFNIPythonCommand,
)
from nipype.interfaces.afni.utils import GCORInputSpec, GCOROutputSpec


class GCOR2(CommandLine):
    """
    Computes the average correlation between every voxel
    and ever other voxel, over any give mask.


    For complete details, see the `@compute_gcor Documentation.
    <https://afni.nimh.nih.gov/pub/dist/doc/program_help/@compute_gcor.html>`_

    Examples
    --------
    >>> from nipype.interfaces import afni
    >>> gcor = afni.GCOR()
    >>> gcor.inputs.in_file = 'structural.nii'
    >>> gcor.inputs.nfirst = 4
    >>> gcor.cmdline
    '@compute_gcor -nfirst 4 -input structural.nii'
    >>> res = gcor.run()  # doctest: +SKIP

    """

    _cmd = "echo GCOR = 45.5; #"
    input_spec = GCORInputSpec
    output_spec = GCOROutputSpec

    def _run_interface(self, runtime):
        runtime = super(GCOR2, self)._run_interface(runtime)

        gcor_line = [
            line.strip()
            for line in runtime.stdout.split("\n")
            if line.strip().startswith("GCOR = ")
        ][-1]
        setattr(self, "_gcor", float(gcor_line[len("GCOR = ") :]))
        return runtime

    def _list_outputs(self):
        return {"out": getattr(self, "_gcor")}


class FakeGCOR(GCOR):

    def _run_interface(self, runtime):
        setattr(self, "_gcor", 57.32)
        return runtime


class Add_Floats(ProcessMIA):
    """ajout de flottants pour tester la sortie qui devrait un flottant"""

    def __init__(self):
        """Dedicated to the attributes initialisation/instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # initialisation of the objects needed for the launch of the brick
        super(Add_Floats, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = [] # no need of third party software!

        # Inputs description
        number1_desc = 'A float'
        number2_desc = 'A float'
    

        # Outputs description
        output_number_desc = 'A float'
        
        # Inputs traits
        self.add_trait("number1",
                       traits.Float(output=False,
                                    desc=number1_desc))

        self.add_trait("number2",
                       traits.Float(output=False,
                                    desc=number2_desc))
        

        # Outputs traits
        self.add_trait("output_number",
                       traits.Float(output=True,
                                    desc=output_number_desc))

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
        super(Add_Floats, self).list_outputs()

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        self.output_number = self.number1 + self.number2
        return


