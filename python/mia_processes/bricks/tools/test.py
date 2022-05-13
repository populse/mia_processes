# -*- coding: utf-8 -*- #

"""tests sandbox"""

# nipype import
from nipype.interfaces.base import traits

# populse_mia imports
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

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
        self.out = self.number1 + self.number2
        return


