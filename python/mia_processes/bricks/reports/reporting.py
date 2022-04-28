"""The reporting library of the mia_processes package.

The purpose of this module is to provide the reporting bricks necessary to
generate automatic report at the end of a pipeline calculation.

:Contains:
    :Class:
        - MRIQC_report

    :Function:
        -
"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# nibabel import
import nibabel as nib

# nipype import
from nipype.interfaces.base import (OutputMultiPath, InputMultiPath, File,
                                    traits, TraitListObject, Undefined,
                                    DictStrStr, Str)
from nipype.interfaces.spm.base import ImageFileSPM

# populse_mia import
from populse_mia import sources_images
from populse_mia.data_manager.project import COLLECTION_CURRENT
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA



# soma-base imports
from soma.qt_gui.qt_backend.Qt import QMessageBox


# Other import
import os
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.lib.pagesizes import A4, landscape, portrait
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                 Table, TableStyle, PageBreak)


class MRIQC_report(ProcessMIA):
    """
        * Generates the report for MRIQC pipeline *

        Please, see the complete documentation for the `MRIQC_report brick in the populse.mia_processes website
        https://populse.github.io/mia_processes/documentation/bricks/reports/MRIQC_report.html

        """

    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(MRIQC_report, self).__init__()

        # Third party software required for the execution of the brick
        self.requirement = []

        # Inputs description
        anat_desc = ('An existing, uncompressed anatomic image file (valid '
                       'extensions: .nii)')

        # Outputs description
        report_desc = 'The generated report (pdf)'

        # Inputs traits
        self.add_trait("anat",
                       ImageFileSPM(copyfile=False,
                                    output=False,
                                    optional=False,
                                    desc=anat_desc))
        # Outputs traits
        self.add_trait("report",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=report_desc))

        self.add_trait("patient_name",
                      traits.Either(Undefined,
                                    Str(),
                                    output=False,
                                    optional=False,
                                    userlevel=1))
        self.patient_name = Undefined
        self.init_default_traits()

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dic),
        if defined here. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(MRIQC_report, self).list_outputs()



        if self.anat and self.anat not in ["<undefined>", traits.Undefined]:
            self.outputs['report'] = 'report.pdf'



        file_position = (self.anat.find(self.project.getName())
                         + len(self.project.getName()) + 1)
        database_filename = self.anat[file_position:]

        if ('PatientName' in
                self.project.session.get_fields_names(COLLECTION_CURRENT)):
            self.patient_name = self.project.session.get_value(COLLECTION_CURRENT,
                                                  database_filename,
                                                  ('PatientName'))

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(MRIQC_report, self).run_process_mia()


        print('self.patient_name: ', self.patient_name)
        header_image = os.path.join(sources_images.__path__[0],
                                    'Logo_populse_mia_LR.jpeg')







        title = '<font size=18> <b> %s report: </b> %s <br/> <br/></font>' % (patient_ref, aujourdhui.strftime('%Y.%m.%d'))     # Title:  patient_ref + report - running date.
        #Canvas creation
        rapport = []

        # First page - cover

        im_entete.hAlign = 'CENTER'
        rapport.append(im_entete)

        rapport.append(Spacer(0 * mm, 10 * mm))  # (width, height)

        rapport.append(Paragraph(nomDuSuperLogiciel,
                                 styles['Center']))

        rapport.append(Spacer(0 * mm, 4 * mm))  # (width, height)  ## In Linux
        # rapport.append(Spacer(0*mm, 4*mm)) # (width, height)  ## In Mac

        line = ReportLine(250)
        line.hAlign = 'CENTER'
        rapport.append(line)

        rapport.append(Spacer(0 * mm, 20 * mm))

        rapport.append(Paragraph(titre,
                                 styles['Center']))

        rapport.append(Spacer(0 * mm, 10 * mm))

        data.setStyle(TableStyle([
            ('LINEABOVE', (0, 0), (-1, 0), 2, colors.black),
            ('LINEBELOW', (0, -1), (-1, -1), 2, colors.black),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE')
        ]))  # Carefull: Table use (col,raw) and 0-base for top left start as usual OR -1-base for lower right start.

        data.hAlign = 'CENTER'
        rapport.append(data)

        rapport.append(Spacer(0 * mm, 10 * mm))

        rapport.append(Paragraph(texteDisclaimer,
                                 styles["Justify"]))

        rapport.append(PageBreak())