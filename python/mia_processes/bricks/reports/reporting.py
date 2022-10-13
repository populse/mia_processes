"""The reporting library of the mia_processes package.

The purpose of this module is to provide the reporting bricks necessary to
generate automatic report at the end of a pipeline calculation.

:Contains:
    :Class:
        - MRIQC_anat_report

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
#from nibabel.processing import resample_from_to

# nipype import
from nipype import info as nipype_info
from nipype.interfaces.base import (DictStrStr, File, InputMultiPath,
                                    OutputMultiPath, Str, traits,
                                    TraitListObject, Undefined,)
from nipype.interfaces.spm.base import ImageFileSPM

# capsul import
from capsul import info as capsul_info

# mia_processes import:
from mia_processes import info as mia_processes_info
from mia_processes.utils import (dict4runtime_update, PageNumCanvas, plot_qi2,
                                 recupCover, ReportLine, slice_planes_plot)
from mia_processes.utils import Report

# populse_mia import
from populse_mia import info as mia_info
from populse_mia import sources_images
from populse_mia.data_manager.project import COLLECTION_CURRENT
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# soma-base imports
from soma.qt_gui.qt_backend.Qt import QMessageBox


# Other import
import json
import numpy as np
import os
import platform
import tempfile
from datetime import datetime

from sys import exit, path, version

#import matplotlib.pyplot as plt
#import scipy.ndimage as ndi



class MRIQC_anat_report(ProcessMIA):
    """
        * Generates the report for anatomical data in MRIQC pipeline *

        Please, see the complete documentation for the `MRIQC_report brick in the populse.mia_processes website
        https://populse.github.io/mia_processes/documentation/bricks/reports/MRIQC_report.html

        """

    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(MRIQC_anat_report, self).__init__()

        # Third party software required for the execution of the brick
        self.requirement = []

        # Inputs description

        IQMs_file_desc = 'A .JSON file containing the IQMs'

        anat_desc = ('An existing, uncompressed anatomical image file (valid '
                     'extensions: .nii)')

        anat_fig_rows_desc = ('The number of lines for the anatomical slice '
                              'planes plot')

        anat_fig_cols_desc = ('The number of columns for the anatomical slice '
                              'planes plot')

        anat_inf_slice_start_desc = ('The first index displayed in anatomical '
                                     'slice planes plot')

        anat_slices_gap_desc = ('Gap between slices in anatomical slice planes '
                                'plot')

        norm_anat_desc = ('An existing, uncompressed normalised anatomical '
                          'image file (valid extensions: .nii)')

        norm_anat_fig_rows_desc = ('The number of lines for the normalised '
                                   'anatomical slice planes plot')

        norm_anat_fig_cols_desc = ('The number of columns for the normalised '
                                   'anatomical slice planes plot')

        norm_anat_inf_slice_start_desc = ('The first index displayed in '
                                          'normalised anatomical slice planes '
                                          'plot')

        norm_anat_slices_gap_desc = ('Gap between slices in normalised '
                                     'anatomical slice planes plot')


        # Outputs description
        report_desc = 'The generated report (pdf)'

        # Inputs traits
        self.add_trait("IQMs_file",
                       File(copyfile=False,
                            output=False,
                            optional=False,
                            desc=IQMs_file_desc))

        self.add_trait("anat",
                       ImageFileSPM(copyfile=False,
                                    output=False,
                                    optional=False,
                                    desc=anat_desc))

        self.add_trait("anat_fig_rows",
                       traits.Int(5,
                                  output=False,
                                  optional=True,
                                  desc=anat_fig_rows_desc))

        self.add_trait("anat_fig_cols",
                       traits.Int(5,
                                  output=False,
                                  optional=True,
                                  desc=anat_fig_cols_desc))

        self.add_trait("anat_inf_slice_start",
                       traits.Either(Undefined,
                                     traits.Int,
                                     default=Undefined,
                                     output=False,
                                     optional=True,
                                     desc=anat_inf_slice_start_desc))

        self.add_trait("anat_slices_gap",
                       traits.Either(Undefined,
                                     traits.Int,
                                     default=Undefined,
                                     output=False,
                                     optional=True,
                                     desc=anat_slices_gap_desc))

        self.add_trait("norm_anat",
                       ImageFileSPM(copyfile=False,
                                    output=False,
                                    optional=False,
                                    desc=norm_anat_desc))

        self.add_trait("norm_anat_fig_rows",
                       traits.Int(5,
                                  output=False,
                                  optional=True,
                                  desc=norm_anat_fig_rows_desc))

        self.add_trait("norm_anat_fig_cols",
                       traits.Int(5,
                                  output=False,
                                  optional=True,
                                  desc=norm_anat_fig_cols_desc))

        self.add_trait("norm_anat_inf_slice_start",
                       traits.Either(Undefined,
                                     traits.Int,
                                     default=Undefined,
                                     output=False,
                                     optional=True,
                                     desc=norm_anat_inf_slice_start_desc))

        self.add_trait("norm_anat_slices_gap",
                       traits.Either(Undefined,
                                     traits.Int,
                                     default=Undefined,
                                     output=False,
                                     optional=True,
                                     desc=norm_anat_slices_gap_desc))

        # Outputs traits
        self.add_trait("report",
                       OutputMultiPath(File(),
                                       output=True,
                                       optional=True,
                                       desc=report_desc))
        # Special parameter used as a messenger for the run_process_mia method
        self.add_trait("dict4runtime",
                       traits.Dict(output=False,
                                   optional=True,
                                   userlevel=1))

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
        super(MRIQC_anat_report, self).list_outputs()

        file_position = (self.anat.find(self.project.getName())
                         + len(self.project.getName()) + 1)
        database_filename = self.anat[file_position:]

        # As we do not have access to the database at the runtime (see #272),
        # we prepare here the data that the run_process_mia method will need
        # via dict4runtime parameter:
        dict4runtime_update(self.dict4runtime, self.project.session,
                            database_filename,
                            'PatientName', 'StudyName', 'AcquisitionDate',
                            'Sex', 'Site', 'Spectro', 'Age')

        # FIXME: Currently, Site and Spectro data is hard-coded. A solution
        #        should be found to retrieve them automatically or to put them
        #        in the input parameters of the brick:
        # Site
        if self.dict4runtime['Site'] in ("", "Undefined"):
            self.dict4runtime['Site'] = 'Grenoble University Hospital - CLUNI'

        # MriScanner
        if self.dict4runtime['Spectro'] in ("", "Undefined"):
            self.dict4runtime['Spectro'] = 'Philips Achieva 3.0T TX'

        # Generate an output name
        if self.anat and self.anat not in ["<undefined>", traits.Undefined]:
            self.outputs['report'] = os.path.join(
                        self.output_directory,
                        "{0}_mriqcReport_{1}.pdf".format(
                            self.dict4runtime['PatientName'] if
                            self.dict4runtime['PatientName'] != "Undefined" else
                            "Undefined_name_ref",
                            datetime.now().strftime('%Y_%m_%d_'
                                                    '%H_%M_%S_%f')[:22]))

        # FIXME: Do we need tags inheritance ?

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(MRIQC_anat_report, self).run_process_mia()

        report = Report(self.dict4runtime, self.IQMs_file, self.anat,
                        self.norm_anat, self.report)
        report.make_report()







        # First page - cover ##################################################
        ######################################################################

        # Third page - IQMs ###################################################
        #######################################################################
        report.append(Paragraph("<font size = 18 ><b>Image parameters"
                                "</b></font>",
                                styles['Center']))
        report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)
        report.append(Spacer(0 * mm, 20 * mm))
        ### Spatial distribution ##############################################
        report.append(Paragraph("<font size = 15 ><b>SPATIAL DISTRIBUTION"
                                "</b></font>",
                                styles['Bullet1']))
        report.append(Spacer(0 * mm, 2 * mm))
        report.append(Paragraph("Information theory to evaluate the spatial "
                                "distribution of information",
                                styles['Left2']))
        report.append(Spacer(0 * mm, 10 * mm))
        # fber
        try:
            report.append(Paragraph(
                            '<font size = 9><sup>%</sup></font><font size = 11>'
                            '<b>FBER</b></font>: ' +
                            str(round(data.get('fber'), 2)),
                            styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                            '<font size = 9><sup>%</sup></font><font size = 11>'
                            '<b>FBER</b></font>: ' +
                            'Not determined',
                            styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))
        # efc
        try:
            report.append(Paragraph(
                            '<font size = 9><sup>#</sup></font><font size = 11>'
                            '<b>EFC</b></font>: ' +
                            str(round(data.get('efc'), 2)),
                            styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                            '<font size = 9><sup>#</sup></font><font size = 11>'
                            '<b>EFC</b></font>: ' +
                            'Not determined',
                            styles['Bullet2']))

        report.append(Spacer(0 * mm, 20 * mm))
        ### Artifacts #########################################################
        report.append(Paragraph("<font size = 15 ><b>ARTIFACTS</b></font>",
                                styles['Bullet1']))
        report.append(Spacer(0 * mm, 2 * mm))
        report.append(Paragraph("Estimates artefacts and signal leakage due to "
                                "rapid movement (e.g. eye movement or blood "
                                "vessel pulsation, etc.)",
                                styles['Left2']))

        report.append(Spacer(0 * mm, 10 * mm))
        # wm2max
        try:
            report.append(Paragraph(
                       '<font size = 9><sup>&</sup></font><font size = 11><b>'
                       'White-matter to maximum intensity ratio </b></font>: ' +
                       str(round(data.get('wm2max'), 2)),
                       styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                       '<font size = 9><sup>&</sup></font><font size = 11><b>'
                       'White-matter to maximum intensity ratio </b></font>: ' +
                       'Not determined',
                       styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))
        # qi_1
        try:
            report.append(Paragraph(
                         '<font size = 9><sup>$</sup></font><font size = 11><b>'
                         'Mortamet’s quality index 1</b></font>: ' +
                         '{:.2e}'.format(data['qi_1']),
                         styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                         '<font size = 9><sup>$</sup></font><font size = 11><b>'
                         'Mortamet’s quality index 1</b></font>: ' +
                         'Not determined',
                         styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))
        # inu_range
        try:
            report.append(Paragraph(
                 '<font size = 9><sup>*</sup></font><font size = 11><b>Bias '
                 'field range (95th percentile - 5th percentile)</b></font>: ' +
                 str(round(data.get('inu_range'), 2)),
                 styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                 '<font size = 9><sup>*</sup></font><font size = 11><b>Bias '
                 'field range (95th percentile - 5th percentile)</b></font>: ' +
                 'Not determined',
                 styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # inu_med
        try:
            report.append(Paragraph(
                    '<font size = 9><sup>*</sup></font><font size = 11><b>Bias '
                    'field median</b></font>: ' +
                    str(round(data.get('inu_med'), 2)),
                    styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                    '<font size = 9><sup>*</sup></font><font size = 11><b>Bias '
                    'field median</b></font>: ' +
                    'Not determined',
                    styles['Bullet2']))

        report.append(Spacer(0 * mm, 100 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        report.append(line)
        report.append(Spacer(0 * mm, 2.5 * mm))
        report.append(Paragraph(
             "<font size = 8><sup>%</sup>Shehzad Z et al., <i>The Preprocessed "
             "Connectomes Project Quality Assessment Protocol - a resource for "
             "measuring the quality of MRI data</i>, Front. Neurosci. "
             "Conference Abstract: Neuroinformatics 2015. Higher values are "
             "better.</font>",
             styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>#</sup>Atkinson et al., <i>Automatic "
            "correction of motion artifacts in magnetic resonance images using "
            "an entropy focus criterion</i>, IEEE Trans Med Imag 16(6):903-910,"
            " 1997. Lower values are better.</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>&</sup>The median intensity within the white "
            "matter mask over the 95% percentile of the full intensity "
            "distribution. Values should be around the interval "
            "[0.6, 0.8].</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>$</sup>Mortamet B et al., <i>Automatic quality"
            " assessment in structural brain magnetic resonance imaging</i>, "
            "Mag Res Med 62(2):365-372, 2009. Lower values are better.</font>",
            styles['Left']))
        report.append(Paragraph(
            "<font size = 8><sup>*</sup>Tustison NJ et al., <i>N4ITK: improved "
            "N3 bias correction</i>, IEEE Trans Med Imag, 29(6):1310-20, 2010. "
            "Median closer to 1 and range closer to 0 are better.</font>",
            styles['Left']))
        report.append(PageBreak())

        # fourth page - IQMs ##################################################
        #######################################################################
        report.append(Paragraph("<font size = 18 ><b>Image parameters"
                                "</b></font>",
                                styles['Center']))
        report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)
        report.append(Spacer(0 * mm, 20 * mm))
        ### Tissues Quality ###################################################
        report.append(Paragraph("<font size = 15 ><b>TISSUES QUALITY</b></font>",
                                styles['Bullet1']))
        report.append(Spacer(0 * mm, 2 * mm))
        report.append(Paragraph("Metrics that do not fall into the above "
                                "categories: statistical properties of tissue "
                                "distributions, volume overlap of tissues, "
                                "image harpness/blurriness, etc.",
                                styles['Left2']))
        report.append(Spacer(0 * mm, 10 * mm))
        # summary_csf_mean
        try:
            report.append(Paragraph(
                 '<font size = 11><b>Mean of the distribution of cerebrospinal '
                 'fluid</b></font>: ' +
                 str(round(data.get('summary_csf_mean'), 2)),
                 styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                 '<font size = 11><b>Mean of the distribution of cerebrospinal '
                 'fluid</b></font>: ' +
                 'Not determined',
                 styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_stdv
        try:
            report.append(Paragraph(
                 '<font size = 11><b>Standard deviation of the distribution of '
                 'cerebrospinal fluid</b></font>: ' +
                 str(round(data.get('summary_csf_stdv'), 2)),
                 styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                 '<font size = 11><b>Standard deviation of the distribution of '
                 'cerebrospinal fluid</b></font>: ' +
                 'Not determined',
                 styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_median
        try:
            report.append(Paragraph(
                             '<font size = 11><b>Median of the distribution of '
                             'cerebrospinal fluid</b></font>: ' +
                             str(round(data.get('summary_csf_median'), 2)),
                             styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                             '<font size = 11><b>Median of the distribution of '
                             'cerebrospinal fluid</b></font>: ' +
                             'Not determined',
                             styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_mad
        try:
            report.append(Paragraph(
                          '<font size = 11><b>Median absolute deviation of the '
                          'distribution of cerebrospinal fluid</b></font>: ' +
                          str(round(data.get('summary_csf_mad'), 2)),
                          styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                          '<font size = 11><b>Median absolute deviation of the '
                          'distribution of cerebrospinal fluid</b></font>: ' +
                          'Not determined',
                          styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_p95
        try:
            report.append(Paragraph(
                     '<font size = 11><b>95% percentile of the distribution of '
                     'cerebrospinal fluid</b></font>: ' +
                     str(round(data.get('summary_csf_p95'), 2)),
                     styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                     '<font size = 11><b>95% percentile of the distribution of '
                     'cerebrospinal fluid</b></font>: ' +
                     'Not determined',
                     styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_p05
        try:
            report.append(Paragraph(
                      '<font size = 11><b>5% percentile of the distribution of '
                      'cerebrospinal fluid</b></font>: ' +
                      str(round(data.get('summary_csf_p05'), 2)),
                      styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                      '<font size = 11><b>5% percentile of the distribution of '
                      'cerebrospinal fluid</b></font>: ' +
                      'Not determined',
                      styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_k
        try:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis '
                'of the distribution of cerebrospinal fluid</b></font>: ' +
                str(round(data.get('summary_csf_k'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis '
                'of the distribution of cerebrospinal fluid</b></font>: ' +
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_csf_n
        try:
            report.append(Paragraph(
                   '<font size = 11><b>Number of voxels in the distribution of '
                   'cerebrospinal fluid</b></font>: ' +
                   '{:.2e}'.format(data['summary_csf_n']),
                   styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                   '<font size = 11><b>Number of voxels in the distribution of '
                   'cerebrospinal fluid</b></font>: ' +
                   'Not determined',
                   styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))
        # summary_gm_mean
        try:
            report.append(Paragraph(
                          '<font size = 11><b>Mean of the distribution of gray '
                          'matter</b></font>: ' +
                          str(round(data.get('summary_gm_mean'), 2)),
                          styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                          '<font size = 11><b>Mean of the distribution of gray '
                          'matter</b></font>: ' +
                          'Not determined',
                          styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_stdv
        try:
            report.append(Paragraph(
                    '<font size = 11><b>Standard deviation of the distribution '
                    'of gray matter</b></font>: ' +
                    str(round(data.get('summary_gm_stdv'), 2)),
                    styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                    '<font size = 11><b>Standard deviation of the distribution '
                    'of gray matter</b></font>: ' +
                    'Not determined',
                    styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_median
        try:
            report.append(Paragraph(
                             '<font size = 11><b>Median of the distribution of '
                             'gray matter</b></font>: ' +
                             str(round(data.get('summary_gm_median'), 2)),
                             styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                             '<font size = 11><b>Median of the distribution of '
                             'gray matter</b></font>: ' +
                             'Not determined',
                             styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_mad
        try:
            report.append(Paragraph(
                          '<font size = 11><b>Median absolute deviation of the '
                          'distribution of gray matter</b></font>: ' +
                          str(round(data.get('summary_gm_mad'), 2)),
                          styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                          '<font size = 11><b>Median absolute deviation of the '
                          'distribution of gray matter</b></font>: ' +
                          'Not determined',
                          styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_p95
        try:
            report.append(Paragraph(
                     '<font size = 11><b>95% percentile of the distribution of '
                     'gray matter</b></font>: ' +
                     str(round(data.get('summary_gm_p95'), 2)),
                     styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                     '<font size = 11><b>95% percentile of the distribution of '
                     'gray matter</b></font>: ' +
                     'Not determined',
                     styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_p05
        try:
            report.append(Paragraph(
                         '<font size = 11><b>5% percentile of the distribution '
                         'of gray matter</b></font>: ' +
                         str(round(data.get('summary_gm_p05'), 2)),
                         styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                         '<font size = 11><b>5% percentile of the distribution '
                         'of gray matter</b></font>: ' +
                         'Not determined',
                         styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_k
        try:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis '
                'of the distribution of gray matter</b></font>: ' +
                str(round(data.get('summary_gm_k'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis '
                'of the distribution of gray matter</b></font>: ' +
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_gm_n
        try:
            report.append(Paragraph(
                   '<font size = 11><b>Number of voxels in the distribution of '
                   'gray matter</b></font>: ' +
                   '{:.2e}'.format(data['summary_gm_n']),
                   styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                   '<font size = 11><b>Number of voxels in the distribution of '
                   'gray matter</b></font>: ' +
                   'Not determined',
                   styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))
        # summary_wm_mean
        try:
            report.append(Paragraph(
                               '<font size = 11><b>Mean of the distribution of '
                               'white matter</b></font>: ' +
                               str(round(data.get('summary_wm_mean'), 2)),
                               styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                               '<font size = 11><b>Mean of the distribution of '
                               'white matter</b></font>: ' +
                               'Not determined',
                               styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_stdv
        try:
            report.append(Paragraph(
                    '<font size = 11><b>Standard deviation of the distribution '
                    'of white matter</b></font>: ' +
                    str(round(data.get('summary_wm_stdv'), 2)),
                    styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                    '<font size = 11><b>Standard deviation of the distribution '
                    'of white matter</b></font>: ' +
                    'Not determined',
                    styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_median
        try:
            report.append(Paragraph(
                                '<font size = 11><b>Median of the distribution '
                                'of white matter</b></font>: ' +
                                str(round(data.get('summary_wm_median'), 2)),
                                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                                '<font size = 11><b>Median of the distribution '
                                'of white matter</b></font>: ' +
                                'Not determined',
                                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_mad
        try:
            report.append(Paragraph(
                          '<font size = 11><b>Median absolute deviation of the '
                          'distribution of white matter</b></font>: ' +
                          str(round(data.get('summary_wm_mad'), 2)),
                          styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                          '<font size = 11><b>Median absolute deviation of the '
                          'distribution of white matter</b></font>: ' +
                          'Not determined',
                          styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_p95
        try:
            report.append(Paragraph(
                        '<font size = 11><b>95% percentile of the distribution '
                        'of white matter</b></font>: ' +
                        str(round(data.get('summary_wm_p95'), 2)),
                        styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                        '<font size = 11><b>95% percentile of the distribution '
                        'of white matter</b></font>: ' +
                        'Not determined',
                        styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_p05
        try:
            report.append(Paragraph(
                         '<font size = 11><b>5% percentile of the distribution '
                         'of white matter</b></font>: ' +
                         str(round(data.get('summary_wm_p05'), 2)),
                         styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                         '<font size = 11><b>5% percentile of the distribution '
                         'of white matter</b></font>: ' +
                         'Not determined',
                         styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_k
        try:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis '
                'of the distribution of white matter</b></font>: ' +
                str(round(data.get('summary_wm_k'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis '
                'of the distribution of white matter</b></font>: ' +
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_wm_n
        try:
            report.append(Paragraph(
                   '<font size = 11><b>Number of voxels in the distribution of '
                   'white matter</b></font>: ' +
                   '{:.2e}'.format(data['summary_wm_n']),
                   styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                   '<font size = 11><b>Number of voxels in the distribution of '
                   'white matter</b></font>: ' +
                   'Not determined',
                   styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))
        # summary_bg_mean
        try:
            report.append(Paragraph(
                               '<font size = 11><b>Mean of the distribution of '
                               'background</b></font>: ' +
                               str(round(data.get('summary_bg_mean'), 2)),
                               styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                               '<font size = 11><b>Mean of the distribution of '
                               'background</b></font>: ' +
                               'Not determined',
                               styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_stdv
        try:
            report.append(Paragraph(
                 '<font size = 11><b>Standard deviation of the distribution of '
                 'background</b></font>: ' +
                 str(round(data.get('summary_bg_stdv'), 2)),
                 styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                 '<font size = 11><b>Standard deviation of the distribution of '
                 'background</b></font>: ' +
                 'Not determined',
                 styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_median
        try:
            report.append(Paragraph(
                             '<font size = 11><b>Median of the distribution of '
                             'background</b></font>: ' +
                             str(round(data.get('summary_bg_median'), 2)),
                             styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                             '<font size = 11><b>Median of the distribution of '
                             'background</b></font>: ' +
                             'Not determined',
                             styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_mad
        try:
            report.append(Paragraph(
                          '<font size = 11><b>Median absolute deviation of the '
                          'distribution of background</b></font>: ' +
                          str(round(data.get('summary_bg_mad'), 2)),
                          styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                          '<font size = 11><b>Median absolute deviation of the '
                          'distribution of background</b></font>: ' +
                          'Not determined',
                          styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_p95
        try:
            report.append(Paragraph(
                     '<font size = 11><b>95% percentile of the distribution of '
                     'background</b></font>: ' +
                     str(round(data.get('summary_bg_p95'), 2)),
                     styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                     '<font size = 11><b>95% percentile of the distribution of '
                     'background</b></font>: ' +
                     'Not determined',
                     styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_p05
        try:
            report.append(Paragraph(
                      '<font size = 11><b>5% percentile of the distribution of '
                      'background</b></font>: ' +
                      str(round(data.get('summary_bg_p05'), 2)),
                      styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                      '<font size = 11><b>5% percentile of the distribution of '
                      'background</b></font>: ' +
                      'Not determined',
                      styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_k
        try:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis '
                'of the distribution of background</b></font>: ' +
                str(round(data.get('summary_bg_k'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>*</sup></font><font size = 11><b>Kurtosis '
                'of the distribution of background</b></font>: ' +
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # summary_bg_n
        try:
            report.append(Paragraph(
                   '<font size = 11><b>Number of voxels in the distribution of '
                   'background</b></font>: ' +
                   '{:.2e}'.format(data['summary_bg_n']),
                   styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                   '<font size = 11><b>Number of voxels in the distribution of '
                   'background</b></font>: ' +
                   'Not determined',
                   styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))
        # fwhm_x
        try:
            report.append(Paragraph(
                 '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of '
                 'the distribution in units of voxels, along x '
                 'dimension</b></font>: ' +
                 str(round(data.get('fwhm_x'), 2)),
                 styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                 '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of '
                 'the distribution in units of voxels, along x '
                 'dimension</b></font>: ' +
                 'Not determined',
                 styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_y
        try:
            report.append(Paragraph(
                 '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of '
                 'the distribution in units of voxels, along y '
                 'dimension</b></font>: ' +
                 str(round(data.get('fwhm_y'), 2)),
                 styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                 '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of '
                 'the distribution in units of voxels, along y '
                 'dimension</b></font>: ' +
                 'Not determined',
                 styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_z
        try:
            report.append(Paragraph(
                 '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of '
                 'the distribution in units of voxels, along z '
                 'dimension</b></font>: ' +
                 str(round(data.get('fwhm_z'), 2)),
                 styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                 '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of '
                 'the distribution in units of voxels, along z '
                 'dimension</b></font>: ' +
                 'Not determined',
                 styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # fwhm_avg
        try:
            report.append(Paragraph(
                 '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of '
                 'the distribution in units of voxels, average</b></font>: ' +
                 str(round(data.get('fwhm_avg'), 2)),
                 styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                 '<font size = 9><sup>#</sup></font><font size = 11><b>FWHM of '
                 'the distribution in units of voxels, average</b></font>: ' +
                 'Not determined',
                 styles['Bullet2']))

        report.append(Spacer(0 * mm, 7 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        report.append(line)
        report.append(Spacer(0 * mm, 2.5 * mm))
        report.append(Paragraph(
            "<font size = 8><sup>*</sup>K is always ≥ -2. If the distribution "
            "is Gaussian, K = 0. If a distribution has less weight on its "
            "center and tails compared to a Gaussian of the same variance, "
            "then K < 0. If the distribution has more weight on its center and "
            "tails, then K > 0.</font>",
            styles['Left']))
        report.append(Paragraph(
              "<font size = 8><sup>#</sup>Forman SD et al., <i>Improved "
              "assessment of significant activation in functional magnetic "
              "resonance imaging(fMRI): use of a cluster - size threshold</i>, "
              "Magn.Reson.Med. 33(5), 636–647, 1995. Lower values are better, "
              "higher values indicate a blurrier image.</font>",
              styles['Left']))
        report.append(PageBreak())

        # Fifth page page - IQMs###############################################
        #######################################################################
        report.append(Paragraph("<font size = 18 ><b>Image parameters</b>"
                                "</font>",
                                styles['Center']))
        report.append(Spacer(0 * mm, 4 * mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)
        report.append(Spacer(0 * mm, 20 * mm))
        ### Tissues Quality ###################################################
        report.append(Paragraph("<font size = 15 ><b>TISSUES QUALITY</b>"
                                "</font>",
                                styles['Bullet1']))
        report.append(Spacer(0 * mm, 2 * mm))
        report.append(Paragraph(
              "Metrics that do not fall into the above categories: statistical "
              "properties of tissue distributions, volume overlap of tissues, "
              "image harpness/blurriness, etc.",
              styles['Left2']))
        report.append(Spacer(0 * mm, 10 * mm))
        # icvs_csf
        try:
            report.append(Paragraph(
                           '<font size = 11><b>Intracranial volume fraction of '
                           'cerebrospinal fluid</b></font>: ' +
                           str(round(data.get('icvs_csf'), 2)),
                           styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                           '<font size = 11><b>Intracranial volume fraction of '
                           'cerebrospinal fluid</b></font>: ' +
                           'Not determined',
                           styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # icvs_gm
        try:
            report.append(Paragraph(
                           '<font size = 11><b>Intracranial volume fraction of '
                           'gray matter</b></font>: ' +
                           str(round(data.get('icvs_gm'), 2)),
                           styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                           '<font size = 11><b>Intracranial volume fraction of '
                           'gray matter</b></font>: ' +
                           'Not determined',
                           styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # icvs_wm
        try:
            report.append(Paragraph(
                           '<font size = 11><b>Intracranial volume fraction of '
                           'white matter</b></font>: ' +
                           str(round(data.get('icvs_wm'), 2)),
                           styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                           '<font size = 11><b>Intracranial volume fraction of '
                           'white matter</b></font>: ' +
                           'Not determined',
                           styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))
        # rpve_csf
        try:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Residual '
                'partial voluming error of cerebrospinal fluid</b></font>: ' +
                str(round(data.get('rpve_csf'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Residual '
                'partial voluming error of cerebrospinal fluid</b></font>: ' +
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # rpve_gm
        try:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Residual '
                'partial voluming error of gray matter</b></font>: ' +
                str(round(data.get('rpve_gm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Residual '
                'partial voluming error of gray matter</b></font>: ' +
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # rpve_wm
        try:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Residual '
                'partial voluming error of white matter</b></font>: ' +
                str(round(data.get('rpve_wm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                '<font size = 9><sup>%</sup></font><font size = 11><b>Residual '
                'partial voluming error of gray matter</b></font>: ' +
                'Not determined',
                styles['Bullet2']))

        report.append(Spacer(0 * mm, 2.5 * mm))
        # tpm_overlap_csf
        try:
            report.append(Paragraph(
                 '<font size = 9><sup>&</sup></font><font size = 11><b>Overlap '
                 'of the tissue probability map and the corresponding ICBM '
                 'nonlinear-asymmetric 2009c template, for cerebrospinal '
                 'fluid</b></font>: ' +
                 str(round(data.get('tpm_overlap_csf'), 2)),
                 styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                 '<font size = 9><sup>&</sup></font><font size = 11><b>Overlap '
                 'of the tissue probability map and the corresponding ICBM '
                 'nonlinear-asymmetric 2009c template, for cerebrospinal '
                 'fluid</b></font>: ' +
                 'Not determined',
                 styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # tpm_overlap_gm
        try:
            report.append(Paragraph(
                 '<font size = 9><sup>&</sup></font><font size = 11><b>Overlap '
                 'of the tissue probability map and the corresponding ICBM '
                 'nonlinear-asymmetric 2009c template, '
                 'for gray matter</b></font>: ' +
                 str(round(data.get('tpm_overlap_gm'), 2)),
                 styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                 '<font size = 9><sup>&</sup></font><font size = 11><b>Overlap '
                 'of the tissue probability map and the corresponding ICBM '
                 'nonlinear-asymmetric 2009c template, '
                 'for gray matter</b></font>: ' +
                 'Not determined',
                 styles['Bullet2']))

        report.append(Spacer(0 * mm, 1 * mm))
        # tpm_overlap_wm
        try:
            report.append(Paragraph(
                 '<font size = 9><sup>&</sup></font><font size = 11><b>Overlap '
                 'of the tissue probability map and the corresponding ICBM '
                 'nonlinear-asymmetric 2009c template, '
                 'for  white matter</b></font>: ' +
                 str(round(data.get('tpm_overlap_wm'), 2)),
                styles['Bullet2']))

        except TypeError:
            report.append(Paragraph(
                 '<font size = 9><sup>&</sup></font><font size = 11><b>Overlap '
                 'of the tissue probability map and the corresponding ICBM '
                 'nonlinear-asymmetric 2009c template, '
                 'for  white matter</b></font>: ' +
                 'Not determined',
                 styles['Bullet2']))

        report.append(Spacer(0 * mm, 147 * mm))
        # Footnote
        line = ReportLine(500)
        line.hAlign = 'CENTER'
        report.append(line)
        report.append(Spacer(0 * mm, 2.5 * mm))
        report.append(Paragraph("<font size = 8><sup>%</sup>Lower values "
                                "are better.</font>",
                                styles['Left']))
        report.append(Paragraph("<font size = 8><sup>&</sup>Higher values "
                                "are better.</font>",
                                styles['Left']))
        report.append(PageBreak())

        # Sixth page - slice planes display - Raw anatomic ####################
        #######################################################################
        report.append(Paragraph("<font size = 18 > <b>MRI axial slice planes "
                                "display</b> </font>",
                                styles['Center']))
        report.append(Spacer(0*mm, 4*mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)
        report.append(Spacer(0 * mm, 10 * mm))
        report.append(Paragraph("<font size = 14 >Raw anatomic images</font>",
                                styles['Center']))
        report.append(Spacer(0*mm, 10*mm))
        report.append(Paragraph('<font size = 9 ><i>"Neurological" '
                                'convention, the left side of the image '
                                'corresponds to the left side of the '
                                'brain.</i></font>',
                                styles['Center']))
        report.append(Spacer(0*mm, 1*mm))
        tmpdir = tempfile.TemporaryDirectory()
        slices_image = slice_planes_plot(
                                      self.anat, self.anat_fig_rows,
                                      self.anat_fig_cols,
                                      inf_slice_start=self.anat_inf_slice_start,
                                      slices_gap=self.anat_slices_gap,
                                      cmap="Greys_r", out_dir=tmpdir.name)
        # reminder: A4 == 210mmx297mm
        slices_image = Image(slices_image,
                             width=7.4803 * inch, height=9.0551 * inch)
        slices_image.hAlign = 'CENTER'
        report.append(slices_image)
        report.append(PageBreak())

        # Seventh page - slice planes display - Normalised anatomic (MNI) #####
        #######################################################################
        report.append(Paragraph("<font size = 18 > <b>MRI axial slice planes"
                                "display</b> </font>",
                                styles['Center']))
        report.append(Spacer(0*mm, 4*mm))
        line = ReportLine(150)
        line.hAlign = 'CENTER'
        report.append(line)
        report.append(Spacer(0 * mm, 10 * mm))
        report.append(Paragraph("<font size = 14 >Normalised anatomic "
                                "(MNI)</font>",
                                styles['Center']))

        report.append(Spacer(0*mm, 10*mm))
        report.append(Paragraph('<font size = 9 > <i> "Neurological" '
                                'convention, the left side of the image '
                                'corresponds to the left side of the '
                                'brain. </i> <br/> </font>',
                                styles['Center']))
        report.append(Spacer(0*mm, 1*mm))
        slices_image = slice_planes_plot(
                                 self.norm_anat,
                                 self.norm_anat_fig_rows,
                                 self.norm_anat_fig_cols,
                                 inf_slice_start=self.norm_anat_inf_slice_start,
                                 slices_gap=self.norm_anat_slices_gap,
                                 cmap="Greys_r", out_dir=tmpdir.name)
        # Currently, we make and save (in derived_data) the qi2 graph, but we
        # don't include it in the report because the result with the test data
        # looks strange.
        _ = plot_qi2(np.asarray(data['histogram_qi2_x_grid']),
                     np.asarray(data['histogram_qi2_ref_pdf']),
                     np.asarray(data['histogram_qi2_fit_pdf']),
                     np.asarray(data['histogram_qi2_ref_data']),
                     int(data['histogram_qi2_cutoff_idx']),
                     out_file=os.path.join(
                         os.path.split(os.path.split(self.anat)[0])[0],
                         'derived_data', 'qi2_plot.svg'))
        # reminder: A4 == 210mmx297mm
        slices_image = Image(slices_image, width=7.4803 * inch, height=9.0551 * inch)
        slices_image.hAlign = 'CENTER'
        report.append(slices_image)
        page.build(report, canvasmaker=PageNumCanvas)
        tmpdir.cleanup()
