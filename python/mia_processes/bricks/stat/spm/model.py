# -*- coding: utf-8 -*- #

"""The library for the SPM fMRI statistical analysis of the mia_processes
package.

The purpose of this module is to customise the main spm statistical analysis
bricks provided by nipype and to correct some things that do not work directly
in populse_mia.

:Contains:
    :Class:
        - EstimateContrast
        - EstimateModel
        - Level1Design

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# mia_processes import
#from .stats_pop_ups import _SessionQuery # currently in ec_dev package

# nibabel import
import nibabel as nib
import nibabel.processing as nibp

# nipype import
from nipype.interfaces import spm
from nipype.interfaces.base import OutputMultiPath, InputMultiPath, File, traits
from nipype.interfaces.spm.base import ImageFileSPM

# populse_db and populse_mia import
from populse_mia.data_manager.project import (COLLECTION_CURRENT,
                                              COLLECTION_INITIAL)
from populse_db.database import FIELD_TYPE_INTEGER
from populse_mia.data_manager.database_mia import (TAG_ORIGIN_USER,
                                                   TAG_UNIT_DEGREE)

# populse_mia imports
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# PyQt5 imports
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QPushButton, QVBoxLayout, QWidget, QHBoxLayout, QLabel, QLineEdit,
    QGroupBox, QMessageBox, QToolButton, QDialog, QDialogButtonBox,
    QApplication)

# soma-base imports
from soma.qt_gui.qt_backend.Qt import QMessageBox

# Other import
import os
import scipy.io
from shutil import copy2
from traits.api import Undefined
import numpy as np


class EstimateContrast(ProcessMIA):
    """Estimate contrasts of interest"""
    def __init__(self):
        """Dedicated to the attributes initialisation / instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(EstimateContrast, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['spm', 'nipype']

        # Inputs description
        spm_mat_file_desc = ('SPM.mat file (a pathlike object or string '
                             'representing a file)')
        session_type_desc = ('Selects the contrast type. One of tcon '
                             '(T-contrast), fcon (F-contrast) or tconsess '
                             '(T-contrast cond/sess based)')
        contrast_name_desc = 'Name of contrasts (a list of string)'
        condition_name_desc = 'Conditions information (a list of list of string)'
        contrast_weight_desc = 'Contrast weights (list of list of float)'
        session_desc = 'Session list (a list of list of float)'
        multi_reg_desc = 'The regressor files(a list of file)'

        # Inputs
        self.add_trait("spm_mat_file",
                       File(output=False,
                            copyfile=True,
                            desc=spm_mat_file_desc))

        self.add_trait("session_type",
                       traits.Enum("tcon",
                                   "fcon",
                                   "tconsess",
                                   output=False,
                                   optional=True,
                                   desc=session_type_desc))

        # self.add_trait("contrast_session_name",
        #                traits.String(output=False,
        #                              optional=False))
        self.add_trait("contrast_name",
                       traits.List(traits.String(),
                                   output=False,
                                   optional=True,
                                   desc=contrast_name_desc))
        self.contrast_name = ['+']

        self.add_trait("condition_name",
                       traits.List(traits.List(traits.String()),
                                   output=False,
                                   optional=True,
                                   desc=condition_name_desc))
        self.condition_name = [['R1_1']]

        self.add_trait("contrast_weight",
                       traits.List(traits.List(traits.Float()),
                                   output=False,
                                   optional=True,
                                   desc=contrast_weight_desc))
        self.contrast_weight = [[1.0]]

        self.add_trait("session",
                       traits.Either(Undefined,
                                     traits.List(traits.Either(
                                                        Undefined,
                                                        traits.List(
                                                            traits.Float()))),
                                     output=False,
                                     optional=True,
                                     desc=session_desc))
        self.session = Undefined

        # self.add_trait("replicate",
        #                traits.Enum("none",
        #                            "repl",
        #                            "replsc",
        #                            "sess",
        #                            "both",
        #                            "bothsc",
        #                            output=False,
        #                            optional=True))

        self.add_trait("multi_reg",
                       traits.Either(Undefined,
                                     traits.List(File()),
                                     output=False,
                                     optional=True,
                                     desc=multi_reg_desc))

        #self.add_trait("multi_reg", traits.List(output=False, optional=True))
        # self.add_trait("contrasts", traits.List([('+', 'T', ['R1_1'], [1])], output=False, optional=True))
        self.add_trait("beta_images", InputMultiPath(File(), output=False, copyfile=False))
        self.add_trait("residual_image", File(output=False, copyfile=False))
        self.add_trait("use_derivs", traits.Bool(output=False, optional=True, xor=['group_contrast']))
        self.add_trait("group_contrast", traits.Bool(output=False, optional=True, xor=['use_derivs']))

        # Outputs
        self.add_trait("con_images", OutputMultiPath(File(), optional=True, output=True))
        self.add_trait("spmT_images", OutputMultiPath(File(), optional=True, output=True))
        self.add_trait("spmF_images", OutputMultiPath(File(), optional=True, output=True))
        self.add_trait("out_spm_mat_file", File(output=True, copyfile=False))

        self.init_default_traits()
        self.init_process('nipype.interfaces.spm.EstimateContrast')

    def _get_contrasts(self, session_type):
        """blabla"""
        #contrast = tuple()
        contrasts = [tuple()]

        # if session_type == 'tcon':
        #     contrast[0] = self.contrast_session_name
        #     contrast[1] = 'T'
        #     contrast[2] = self.contrast_weight
        #     contrast[3] = self.replicate
        #
        # elif session_type == 'fcon':
        #     contrast[0] = self.contrast_session_name
        #     contrast[1] = 'F'
        #     contrast[2] = self.contrast_weight
        #     contrast[3] = self.replicate
        #
        # elif session_type == 'tconsess':
        #     pass
        if session_type == 'tcon':
            stat = 'T'

        elif session_type == 'fcon':
            stat = 'F'

        else:
            stat = None

        if stat is not None:

            if self.session is Undefined:
                contrasts = [(cont_name, stat, condition, cont_weight)
                             for cont_name, condition, cont_weight in
                                 zip(self.contrast_name, self.condition_name,
                                     self.contrast_weight)]
            else:
                contrasts = [(cont_name, stat, condition, cont_weight, sess) if
                             sess is not Undefined else
                             (cont_name, stat, condition, cont_weight)
                             for cont_name, condition, cont_weight, sess in
                                 zip(self.contrast_name, self.condition_name,
                                     self.contrast_weight, self.session)]

        return contrasts

    def list_outputs(self, is_plugged=None):
        """blabla"""
        super(EstimateContrast, self).list_outputs()

        # Old version calling Nipype's method
        '''process = spm.EstimateContrast()
        if not self.spm_mat_file:
            return {}, {}
        else:
            process.inputs.spm_mat_file = self.spm_mat_file

        if not self.contrasts:
            return {}, {}
        else:
            process.inputs.contrasts = self.contrasts

        if not self.beta_images:
            return {}, {}
        else:
            process.inputs.beta_images = self.beta_images

        if not self.residual_image:
            return {}, {}
        else:
            process.inputs.residual_image = self.residual_image

        outputs = process._list_outputs()
        outputs["out_spm_mat_file"] = outputs.pop("spm_mat_file")
        return outputs, {}'''

        # Own list_outputs to avoid to read in the SPM.mat file
        # if not self.spm_mat_file:
        #     return {}, {}
        #
        # if not self.contrasts:
        #     return {}, {}
        #
        # if not self.beta_images:
        #     return {}, {}
        #
        # if not self.residual_image:
        #     return {}, {}

        #if self.outputs:
        #    self.outputs = {}

        #if self.inheritance_dict:
        #    self.inheritance_dict = {}

        if ((self.spm_mat_file) and
                (self.spm_mat_file not in ['<undefined>', Undefined]) and
                          self.multi_reg not in [Undefined, '<undefined>', []]):
            # The management of self.process.output_directory could be delegated
            # to the populse_mia.user_interface.pipeline_manager.process_mia
            # module. We can't do it at the moment because the
            # sync_process_output_traits() of the capsul/process/nipype_process
            # module raises an exception in nipype if the mandatory parameters
            # are not yet defined!
            if self.output_directory:
                self.process.output_directory = self.output_directory

            else:
                print('No output_directory was found...!\n')

            _, spm_mat_file = os.path.split(self.spm_mat_file)
            self.outputs['out_spm_mat_file'] = os.path.join(
                                                          self.output_directory,
                                                          spm_mat_file)

            # Counting the number of spmT and con files to create
            #if self.multi_reg not in [Undefined, '<undefined>', []]:
            nb_spmT = 0

            for reg_file in self.multi_reg:
                if os.path.splitext(reg_file)[1] != '.txt':
                    nb_spmT += 1

            spmT_files = []
            con_files = []

            for i in range(nb_spmT):
                i += 1
                spmT_files.append(os.path.join(self.output_directory,
                                               'spmT_{:04d}.nii'.format(i)))
                con_files.append(os.path.join(self.output_directory,
                                              'con_{:04d}.nii'.format(i)))
                #spmT_files.append(os.path.join(path, 'spmT_{:04d}.nii'.format(i)))
                #con_files.append(os.path.join(path, 'con_{:04d}.nii'.format(i)))

            if spmT_files:
                self.outputs['spmT_images'] = spmT_files
                self.outputs['con_images'] = con_files

        if self.outputs:
            # FIXME: Quick and dirty. This brick was written quickly and will
            #        need to be reworked to cover all cases. At that time,
            #        inheritance will also need to be reviewed. Currently we
            #        take the first element of the lists, but in the general
            #        case there can be several elements in the lists
            self.inheritance_dict[self.outputs['out_spm_mat_file']] = self.spm_mat_file

            if spmT_files:
                self.inheritance_dict[self.outputs['spmT_images'][0]] = self.spm_mat_file
                self.inheritance_dict[self.outputs['con_images'][0]] = self.spm_mat_file

        # What about spmF images ?
        #return outputs

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """blabla"""

        super(EstimateContrast, self).run_process_mia()

        #self.process.spm_mat_file = os.path.abspath(self.spm_mat_file)
        self.process.spm_mat_file = self.spm_mat_file
        self.process.contrasts = self._get_contrasts(self.session_type)
        self.process.beta_images = self.beta_images
        #self.process.residual_image = os.path.abspath(self.residual_image)
        self.process.residual_image = self.residual_image

        if self.use_derivs is not None:
            self.process.use_derivs = self.use_derivs

        else:

            if self.group_contrast is not None:
                self.process.group_contrast = self.group_contrast

            else:
                self.process.use_derivs = False
        #self.process.run()

        return self.process.run(configuration_dict={})

class EstimateModel(ProcessMIA):
    """
 - EstimateModel (User_processes.stats.spm.stats.EstimateModel) <=> Model estimation (SPM12 name).
*** Estimation of model parameters using classical (ReML - Restricted Maximum Likelihood) or Bayesian algorithms.***
    * Input parameters:
        # spm_mat_file <=> spmmat: The SPM.mat file that contains the design specification (a file object).
            <ex. /home/ArthurBlair/data/raw_data/SPM.mat>
        # estimation_method <=> method: Estimation  procedures for fMRI models (A dictionary with keys which are ‘Classical’ or ‘Bayesian’ or ‘Bayesian2’ and with values which are 1).
                                        - Classical: Restricted Maximum Likelihood (ReML) estimation of first or second level models
                                        - Bayesian: Bayesian estimation of first level models (not yet fully implemented)
                                        - Bayesian2: Bayesian estimation of second level models (not yet fully implemented)
             <ex. {'Classical': 1}>
        # write_residuals <=> write_residuals: Write images of residuals to disk. This is only implemented for classical inference (a boolean)
             <ex. True>
        # flags: Additional arguments (a dictionary with keys which are any value and with values which are any value)
             <ex. {}>
        # version: Version of spm (a string)
             <ex. spm12>
        # tot_reg_num: The total number of estimated regression coefficients (an integer).
             <ex. 8>
    * Output parameters:
        # out_spm_mat_file: The SPM.mat file containing specification of the design and estimated model parameters (a pathlike object or string representing an existing file).
             <ex. /home/ArthurBlair/data/raw_data/SPM.mat>
        # beta_images: Images of estimated regression coefficients beta_000k.img where k indexes the kth regression coefficient (a list of items which are a pathlike object or string representing an existing file).
             <ex. ['/home/ArthurBlair/data/raw_data/beta_0001.nii',
                   '/home/ArthurBlair/data/raw_data/beta_0002.nii',
                   '/home/ArthurBlair/data/raw_data/beta_0003.nii',
                   '/home/ArthurBlair/data/raw_data/beta_0004.nii',
                   '/home/ArthurBlair/data/raw_data/beta_0005.nii',
                   '/home/ArthurBlair/data/raw_data/beta_0006.nii',
                   '/home/ArthurBlair/data/raw_data/beta_0007.nii',
                   '/home/ArthurBlair/data/raw_data/beta_0008.nii']>
        # mask_image: The mask.img image indicating which voxels were included in the analysis (a pathlike object or string representing an existing file).
             <ex. /home/ArthurBlair/data/raw_data/mask.nii>
        # residual_image: The ResMS.img image of the variance of the error (a pathlike object or string representing an existing file).
             <ex. /home/ArthurBlair/data/raw_data/ResMS.nii >
        # residual_images: The individual error Res_000k.img images (a list of items which are a pathlike object or string representing an existing file) where k indexes the kth dynamic (fourth dimensional points of the fuctional). These images are generated only if write_residuals is True.
             <ex. ['/home/ArthurBlair/data/raw_data/Res_0001.nii',
                   '/home/ArthurBlair/data/raw_data/Res_0002.nii',
                   '/home/ArthurBlair/data/raw_data/Res_0003.nii',
                    ...,
                   '/home/ArthurBlair/data/raw_data/Res_0238.nii',
                   '/home/ArthurBlair/data/raw_data/Res_0239.nii',
                   '/home/ArthurBlair/data/raw_data/Res_0240.nii']>
        # RPVimage: The RPV.img image of the estimated resolution elements per voxel (a pathlike object or string representing an existing file).
             <ex. /home/ArthurBlair/data/raw_data/RPV.nii>
        # labels:
             <ex. >
        # SDerror:
             <ex. >
        # ARcoef: Images of the autoregressive (AR) coefficient (a list of items which are a pathlike object or string representing an existing file).
             <ex. >
        # Cbetas:
             <ex. >
        # SDbetas:
             <ex. >

    """
    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(EstimateModel, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['spm', 'nipype']

        # Inputs description
        spm_mat_file_desc = ('The SPM.mat file that contains the design '
                             'specification (a file object)')
        estimation_method_desc = ('A dictionary of either Classical: 1, '
                                  'Bayesian: 1, or Bayesian2: 1')
        write_residuals_desc = 'Write individual residual images (a boolean)'
        flags_desc = ('Additional arguments (a dictionary with keys which are '
                      'any value and with values which are any value)')
        version_desc = 'Version of spm (a string)'
        tot_reg_num_desc = ('The total number of estimated regression'
                            'coefficients (an integer)')

        # Outputs description
        out_spm_mat_file_desc = ('The file containing specification of the '
                                 'design and estimated model parameters (a '
                                 'pathlike object or string representing an '
                                 'existing file)')
        beta_images_desc = ('Images of estimated regression coefficients (a '
                            'list of items which are a pathlike object or '
                            'string representing an existing file)')
        mask_image_desc = ('The image indicating which voxels were included in '
                           'the analysis (a pathlike object or string '
                           'representing an existing file)')
        residual_image_desc = ('The image of the variance of the error (a '
                               'pathlike object or string representing an '
                               'existing file)')
        residual_images_desc = ('The individual error images (a list of items '
                                'which are a pathlike object or string '
                                'representing an existing file)')
        RPVimage_desc = ('The image of the estimated resolution elements per '
                         'voxel (a pathlike object or string representing an '
                         'existing file).')
        labels_desc = 'blabla'
        SDerror_desc = 'blabla'
        ARcoef_desc = ('Images of the autoregressive (AR) coefficient (a list '
                       'of items which are a pathlike object or string '
                       'representing an existing file).')
        Cbetas_desc = 'blabla'
        SDbetas_desc = 'blabla'

        # Inputs traits
        self.add_trait("spm_mat_file",
                       File(output=False,
                            optional=False,
                            desc=spm_mat_file_desc))

        self.add_trait("estimation_method",
                       traits.Dict(
                           traits.Enum("Classical", "Bayesian", "Bayesian2"),
                           traits.Enum(1),
                           usedefault=True,
                           output=False,
                           optional=True,
                           desc=estimation_method_desc))
        self.estimation_method = {'Classical': 1}

        self.add_trait("write_residuals",
                       traits.Bool(usedefault=True,
                                   output=False,
                                   optional=True,
                                   desc=write_residuals_desc))

        self.add_trait("flags",
                       traits.Dict(output=False,
                                   optional=True,
                                   desc=flags_desc))

        self.add_trait("version",
                       traits.String("spm12",
                                     usedefault=True,
                                     output=False,
                                     optional=True,
                                     desc=version_desc))

        self.add_trait("tot_reg_num",
                       traits.Either((traits.Int(), Undefined),
                                     output=False,
                                     optional=True,
                                     desc=tot_reg_num_desc))
        self.tot_reg_num = Undefined

        # Output traits
        self.add_trait("out_spm_mat_file",
                       File(output=True,
                            optional=False,
                            desc=out_spm_mat_file_desc))

        self.add_trait("beta_images",
                       OutputMultiPath(ImageFileSPM(),
                                       output=True,
                                       optional=True,
                                       desc=beta_images_desc))

        self.add_trait("mask_image",
                       ImageFileSPM(output=True,
                                    optional=True,
                                    desc=mask_image_desc))

        self.add_trait("residual_image",
                       ImageFileSPM(output=True,
                                    optional=True,
                                    desc=residual_image_desc))

        self.add_trait("residual_images",
                       OutputMultiPath(ImageFileSPM(),
                                       output=True,
                                       optional=True,
                                       desc=residual_images_desc))

        self.add_trait("RPVimage",
                       ImageFileSPM(output=True,
                                    optional=True,
                                    desc=RPVimage_desc))

        self.add_trait("labels",
                       ImageFileSPM(output=True,
                                    optional=True,
                                    desc=labels_desc))

        self.add_trait("SDerror",
                       OutputMultiPath(ImageFileSPM(),
                                       output=True,
                                       optional=True,
                                       desc=SDerror_desc))

        self.add_trait("ARcoef",
                       OutputMultiPath(ImageFileSPM(),
                                       output=True,
                                       optional=True,
                                       desc=ARcoef_desc))

        self.add_trait("Cbetas",
                       OutputMultiPath(ImageFileSPM(),
                                       output=True,
                                       optional=True,
                                       desc=Cbetas_desc))

        self.add_trait("SDbetas",
                       OutputMultiPath(ImageFileSPM(),
                                       output=True,
                                       optional=True,
                                       desc=SDbetas_desc))

        self.init_default_traits()
        self.init_process('nipype.interfaces.spm.EstimateModel')

    def _get_dbFieldValue(self, document, field):
        """Return, for a document, the field value from the database.

        :param document: the absolute path of the document.
        :param : the field name.
        :returns: the value of the field for the document in the database.
        """
        complete_path = document
        file_position = (complete_path.find(self.project.getName()) + len(self.project.getName()) + 1)
        database_filename = complete_path[file_position:]
        return self.project.session.get_value(COLLECTION_CURRENT,
                                              database_filename,
                                              field)

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dict),
        if defined here. In order not to include an output in the database,
        this output must be a value of the optional key 'notInDb' of the
        self.outputs dictionary. To work properly this method must return
        self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(EstimateModel, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        #if self.outputs:
        #    self.outputs = {}

        #if self.inheritance_dict:
        #    self.inheritance_dict = {}

        # Old version calling Nipype's method
        '''
        process = spm.EstimateModel()
        if not self.spm_mat_file:
            return {}, {}
        else:
            process.inputs.spm_mat_file = self.spm_mat_file
        if not self.estimation_method:
            return {}, {}
        else:
            process.inputs.estimation_method = self.estimation_method

        process.inputs.write_residuals = self.write_residuals
        process.inputs.flags = self.flags

        outputs = process._list_outputs()
        outputs["out_spm_mat_file"] = outputs.pop("spm_mat_file")
        '''

        # Own list_outputs to avoid reading in the SPM.mat file
        if ( (self.spm_mat_file) and
              (self.spm_mat_file not in ['<undefined>', Undefined]) and
              (self.estimation_method) and
              (self.estimation_method not in ['<undefined>', Undefined]) ):
            im_form = "nii" if "12" in self.version else "img"

            # The management of self.process.output_directory could be delegated
            # to the populse_mia.user_interface.pipeline_manager.process_mia
            # module. We can't do it at the moment because the
            # sync_process_output_traits() of the capsul/process/nipype_process
            # module raises an exception in nipype if the mandatory parameters
            # are not yet defined!
            if self.output_directory:
                self.process.output_directory = self.output_directory

            else:
                print('No output_directory was found...!\n')

            _, spm_mat_file = os.path.split(self.spm_mat_file)
            self.outputs['out_spm_mat_file'] = os.path.join(self.output_directory,
                                                            spm_mat_file)
            # Detecting the number of beta files to create
            betas = []
            # Those lines are false, we cannot read in the
            # self.spm_mat_file because it is not created yet! (1)
            """
            try:
                spm_matFile = scipy.io.loadmat(self.spm_mat_file)
                nb_reg = spm_matFile['SPM'][0,0]['Sess'][0,0]['col'].shape[1]
                print('\n EstimateModel brick: {} regressors found ...\n'.format(nb_reg))
                nb_reg += 1  # Adding the constant value

            except Exception as e:
                print('\nEstimateModel brick exception; {0}: {1}\n'.format(e.__class__, e))
            """
            """
            if (  (not nb_reg) and (self.multi_reg) and
              (self.multi_reg not in ['<undefined>', Undefined])  ):

                for reg_file in self.multi_reg:
                    # Those lines are false, we cannot read in all files in
                    # self.multi_reg because they are not created yet! (2)
                    '''
                    if os.path.splitext(reg_file)[1] == '.txt':

                        with open(reg_file) as f:
                            first_line = f.readline()
                            print('FIRST LINE', first_line)
                            nb_reg += len(first_line.split()) + 1  # Adding the constant value

                    # and for other reg_file?
                    '''

                    # The current workaround to try to detect the number of
                    # beta files to be created!
                    # TODO: find a better way
                    if os.path.basename(reg_file)[0:3] == 'rp_':
                        nb_reg += 6

                    if os.path.basename(reg_file)[0:22] == 'regressor_physio_EtCO2':
                        nb_reg += 1

                if nb_reg:
                    nb_reg += 1 # Adding the constant value
                """

            if self.tot_reg_num in ['<undefined>', Undefined]: #and
                #self._get_dbFieldValue(self.spm_mat_file, 'Regress num') is not None):
                self.tot_reg_num = self._get_dbFieldValue(self.spm_mat_file, 'Regress num')

            # Bayesian and Bayesian2 are not yet fully implemented
            if (('Bayesian' in self.estimation_method.keys()) or
                  ('Bayesian2' in self.estimation_method.keys())):
                self.outputs['labels'] = os.path.join(self.output_directory, "labels.{}".format(im_form))
                '''
                outputs['SDerror'] = glob(os.path.join(path, 'Sess*_SDerror*'))
                outputs['ARcoef'] = glob(os.path.join(path, 'Sess*_AR_*'))
                '''

            if 'Classical' in self.estimation_method.keys():
                self.outputs['residual_image'] = os.path.join(self.output_directory, "ResMS.{}".format(im_form))
                self.outputs['RPVimage'] = os.path.join(self.output_directory, "RPV.{}".format(im_form))
                self.outputs['mask_image'] = os.path.join(self.output_directory, "mask.{}".format(im_form))

                if self.write_residuals:
                    nb_dyn = self._get_dbFieldValue(self.spm_mat_file, 'Dynamic Number')

                    if nb_dyn:
                        nb_dyn = sum(nb_dyn)
                        ress = []

                        for i in range(nb_dyn):
                            i += 1
                            ress.append('Res_{:04d}.{}'.format(i, im_form))

                        if ress:
                            self.outputs['residual_images'] = [os.path.join(self.output_directory, res) for res in ress]

                    else:
                        print('- The number of dynamics could not be determined '
                              'automatically. It is not possible to safely '
                              'create the residual_images output parameter ...')

                for i in range(int(0 if self.tot_reg_num in ['<undefined>', Undefined, None] else self.tot_reg_num)):
                    i += 1
                    betas.append('beta_{:04d}.{}'.format(i, im_form))

                if betas:
                    self.outputs['beta_images'] = [os.path.join(self.output_directory, beta) for beta in betas]

                else:
                    print('- No beta image were found to be added to the database ...')

            if ((self.tot_reg_num in ['<undefined>', Undefined, None]) or
                                         (self.write_residuals and not nb_dyn)):
                self.outputs = {}

        if self.outputs:

            for key, value in self.outputs.items():

                if key == "out_spm_mat_file":
                    self.inheritance_dict[value] = self.spm_mat_file

                if key == "beta_images":

                    for fullname in value:
                        self.inheritance_dict[fullname] = self.spm_mat_file

                if key == "mask_image":
                    self.inheritance_dict[value] = self.spm_mat_file

                if key == "residual_image":
                    self.inheritance_dict[value] = self.spm_mat_file

                if key == "residual_images":

                    for fullname in value:
                        self.inheritance_dict[fullname] = self.spm_mat_file

                if key == "RPVimage":
                    self.inheritance_dict[value] = self.spm_mat_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(EstimateModel, self).run_process_mia()
        self.process.spm_mat_file = os.path.abspath(self.spm_mat_file)
        self.process.estimation_method = self.estimation_method
        self.process.write_residuals = self.write_residuals
        self.process.flags = self.flags

        # Removing the image files to avoid a bug
        # I do not observe the bug now, i comment the following lines:
        # TODO: In fact, self.outputs is == {} at this point(see issue #272).
        #       If necessary, we can use the dict4runtine object?
        #
        # for key, value in self.outputs.items():

        #   if key not in ["out_spm_mat_file"]:

        #        if value not in ["<undefined>", Undefined]:

        #            if type(value) in [list, traits.TraitListObject, traits.List]:

        #                for element in value:

        #                    if os.path.isfile(element):
        #                        os.remove(element)

        #           else:

        #                if os.path.isfile(value):
        #                    os.remove(value)

        #self.process.run()
        return self.process.run(configuration_dict={})

class Level1Design(ProcessMIA):
    """
    *Specification of the GLM design matrix, fMRI data files and filtering*

    Please, see the complete documentation for the Level1Design brick in the
    https://gricad-gitlab.univ-grenoble-alpes.fr/condamie/ec_dev web site

    """

    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Level1Design, self).__init__()

        # Third party pieces of software required for the execution of the brick
        self.requirement = ['spm', 'nipype']

        # Inputs description
        timing_units_desc = "One of 'scans' or 'secs'"
        interscan_interval_desc = 'TR in secs (a float)'
        microtime_resolution_desc = ('The number of time-bins per scan used '
                                     'to build regressors (an integer)')
        microtime_onset_desc = ('The reference time-bin at which the '
                                'regressors are resampled to coincide with '
                                'data acquisition (an integer)')
        sess_scans_desc = ('The fMRI scan for each session (a list of items '
                           'which are a pathlike object or string representing '
                           'an existing file)')
        sess_cond_names_desc = ('The name of each condition (list of items '
                                'which are a list of items which are a '
                                'string)')
        sess_cond_onsets_desc = ('The onset times (in seconds or in scans) of '
                                 'the epochs or events within each condition '
                                 '(a list of items which are a list of items '
                                 'which are a list of items which are a '
                                 'float)')
        sess_cond_durations_desc = ('The duration times (in seconds or in '
                                    'scans) of the epochs within each '
                                    'condition (a list of items which are a '
                                    'list of items which are a list of items '
                                    'which are a float)')
        sess_cond_tmod_desc = ('Allows for the characterisation of linear or '
                               'nonlinear time effects (a list of items which '
                               'are a list of items which are 0 or 1 or 2 or 3 '
                               'or 4 or 5 or 6)')
        sess_cond_pmod_names_desc = (
            'The name of the parametric modulation (a '
            'list of items which are a list of items '
            'which are a list of items which are a '
            'string)')
        sess_cond_pmod_values_desc = ('The values used for the parametric '
                                      'modulation, one for each occurence of '
                                      'the event (a list of items which are a '
                                      'list of items which are a list of items '
                                      'which are a list of items which are a '
                                      'float)')
        sess_cond_pmod_polys_desc = ('The polynomial expansion used for the '
                                     'parametric modulation (a list of items '
                                     'which are a list of items which are a '
                                     'list of items which are 1 or 2 or 3 or 4 '
                                     'or 5 or 6)')
        sess_cond_orth_desc = ('Orthogonalise regressors within trial types '
                               '(a list of items which are a list of items '
                               'which are 0 or 1)')
        sess_multi_desc = ('A *.mat file containing details of the multiple '
                           'experimental conditions for each session (a list '
                           'of items which are a filename)')
        sess_regress_desc = ("Additional columns included in the design "
                             "matrix, which may model effects that would not "
                             "be convolved with the haemodynamic response (a "
                             "list of items which are a list of items which "
                             "are a dictionary with keys which are 'name' or "
                             "'val' and with values which are a string or a "
                             "list of float)")
        sess_multi_reg_desc = ('The .mat/.txt file(s) containing details of '
                               'multiple regressors (a list of items which are '
                               'a filename)')
        sess_hpf_desc = 'High-pass filter (a list of items which are a float)'
        factor_info_desc = ("A list of items which are a dictionary for each "
                            "factor with keys which are 'name' or 'levels' "
                            "and with values which are a string (name of the "
                            "factor) or an integer (number of levels for the "
                            "factor)")
        bases_desc = ("To define basic functions for modeling hemodynamic "
                      "response (a 'none' string or a dictionary with keys "
                      "which are 'hrf' or 'fourier' or 'fourier_han' or "
                      "'gamma' or 'fir' and with values which are a dictionary "
                      "with keys which are 'derivs' or 'length' or 'order' and "
                      "with values which are a list or a float or an integer)")
        volterra_expansion_order_desc = ('One of 1 or 2 (1: do not model '
                                         'interactions, 2: model interactions)')
        global_intensity_normalization_desc = (
            "Global intensity normalisation "
            "with scaling or not (one of "
            "'none' or 'scaling')")
        mask_threshold_desc = ('Masking threshold, defined as proportion of '
                               'globals (a float)')
        mask_image_desc = ('Image for explicitly masking the analysis (a '
                           'pathlike object or string representing a file)')
        model_serial_correlations_desc = ('one of AR(1), or FAST or none '
                                          '(AR(1): autoregressive model, '
                                          'FAST: available in SPM12, '
                                          'none: serial correlation is '
                                          'ignored)')

        # Outputs description
        spm_mat_file_desc = ('SPM.mat file (a pathlike object or string '
                             'representing a file')

        # Inputs traits
        self.add_trait("timing_units",
                       traits.Enum('scans',
                                   'secs',
                                   usedefault=True,
                                   output=False,
                                   optional=True,
                                   desc=timing_units_desc))

        # Plan to retrieve this parameter automatically from the database?
        self.add_trait("interscan_interval",
                       traits.Either(traits.Float(), Undefined,
                                     usedefault=True,
                                     output=False,
                                     optional=False,
                                     desc=interscan_interval_desc))
        self.interscan_interval = Undefined

        self.add_trait("microtime_resolution",
                       traits.Int(16,
                                  usedefault=True,
                                  output=False,
                                  optional=True,
                                  desc=microtime_resolution_desc))

        # In study without slice-timing correction, as cevastoc32, it should be 8?
        self.add_trait("microtime_onset",
                       traits.Int(8,
                                  usedefault=True,
                                  output=False,
                                  optional=True,
                                  desc=microtime_onset_desc))
        #TODO: mictotime onset (fmri_spec.timing.fmri_t0) = 1 in Amigo

        self.add_trait("sess_scans",
                       InputMultiPath(File(exists=True),
                                      usedefault=True,
                                      output=False,
                                      desc=sess_scans_desc))

        self.add_trait("sess_cond_names",
                       traits.List(traits.Either(traits.List(traits.String()),
                                                 None),
                                   value=[[]],
                                   output=False,
                                   optional=True,
                                   desc=sess_cond_names_desc))

        self.add_trait("sess_cond_onsets",
                       traits.List(traits.Either(traits.List(traits.List(
                           traits.Float())),
                           None),
                           value=[[[]]],
                           output=False,
                           optional=True,
                           desc=sess_cond_onsets_desc))

        self.add_trait("sess_cond_durations",
                       traits.List(traits.Either(traits.List(traits.List(
                           traits.Float())),
                           None),
                           value=[[[]]],
                           output=False,
                           optional=True,
                           desc=sess_cond_durations_desc))

        self.add_trait("sess_cond_tmod",
                       traits.List(traits.Either(traits.List(traits.Enum(0, 1,
                                                                         2, 3,
                                                                         4, 5,
                                                                         6)),
                                                 None),
                                   value=[[0]],
                                   output=False,
                                   optional=True,
                                   desc=sess_cond_tmod_desc))

        self.add_trait("sess_cond_pmod_names",
                       traits.List(traits.Either(
                           traits.List(traits.Either(
                               traits.List(
                                   traits.String()),
                               None)),
                           None),
                           value=[[[]]],
                           output=False,
                           optional=True,
                           desc=sess_cond_pmod_names_desc))

        self.add_trait("sess_cond_pmod_values",
                       traits.List(traits.Either(
                           traits.List(traits.Either(
                               traits.List(
                                   traits.List(
                                       traits.Float())),
                               None)),
                           None),
                           value=[[[[]]]],
                           output=False,
                           optional=True,
                           desc=sess_cond_pmod_values_desc))

        self.add_trait("sess_cond_pmod_polys",
                       traits.List(traits.Either(
                           traits.List(traits.Either(
                               traits.List(
                                   traits.Enum(1, 2,
                                               3, 4,
                                               5, 6)),
                               None)),
                           None),
                           value=[[[]]],
                           output=False,
                           optional=True,
                           desc=sess_cond_pmod_polys_desc))

        self.add_trait("sess_cond_orth",
                       traits.List(
                           traits.Either(traits.List(traits.Enum(0, 1)),
                                         None),
                           value=[[]],
                           output=False,
                           optional=True,
                           desc=sess_cond_orth_desc))

        self.add_trait("sess_multi",
                       traits.List(traits.File(),
                                   value=[],
                                   output=False,
                                   optional=True,
                                   desc=sess_multi_desc))

        self.add_trait("sess_regress",
                       traits.List(traits.List(traits.Dict(
                           traits.Enum("name", "val"),
                           traits.Union(
                               traits.Str,
                               traits.List(traits.Float())))),
                           value=[[]],
                           output=False,
                           optional=True,
                           desc=sess_regress_desc))

        self.add_trait("sess_multi_reg",
                       InputMultiPath(traits.Either(traits.List(traits.File()),
                                                    None),
                                      value=[[]],
                                      output=False,
                                      optional=True,
                                      desc=sess_multi_reg_desc))

        self.add_trait("sess_hpf",
                       traits.List(traits.Float(),
                                   value=[427.2],
                                   usedefault=True,
                                   output=False,
                                   optional=True,
                                   desc=sess_hpf_desc))
        #TODO: 427.2 corresponds to the value used in Amigo
        # Duration * TR * 3.56 = 427.2 if Duration == 40; TR = 3s; why 3.56 ?
        # I was expecting rather to:
        # (time between first block start - second block start) * TR *2 = 80 *3 *2 = 480
        # Can we code an automatic recovery procedure for the hpf parameter ?
        # (in this case, the user would have to declare additional tags in the
        # database, like the block duration !)

        self.add_trait("factor_info",
                       traits.List(traits.Dict(traits.Enum("name", "levels"),
                                               traits.Either(traits.Str,
                                                             traits.Int)),
                                   usedefault=True,
                                   output=False,
                                   optional=True,
                                   desc=factor_info_desc))

        # traits.Union available from traits 6.1
        self.add_trait("bases",
                       traits.Union(traits.Dict(
                           traits.Enum("hrf",
                                       "fourier",
                                       "fourier_han",
                                       "gamma",
                                       "fir"),
                           traits.Dict(
                               traits.Enum("derivs",
                                           "length",
                                           "order"),
                               traits.Union(traits.Enum([0, 0],
                                                        [1, 0],
                                                        [1, 1]),
                                            traits.Int,
                                            traits.Float))),
                           traits.Enum(['none']),
                           usedefault=True,
                           output=False,
                           optional=True,
                           desc=bases_desc))
        self.bases = {"hrf": {"derivs": [0, 0]}}

        self.add_trait("volterra_expansion_order",
                       traits.Enum(1,
                                   2,
                                   usedefault=True,
                                   output=False,
                                   optional=True,
                                   desc=volterra_expansion_order_desc))

        self.add_trait("global_intensity_normalization",
                       traits.Enum('none',
                                   'scaling',
                                   usedefault=True,
                                   output=False,
                                   optional=True,
                                   desc=global_intensity_normalization_desc))

        self.add_trait("mask_threshold",
                       traits.Float(0.8,
                                    usedefault=True,
                                    output=False,
                                    optional=True,
                                    desc=mask_threshold_desc))

        self.add_trait("mask_image",
                       ImageFileSPM(output=False,
                                    optional=True,
                                    desc=mask_image_desc))

        self.add_trait("model_serial_correlations",
                       traits.Enum('AR(1)',
                                   'FAST',
                                   'none',
                                   usedefault=True,
                                   output=False,
                                   optional=True,
                                   desc=model_serial_correlations_desc))

        # Output traits
        self.add_trait("spm_mat_file",
                       File(output=True,
                            desc=spm_mat_file_desc))

        # Special parameter used as a messenger for the run_process_mia method
        self.add_trait("dict4runtime",
                       traits.Dict(output=False,
                                   optional=True,
                                   userlevel=1))

        self.init_default_traits()
        self.init_process('nipype.interfaces.spm.Level1Design')

    def _get_conditions(self, idx_session, idx_cond):
        """Generate the condition dictionary.

        One condition dictionary is generated for each condition.
        The condition dictionary has the following keys:
        - name (a string)
        - onset (a list of float)
        - duration (a list of float)
        - tmod (an integer between 0 and 6)
        - pmod a dictionary
        - orth (0 or 1)

        :param idx_session: session index.
        :param idx_cond: condition index
        :returns: conditions dictionary
        """
        cond = dict()
        cond["name"] = self.sess_cond_names[idx_session][idx_cond]
        onset = []

        for (i,
             value) in enumerate(self.sess_cond_onsets[idx_session][idx_cond]):
            onset.insert(i, np.float64(value))

        cond["onset"] = onset

        if ((len(self.sess_cond_durations[idx_session][idx_cond]) ==
             len(self.sess_cond_onsets[idx_session][idx_cond])) or
                (len(self.sess_cond_durations[idx_session][idx_cond]) == 1)):
            duration = []

            for (i,
                 value) in enumerate(self.sess_cond_durations[
                                         idx_session][idx_cond]):
                duration.insert(i, np.float64(value))
            cond["duration"] = duration

        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle("User_processes_ECdev - Level1Design Error!")
            msg.setText("Warning: The number of values in the "
                        "sess_cond_durations parameter does not correspond to "
                        "the number of values in the sess_cond_onsets "
                        "parameter, for the session '{0}' and the "
                        "condition '{1}'! Pease, check your "
                        "settings ...".format(idx_session, cond["name"]))
            msg.setStandardButtons(QMessageBox.Close)
            msg.buttonClicked.connect(msg.close)
            msg.exec()
            return False

        try:
            cond["tmod"] = self.sess_cond_tmod[idx_session][idx_cond]

        except IndexError:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle("User_processes_ECdev - Level1Design Error!")
            msg.setText("Warning: It seems that not all sess_cond_tmod "
                        "parameter values are defined for all conditions. "
                        "Each condition must have a value for the "
                        "sess_cond_tmod parameter. Please, check the "
                        "sess_cond_tmod parameter and try again to initialize "
                        "the pipeline ...!")
            msg.setStandardButtons(QMessageBox.Close)
            msg.buttonClicked.connect(msg.close)
            msg.exec()
            return False

        if ((self.sess_cond_pmod_names[idx_session]) and
                (self.sess_cond_pmod_names[idx_session][idx_cond])):
            pmods = []

            for (i, val) in enumerate(
                    self.sess_cond_pmod_names[idx_session][idx_cond]):
                pmod = dict()
                pmod["name"] = val
                param = []

                for (j, val) in enumerate(
                        self.sess_cond_pmod_values[idx_session][idx_cond][i]):
                    param.insert(j, np.float64(val))

                pmod["param"] = param
                pmod["poly"] = self.sess_cond_pmod_polys[idx_session][
                    idx_cond][i]
                pmods.append(pmod)

            cond["pmod"] = pmods

        try:
            cond["orth"] = self.sess_cond_orth[idx_session][idx_cond]

        except IndexError:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setWindowTitle("User_processes_ECdev - Level1Design Error!")
            msg.setText("Warning: It seems that not all sess_cond_orth "
                        "parameter values are defined for all conditions. "
                        "Each condition must have a value for the "
                        "sess_cond_orth parameter. Please, check the "
                        "sess_cond_orth parameter and try again to initialize "
                        "the pipeline ...!")
            msg.setStandardButtons(QMessageBox.Close)
            msg.buttonClicked.connect(msg.close)
            msg.exec()
            return False

        return cond

    def _get_session_info(self, idx_session):
        """Generate the session_info parameter for each session.

        The 'session_info' dictonary is equivalent to the 'sess'
        structure array spm parameter, for each session. It therefore includes
        the 'scans', 'cond', 'multi', 'regress', 'multi_reg' and 'hpf' keys,
        like for the 'sess' parameter of spm.

        :param idx_session: session index.
        :returns: session_info, beta_sess parameters and condition number
        """
        session_info = dict()  # The information for this session
        beta_sess = 0  # The regressors number for this session
        cond_nb = 0  # The condition number for this session

        # scans
        session_info['scans'] = self.sess_scans[idx_session]

        # cond
        if self.sess_cond_names[idx_session]:
            cond = []

            # idx_cond: condition index
            for idx_cond in range(len(self.sess_cond_names[idx_session])):
                condition = self._get_conditions(idx_session, idx_cond)
                cond.append(condition)
                beta_sess += 1
                cond_nb += 1

                if 'pmod' in condition.keys():

                    for i in condition['pmod']:
                        beta_sess += i['poly']

                if 'tmod' in condition.keys():
                    beta_sess += condition['tmod']

            session_info['cond'] = cond

        # multi
        if self.sess_multi and self.sess_multi[idx_session]:
            session_info['multi'] = {self.sess_multi[idx_session]}
            mat = scipy.io.loadmat(self.sess_multi[idx_session],
                                   squeeze_me=True)

            if 'names' in mat:
                beta_sess += mat['names'].shape[0]
                cond_nb += mat['names'].shape[0]

            if 'pmod' in mat:

                for i in range(len(mat['pmod'])):

                    if isinstance(mat['pmod'][i][2], int):
                        beta_sess += mat['pmod'][i][2]

                    elif isinstance(mat['pmod'][i][2], np.ndarray):

                        for j in range(len(mat['pmod'][i][2])):
                            beta_sess += mat['pmod'][i][2][j]

            if 'tmod' in mat:

                for i in range(len(mat['tmod'])):
                    beta_sess += mat['tmod'][i]

        # Don't understand what's the point of multiplier attribute.
        # Moreover, it does not seem to be used afterwards!
        # multiplier = 1

        # if 'hrf' in self.bases.keys():
        #    if 'deriv' in self.bases['hrf'].keys():
        #        multiplier = sum(self.bases['hrf']['deriv']) + 1

        # regress
        if ((self.sess_regress not in ['<undefined>', Undefined]) and
                (not all([not elem for elem in self.sess_regress])) and
                (self.sess_regress[idx_session])):

            for regressor in self.sess_regress[idx_session]:
                val = []

                for j, value in enumerate(regressor['val']):
                    val.insert(j, np.float64(value))

                regressor['val'] = val
                beta_sess += 1

            session_info['regress'] = self.sess_regress[idx_session]

        # multi_reg
        if ((self.sess_multi_reg not in ['<undefined>', Undefined]) and
                (not all([not elem for elem in self.sess_multi_reg])) and
                (self.sess_multi_reg[idx_session])):
            session_info['multi_reg'] = []

            for reg_file in self.sess_multi_reg[idx_session]:

                if os.path.basename(reg_file)[0:3] == 'rp_':
                    beta_sess += 6

                if os.path.splitext(reg_file)[1] == '.mat':
                    mat = scipy.io.loadmat(reg_file)
                    beta_sess += mat['R'].shape[1]

                reg_to_add = [{reg_file}]
                session_info['multi_reg'].append(reg_to_add)

        # hpf
        if self.sess_hpf[idx_session]:
            session_info['hpf'] = self.sess_hpf[idx_session]

        beta_sess += 1
        return session_info, beta_sess, cond_nb

    def _get_dbFieldValue(self, document, field):
        """Return, for a document, the field value from the database.

        :param document: the absolute path of the document.
        :param : the field name.
        :returns: the value of the field for the document in the database.
        """
        file_position = (document.find(self.project.getName()) +
                         len(self.project.getName()) +
                         1)
        database_filename = document[file_position:]

        return self.project.session.get_value(COLLECTION_CURRENT,
                                              database_filename,
                                              field)

    def _set_dbFieldValue(self, document, tag_to_add):
        """blabla
        """
        file_position = (document.find(self.project.getName()) +
                         len(self.project.getName()) +
                         1)
        database_filename = document[file_position:]

        field_names = self.project.session.get_fields_names(COLLECTION_CURRENT)

        if tag_to_add["name"] not in field_names:
            (self.project.session.add_field)(
                COLLECTION_CURRENT,
                tag_to_add["name"],
                tag_to_add["field_type"],
                tag_to_add["description"],
                tag_to_add["visibility"],
                tag_to_add["origin"],
                tag_to_add["unit"],
                tag_to_add["default_value"],
            )

        if tag_to_add["name"] not in (
                self.project.session.get_fields_names)(COLLECTION_INITIAL):
            (self.project.session.add_field)(
                COLLECTION_INITIAL,
                tag_to_add["name"],
                tag_to_add["field_type"],
                tag_to_add["description"],
                tag_to_add["visibility"],
                tag_to_add["origin"],
                tag_to_add["unit"],
                tag_to_add["default_value"],
            )

        if self.project.session.get_document(COLLECTION_CURRENT, database_filename):
            print("Path {0} already in database.".format(database_filename))

        else:
            self.project.session.add_document(COLLECTION_CURRENT, database_filename)
            self.project.session.add_document(COLLECTION_INITIAL, database_filename)

        self.project.session.set_values(COLLECTION_CURRENT, database_filename, {tag_to_add["name"]: tag_to_add['value']})
        self.project.session.set_values(COLLECTION_INITIAL, database_filename, {tag_to_add["name"]: tag_to_add['value']})

    def list_outputs(self, is_plugged=None):
        """Dedicated to the initialisation step of the brick.

        The main objective of this method is to produce the outputs of the
        bricks (self.outputs) and the associated tags (self.inheritance_dict),
        The optional self.inheritance_dict can have two structures. On the one
        hand it can be a dictionary whose keys are the documents to inherit
        metadata and the values the documents used for the inheritance of
        these metadata. On the other hand, it can be a dictionary whose keys
        are the documents to inherit metadata and the values are dictionary
        with two keys; parent (for usual inheritance) or own_tags (to add a
        new tag or modify an existing one). In order not to include an output
        in the database, the name of the plug related to this output must be
        an element of the list corresponding to the value of the optional key
        "notInDb" of the self.outputs dictionary. To work properly this method
        must return self.make_initResult() object.

        :param is_plugged: the state, linked or not, of the plugs.
        :returns: a dictionary with requirement, outputs and inheritance_dict.
        """
        # Using the inheritance to ProcessMIA class, list_outputs method
        super(Level1Design, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if ((self.sess_scans) and
                (self.sess_scans not in ['<undefined>', Undefined]) and
                (self.sess_scans[0] not in ['<undefined>', Undefined])):
            # The management of self.process.output_directory could be delegated
            # to the populse_mia.user_interface.pipeline_manager.process_mia
            # module. We can't do it at the moment because the
            # sync_process_output_traits() of the capsul/process/nipype_process
            # module raises an exception in nipype if the mandatory parameters
            # are not yet defined!
            if self.output_directory:
                self.process.output_directory = self.output_directory

            else:
                print('No output_directory was found...!\n')

            self.outputs['spm_mat_file'] = os.path.join(self.output_directory,
                                                        'SPM.mat')

            sessions = []  # The total session_info parameter for nipype
            beta = 0  # The total number of regressors in the GLM design matrix
            cond_totNb = 0  # The total number of condition in cond and multi

            # idx_session: session index
            for idx_session in range(len(self.sess_scans)):
                # session_info: All parameters by session (Scans, Conditions,
                # Multiple conditions, Regressors, Multiple regressors,
                # High-pass filter)
                # beta_sess: Number of regressors by session
                session_info, beta_sess, cond_nb = self._get_session_info(
                    idx_session)
                sessions.append(session_info)
                beta += beta_sess
                cond_totNb += cond_nb

            if self.volterra_expansion_order == 2:
                beta += int((cond_totNb + 1) * cond_totNb / 2)

            #self.sessions = sessions  # The session_info for nipype
            self.dict4runtime['sessions'] = sessions

            ## Some tests to check that the definition of the parameters for
            ## the self.sessions went well:
            # - If self.sessions have no key in ['cond', 'multi', 'regress',
            #   'multi_reg'], or a condition is False, we can think there is
            #   an initialisation issue ...
            check = True  # Flag to stop the check processes
            init_res = []  # Initialisation result; True: Ok, False: Fail

            #print('n self.sessions: ', self.sessions)
            print('n sessions: ', sessions)  # To Remove after debug

            for i in sessions:

                for j in ['cond', 'multi', 'regress', 'multi_reg']:

                    if j in i.keys(): init_res.append(True)

                if i.get('cond') is not None and i.get('cond') is False:
                    init_res.append(False)

                if (not True in init_res or False in init_res) and (check):
                    check = False
                    self.outputs = {}
                    print('\nThere seems to be a problem in the definition of '
                          'session parameters. The initialisation failed, '
                          'please, check your settings  ...')

            # - If sess_multi_reg is plugged and sessions have no multi_reg
            #   key, we can think there is an initialisation issue ...

            #print('n self.sess_multi_reg: ', self.sess_multi_reg) # To Remove after debug

            #if is_plugged['sess_multi_reg'] and check:
            if self.sess_multi_reg != [[]] and check:
                init_res = []
                [init_res.append(True) for i in sessions
                 if 'multi_reg' in i.keys()]

                if not any(init_res):
                    self.outputs = {}
                    print('\nThe sess_multi_reg plug is linked to a node. '
                          'However, no parameters for multi_reg have been '
                          'defined in the nipype session_info during '
                          'initialisation. This leads to an initialisation '
                          'failure. Please, unplug the sess_multi_reg plug if '
                          'not needed or check your settings ...')

            # - If sess_multi is plugged and sessions have no multi key,
            #   we can think there is an initialisation issue ...

            print('n self.sess_multi: ', self.sess_multi) # To Remove after debug

            #if is_plugged['sess_multi'] and check:
            if self.sess_multi != [] and check:
                init_res = []
                [init_res.append(True) for i in sessions
                 if 'multi' in i.keys()]

                if not any(init_res):
                    self.outputs = {}
                    print(
                        '\nThe sess_multi plug is linked to a node. However, '
                        'no parameters for sess_multi have been defined in '
                        'the nipype session_info during initialisation. This '
                        'leads to an initialisation failure. Please, unplug '
                        'the sess_multi plug if not needed or check your '
                        'settings ...')

        if self.outputs:
            self.inheritance_dict[self.outputs['spm_mat_file']] = dict()
            # Currently, spm_mat_file will only inherit the first scan if there
            # are several scans in self.scans. This would require some thought
            # for this particular case!
            self.inheritance_dict[self.outputs['spm_mat_file']]['parent'] = (
                self.sess_scans[0])
            self.inheritance_dict[self.outputs['spm_mat_file']][
                'own_tags'] = []
            tag_to_add = dict()
            tag_to_add['name'] = 'Regress num'
            tag_to_add['field_type'] = FIELD_TYPE_INTEGER
            tag_to_add['description'] = 'Total number of regressors'
            tag_to_add['visibility'] = True
            tag_to_add['origin'] = TAG_ORIGIN_USER
            tag_to_add['unit'] = None
            tag_to_add['default_value'] = None
            tag_to_add['value'] = beta
            self.inheritance_dict[self.outputs['spm_mat_file']][
                'own_tags'].append(tag_to_add)
            # FIXME: In the latest version of mia, indexing of the database with
            #        particular tags defined in the processes is done only at
            #        the end of the initialisation of the whole pipeline. So we
            #        cannot use the value of these tags in other processes of
            #        the pipeline at the time of initialisation
            #        (see populse_mia #290). Until better I use a quick and
            #        dirty hack with the _set_dbFieldValue() method !
            self._set_dbFieldValue(self.outputs['spm_mat_file'], tag_to_add)
            dyn_num = 0

            for scan in self.sess_scans:
                dimensions = self._get_dbFieldValue(scan,
                                                    'Dataset dimensions (Count, X,Y,Z,T...)')

                if ((dimensions is not None) and
                        (isinstance(dimensions, list)) and
                        (len(dimensions) == 5)):
                    dyn_num += dimensions[
                        4]  # to check with EstimateModel brick; total sum of the dynamics in self.sess_scans -> strange !
                    tag_to_add = dict()
                    tag_to_add['name'] = 'Dynamic Number'
                    tag_to_add['field_type'] = FIELD_TYPE_INTEGER
                    tag_to_add['description'] = ('Total number of dynamics in '
                                                 'the functionals')
                    tag_to_add['visibility'] = True
                    tag_to_add['origin'] = TAG_ORIGIN_USER
                    tag_to_add['unit'] = None
                    tag_to_add['default_value'] = None
                    tag_to_add['value'] = dyn_num
                    self.inheritance_dict[self.outputs['spm_mat_file']
                    ]['own_tags'].append(tag_to_add)
                    # FIXME: In the latest version of mia, indexing of the database with
                    #        particular tags defined in the processes is done only at
                    #        the end of the initialisation of the whole pipeline. So we
                    #        cannot use the value of these tags in other processes of
                    #        the pipeline at the time of initialisation
                    #        (see populse_mia #290). Until better I use a quick and
                    #        dirty hack with the _set_dbFieldValue() method !
                    self._set_dbFieldValue(self.outputs['spm_mat_file'],
                                           tag_to_add)

                else:
                    print('\nWarning! The dynamics number for at least one '
                          'scan could not be found in the database. This can '
                          'cause a running issue in next bricks (for example '
                          'the EstimateModel brick) ... !')

            if ((self.interscan_interval is Undefined) and ('RepetitionTime' in
                                                            self.project.session.get_fields_names(
                                                                COLLECTION_CURRENT))):

                if self._get_dbFieldValue(self.sess_scans[0], 'RepetitionTime') is not None:
                    self.interscan_interval = self._get_dbFieldValue(
                                                    self.sess_scans[0],
                                                    'RepetitionTime')[0] / 1000

                else:
                    self.outputs = {}
                    print('\nThe interscan_interval parameter (repetition time '
                          'in seconds) could not be determined automatically '
                          'for the Level1Design brick. Please set this value '
                          'and launch the calculation again!\n')

            elif self.interscan_interval is Undefined:
                self.outputs = {}
                print('\nThe interscan_interval parameter (repetition time '
                      'in seconds) could not be determined automatically '
                      'for the Level1Design brick. Please set this value '
                      'and launch the calculation again!\n')

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Level1Design, self).run_process_mia()
        # Removing the spm_mat_file to avoid a bug (nipy/nipype Issues #2612)
        #cur_dir = os.getcwd()
        #out_file = os.path.join(cur_dir, 'SPM.mat')

        if self.output_directory:
            out_file = os.path.join(self.output_directory, "SPM.mat")

            if os.path.isfile(out_file):
                os.remove(out_file)

        else:
            print('No output_directory was found...!\n')

        #if os.path.isfile(out_file):
        #    os.remove(out_file)

        self.process.timing_units = self.timing_units
        self.process.interscan_interval = self.interscan_interval
        self.process.microtime_resolution = self.microtime_resolution
        self.process.microtime_onset = self.microtime_onset
        #self.process.session_info = self.sessions
        self.process.session_info = self.dict4runtime['sessions']
        self.process.factor_info = self.factor_info
        self.process.bases = self.bases
        self.process.volterra_expansion_order = self.volterra_expansion_order
        (self.process.
         global_intensity_normalization) = self.global_intensity_normalization
        self.process.mask_threshold = self.mask_threshold

        # Only one mask can be defined in spm. If more than one is given, only
        # the first one will be used for all sessions ...
        if self.mask_image != Undefined:

            if type(self.mask_image) in [list,
                                         traits.TraitListObject,
                                         traits.List]:
                self.process.mask_image = os.path.abspath(self.mask_image[0])

            else:
                self.process.mask_image = os.path.abspath(self.mask_image)

        self.process.model_serial_correlations = self.model_serial_correlations

        #self.process.run()
        return self.process.run(configuration_dict={})

        # Copying the generated SPM.mat file in the data directory
        #if ((self.sess_scans) and
        #        (self.sess_scans not in ['<undefined>', Undefined]) and
        #        (self.sess_scans[0] not in ['<undefined>', Undefined])):
        #    scan_image = os.path.abspath(self.sess_scans[0])
        #    scan_folder, _ = os.path.split(scan_image)

        #copy2(out_file, scan_folder)

        #if os.path.isfile(out_file):
        #    os.remove(out_file)
