# mia_processes import
from mia_processes.process_mia import Process_Mia
from aslpipeline.asl_processes import ASL_processes
from aslpipeline.asl_processes import Anat

# nipype import
from nipype.interfaces.base import File, traits, OutputMultiPath

# populse_mia import
from populse_mia.data_manager.project import (
    COLLECTION_CURRENT, COLLECTION_INITIAL, COLLECTION_BRICK, TAG_CHECKSUM,
    TAG_FILENAME, TAG_BRICKS, BRICK_NAME)
from populse_mia.data_manager.database_mia import (
    TAG_ORIGIN_BUILTIN, TAG_ORIGIN_USER)

# Other import
import os
import nibabel as nib
import numpy as np
import nipype.interfaces.spm as spm
import ast # for tools

class tools:
    def add_new_tag(self, new_tag_name, new_default_value, tag_type,
                      new_tag_description, new_tag_unit):

        # We add the tag and a value for each scan in the Database
        self.project.session.add_field(COLLECTION_CURRENT, new_tag_name,
                                       tag_type, new_tag_description, True,
                                       TAG_ORIGIN_USER, new_tag_unit,
                                       new_default_value)
        self.project.session.add_field(COLLECTION_INITIAL, new_tag_name,
                                       tag_type, new_tag_description, True,
                                       TAG_ORIGIN_USER, new_tag_unit,
                                       new_default_value)

        old_list = ast.literal_eval(new_default_value)
        list_to_return = []
        for old_element in old_list:
            list_to_return.append(int(
                old_element))

        for scan in self.project.session.get_documents(COLLECTION_CURRENT):
            self.project.session.add_value(
                COLLECTION_CURRENT, getattr(scan, TAG_FILENAME),
                new_tag_name, list_to_return)
            self.project.session.add_value(
                COLLECTION_INITIAL, getattr(scan, TAG_FILENAME),
                new_tag_name, list_to_return)


class ASL_Rescaling(Process_Mia, ASL_processes): #PROBLEM: with database // cant change the name, does not work

    def __init__(self):
        super(ASL_Rescaling, self).__init__()
        # Inputs
        self.add_trait("ASL_file", traits.File(output=False, optional = False))
        self.add_trait("method", traits.Str('whitepaper', output = False, optional = True))
        self.add_trait("M0_image_file", traits.File(output = False, optional = True))
        self.add_trait("T1_map_file", traits.File(output=False,optional = True))
        self.add_trait("T1_blood", traits.Float(1.65, output = False, optional = True)) # seconds
        self.add_trait("T2w_blood", traits.Float(0.1, output = False, optional = True)) # seconds
        self.add_trait("rho_blood", traits.Float(0.8, output = False, optional = True))
        self.add_trait("T1_CSF", traits.Float(4.0, output = False, optional = True)) # seconds
        self.add_trait("T2w_CSF", traits.Float(0.1, output = False, optional = True)) # seconds
        self.add_trait("rho_CSF", traits.Float(1.0, output = False, optional = True))
        self.add_trait("rho_brain", traits.Float(1040, output = False, optional = True)) #kg/m^3
        self.add_trait("lambda_", traits.Float(0.0009, output = False, optional = True)) #kg/m^3
        self.add_trait("label_eff", traits.Float(output = False, optional = True))
        self.add_trait("transit_delay", traits.Float(1.5, output = False, optional = True)) #seconds
        self.add_trait("tAcq_slice", traits.Float(0.0,output=False, optional = True))
        self.add_trait("back_supp", traits.Int(output=False, optional = True))
        self.add_trait("TI1", traits.Float(1.8,output = False,optional = True))
        self.add_trait("PLD", traits.Float(2.0,output = False,optional = True))
        self.add_trait("TRm0", traits.Float(10.0,output = False,optional = True))
        self.add_trait("T1t", traits.Float(0.9,output = False,optional = True))
        self.add_trait("ASL_mode", traits.Str(output = False, optional = True))
        self.add_trait("output_prefix", traits.Str('t',output = False,optional = True))
        self.add_trait("output_datatype", traits.Str('int16',output = False,optional = True))

        # Outputs
        self.add_trait("ASL_output", traits.File(output=True))
        self.add_trait("M0_output", traits.File(output=True))

    def list_outputs(self):
        super(ASL_Rescaling, self).list_outputs()
        if not self.ASL_file or self.ASL_file in ["<undefined>", traits.Undefined]:
            return {}, {}
        elif not self.M0_image_file:
            path, filename = os.path.split(self.ASL_file)
            filename1 = self.output_prefix + "ASL_" + filename
            filename2 = self.output_prefix + "M0_" + filename
            self.ASL_output = os.path.join(path, filename1)
            self.M0_output = os.path.join(path, filename2)
            print(self.ASL_output, self.M0_output, self.ASL_file)
            return {"ASL_output": self.ASL_output,"M0_output": self.M0_output},\
                   {self.ASL_output: self.ASL_file, self.M0_output: self.ASL_file}
        else:
            path, filename = os.path.split(self.ASL_file)
            path2, filename2 = os.path.split(self.M0_image_file)
            filename = self.output_prefix + "ASL_" + filename
            filename2 = self.output_prefix + "M0_" + filename2
            self.ASL_output = os.path.join(path, filename)
            self.M0_output = os.path.join(path2, filename2)
            return {"ASL_output": self.ASL_output,"M0_output": self.M0_output},\
                   {self.ASL_output: self.ASL_file, self.M0_output: self.M0_image_file}


    def run_process_mia(self):
        super(ASL_Rescaling, self).run_process_mia()

        # get tags
        ASLPathName = os.path.relpath(self.ASL_file, self.project.folder)
        if self.M0_image_file:
            try:
                M0PathName = os.path.relpath(self.M0_image_file, self.project.folder)
                TEm0 = self.project.session.get_value(COLLECTION_INITIAL, ASLPathName, "EchoTime")[0] * 1e-3
                TE = self.project.session.get_value(COLLECTION_INITIAL, M0PathName, "EchoTime")[0] * 1e-3
            except TypeError:
                print("There is no echo value in the data")
                return
        else:
            try:
                TE = self.project.session.get_value(COLLECTION_INITIAL, ASLPathName, "EchoTime")[0] * 1e-3
            except TypeError:
                print("There is no echo value")
                return

        ScanMode = self.project.session.get_value(COLLECTION_INITIAL, ASLPathName, "ScanMode")

        # image data
        asl_image = nib.load(self.ASL_file)
        asl_image_data = asl_image.get_fdata()

        self.params = {'ASL_file': self.ASL_file,
                       'method': self.method,
                       'M0_image_file': self.M0_image_file,
                       'T1_map_file': self.T1_map_file,
                       'T1_blood': self.T1_blood,
                       'T2w_blood': self.T2w_blood,
                       'rho_blood': self.rho_blood,
                       'T1_CSF': self.T1_CSF,
                       'T2w_CSF': self.T2w_CSF,
                       'rho_CSF': self.rho_CSF,
                       'rho_brain': self.rho_brain,
                       'lambda_': self.lambda_,
                       'label_eff': self.label_eff,
                       'transit_delay': self.transit_delay,
                       'tAcq_slice': self.tAcq_slice,
                       'image':asl_image,
                       'TI1': self.TI1,
                       'PLD':self.PLD,
                       'TRm0':self.TRm0,
                       'T1t': self.T1t,
                       'ASL_mode': self.ASL_mode,
                       'back_supp': self.back_supp,
                       'TE': TE,
                       'scan_mode': ScanMode}

        Processes = ASL_processes()

        # check dimensions
        if len(asl_image_data.shape) <= 3:
            print('Input size:',len(asl_image_data.shape),'is not correct, it should be 5 dimensions.')
            return
        else:
            asl_image_data = Processes.fiveD_to_fourD(asl_image_data)
        SFasl, SFm0_1, SFm0_2end = Processes.Rescaling_process(self.params)

        if not self.M0_image_file:
            m0_image_data = np.array(asl_image_data[:,:,:,0:2]) # split the images
            asl_image_data = np.array(asl_image_data[:,:,:,2:])


        else:
            m0_image = nib.load(self.M0_image_file)
            m0_image_data = m0_image.get_fdata()
            m0_image_data = Processes.fiveD_to_fourD(m0_image_data)

        # get scaled images (ASL and M0)
        SFm0_1 = np.tile(SFm0_1, [1, 1, m0_image_data.shape[2], 1])
        SFm0_2end = np.tile(SFm0_2end,[1,1,m0_image_data.shape[2], m0_image_data.shape[3]-1])
        SFm0 = np.concatenate((SFm0_1, SFm0_2end), axis = 3)
        scaled_m0_image = Processes.Apply_SF(m0_image_data,SFm0)
        scaled_asl_image = Processes.Apply_SF(asl_image_data,SFasl)

        # create de image object
        hdr = asl_image.header.set_data_shape(scaled_asl_image.shape)
        M0_hdr = asl_image.header.set_data_shape(scaled_m0_image.shape)
        asl_image = nib.Nifti1Image(scaled_asl_image,asl_image.affine, hdr)
        m0_image = nib.Nifti1Image(scaled_m0_image,asl_image.affine, hdr)

        # change data type
        asl_image.set_data_dtype(self.output_datatype)
        m0_image.set_data_dtype(self.output_datatype)

        # save image
        nib.save(asl_image, self.ASL_output)
        nib.save(m0_image, self.M0_output)


class ASL_Realign(Process_Mia, ASL_processes, tools):

    def __init__(self):
        super(ASL_Realign, self).__init__()
        # Inputs
        self.add_trait("ASL_file", traits.File(output = False, optional = False))
        self.add_trait("M0_mean_file", traits.File(output = False, optional = False))
        self.add_trait("quality", traits.Float(0.9, output = False, optional = True)) #this value is in the tags, should I take the tag??
        self.add_trait("FWHM", traits.Float(5.0, output = False, optional = True))
        self.add_trait("separation", traits.Float(4.0, output = False, optional = True))
        self.add_trait("interpolation", traits.Int(2, output = False, optional = True))
        self.add_trait("wrapping", traits.List([False, False, False],output = False, optional = True))
        self.add_trait("register_to_mean", traits.Bool(True, output = False, optional = True))
        self.add_trait("weighting", traits.Str('', output = False, optional = True))
        # self.add_trait("save_graphics", traits.Bool(True, output = False, optional = True))
        self.add_trait("motion_threshold", traits.Float(1.0, output = False, optional = True))
        self.add_trait("exclude_fname", traits.Str('excluded_frames.nii', output = False, optional = True)) #nifty?
        self.add_trait("output_prefix", traits.Str('r', output = False, optional = True))
        self.add_trait("write_mask", traits.Bool(True, output = False, optional = True))
        self.add_trait("save_params", traits.Bool(True, output = False, optional = True)) # delete when we know how to handle

        # Outputs
        self.add_trait("ASL_realigned", traits.File(output=True))
        self.add_trait("ASL_mean", traits.File(output=True))
        self.add_trait("graphics_fname", traits.File(output=True))
        self.add_trait("realignment_params", traits.File(output=True, optional = True))

    def list_outputs(self):
        super(ASL_Realign, self).list_outputs()
        if not self.ASL_file or self.ASL_file in ["<undefined>", traits.Undefined]:
            return {}, {}

        if not self.M0_mean_file or self.M0_mean_file in ["<undefined>", traits.Undefined]:
            return {}, {}
        prefix = self.output_prefix
        path, filename = os.path.split(self.ASL_file)
        filename, file_extension = os.path.splitext(filename)
        output_dict = {}
        inheritance_dict = {}
        self.write_which = [0, 0]
        self.realignment_params = ''
        # TODO: Make parameters optional and handle output initialization correctly
        # if self.ASL_realigned != 'test': # and self.ASL_realigned not in ["<undefined>", traits.Undefined]:
        filename1 = prefix + filename + file_extension
        self.write_which[0] = 2
        self.ASL_realigned = os.path.join(path, filename1)
        output_dict['ASL_realigned'] = self.ASL_realigned
        inheritance_dict[self.ASL_realigned] = self.ASL_file
        # else:
        #     print('NOT initializing ASL_realigned parameter. Value of ASL_realigned:')
        #     print(self.ASL_realigned)

        # if self.ASL_mean != 'test': # and self.ASL_mean not in ["<undefined>", traits.Undefined]:
        filename2 = 'mean' + filename + file_extension
        self.write_which[1] = 1
        self.ASL_mean = os.path.join(path, filename2)
        output_dict['ASL_mean'] = self.ASL_mean
        inheritance_dict[self.ASL_mean] = self.ASL_file
        # else:
        #     print('NOT initializing ASL_mean parameter. Value of ASL_mean:')
        #     print(self.ASL_mean)

        # if self.graphics_fname and self.graphics_fname not in ["<undefined>", traits.Undefined]:
        filename3 = prefix + 'plot' + filename + '.png'
        self.save_graphics = True
        self.graphics_fname = os.path.join(path, filename3)
        output_dict['graphics_fname'] = self.graphics_fname
        inheritance_dict[self.graphics_fname] = self.ASL_file
        # else:
        #     print('NOT initializing graphics_fname parameter. Value of graphics_fname:')
        #     print(self.graphics_fname)

        # if self.realignment_params and self.realignment_params not in ["<undefined>", traits.Undefined]:
        filename4 = 'rp_' + filename + '.txt'
        self.realignment_params = os.path.join(path, filename4)
        output_dict['realignment_params'] = self.realignment_params
        inheritance_dict[self.realignment_params] = self.ASL_file
        # else:
        #     print('NOT initializing realignment_params parameter. Value of realignment_params:')
        #     print(self.realignment_params)

        return output_dict, inheritance_dict

    def run_process_mia(self):
        super(ASL_Realign, self).run_process_mia()
        im = nib.load(self.ASL_file)
        img_data = im.get_fdata()

        params = {'ASL_file': self.ASL_file,
                  'quality': self.quality,
                  'FWHM': self.FWHM,
                  'separation': self.separation,
                  'interpolation': self.interpolation,
                  'register_to_mean': self.register_to_mean,
                  'weighting': self.weighting,
                  'save_graphics': self.save_graphics,
                  'motion_threshold': self.motion_threshold,
                  'exclude_fname': self.exclude_fname,
                  'wrapping': self.wrapping,
                  'write_which': self.write_which,
                  'image': im,
                  'out_prefix': self.output_prefix,
                  'write_mask': self.write_mask,
                  'graphics_fname': self.graphics_fname,
                  'params_fname': self.realignment_params,
                  'save_params': self.save_params}

        Processes = ASL_processes()
        exclude = Processes.Realigning_process(params)

        # remove realignment params if not desired
        if not params['save_params']:
            params_path = os.path.relpath(self.realignment_params, self.project.folder)
            self.project.session.remove_document(COLLECTION_INITIAL, params_path)
            self.project.session.remove_document(COLLECTION_CURRENT, params_path)

        # add new tag to all the scans if the tag does not exist
        ASL_path = os.path.relpath(self.ASL_realigned, self.project.folder)
        if 'Excluded' not in self.project.session.get_fields_names(COLLECTION_INITIAL):
            new_tag_name = 'Excluded'
            new_default_value = '[]'
            tag_type = 'list_int'
            new_tag_description = 'nothing'
            new_tag_unit = None
            self.add_new_tag(new_tag_name, new_default_value, tag_type,
                             new_tag_description, new_tag_unit)

        exclude = [3,6,47,2]
        # change value if there are excluded frames
        if exclude != []:
            document = ASL_path
            field = 'Excluded'
            new_value = exclude
            self.project.session.set_value(COLLECTION_INITIAL, document, field, new_value, flush=True)
            self.project.session.set_value(COLLECTION_CURRENT, document, field, new_value, flush=True)

class ASL_Outlier_Detection(Process_Mia, ASL_processes):
    """
    Detects outliers.
    Methods:
        - list_outputs: generates the outputs of the process
        - run_process_mia: runs Outlier_detection_process
    """

    def __init__(self):
        super(ASL_Outlier_Detection, self).__init__()

        # inputs
        self.add_trait("im_file", OutputMultiPath(traits.Either(
            traits.List(traits.File(output = False, optional = False, exists = True)),
            traits.File(output = False, optional = False, exists = True))))
        self.add_trait("exclude_in", traits.File(output = False, optional = True))
        self.add_trait("exclude_out", traits.File(output = False, optional = True))
        self.add_trait("mask_file", traits.File('', output = False, optional = True))
        self.add_trait("mask_threshold", traits.Float(0.5, output = False, optional = True))
        self.add_trait("plot_diagnostics", traits.Bool(True, output = False, optional = True))
        self.add_trait("outlier_threshold", traits.Float(0.0005, output = False, optional = True))

        # outputs
        self.add_trait("output", traits.File(output = True))

    def list_outputs(self):
        super(ASL_Outlier_Detection, self).list_outputs()
        if not self.im_file or self.im_file in ["<undefined>", traits.Undefined]:
            return {}, {}
        path, filename = os.path.split(self.im_file)
        filename = "outlier_" + filename
        self.output = os.path.join(path, filename)
        return {"files_excl": self.im_file}, {self.output: self.im_file}

    def run_process_mia(self):
        super(ASL_Outlier_Detection, self).run_process_mia()

        if isinstance(self.im_file, list):
            print('Reducing list of files to single file.')
            self.im_file = self.im_file[0]

        params = {'im_file': self.im_file,
                 'exclude_in': self.exclude_in, #ASK JAN: does this to params make sense in this context? Because we are adding TAGS
                 'exclude_out': self.exclude_out,
                 'mask_file': self.mask_file,
                 'mask_threshold': self.mask_threshold,
                 'plot_diagnostics': self.plot_diagnostics,
                 'outlier_threshold': self.outlier_threshold}

        Processes = ASL_processes()
        exclude = Processes.Outlier_detection_process(params)

        # check if Exclude tag exists, if not, create it
        ASL_path = os.path.relpath(self.im_file, self.project.folder)
        if 'Excluded' not in self.project.session.get_fields_names(COLLECTION_INITIAL):
            new_tag_name = 'Excluded'
            new_default_value = '[]'
            tag_type = 'list_int'
            new_tag_description = 'nothing'
            new_tag_unit = None
            self.add_new_tag(new_tag_name, new_default_value, tag_type,
                             new_tag_description, new_tag_unit)
        exclude = [3,45]
        # change value if there are excluded frames
        if exclude != []:
            excl_tag = self.project.session.get_value(COLLECTION_INITIAL, ASL_path, "Excluded")
            document =  ASL_path
            field = 'Excluded'
            new_value = list(set(exclude + excl_tag))
            self.project.session.set_value(COLLECTION_INITIAL, document, field, new_value, flush=True)
            self.project.session.set_value(COLLECTION_CURRENT, document, field, new_value, flush=True)

class ASL_Design_Matrix_Generation(Process_Mia, ASL_processes):
    """
    Methods:
        - list_outputs: generates an object with the design matrix
        - run_process_mia: runs Generate_design_matrix_process
    """

    def __init__(self):
        super(ASL_Design_Matrix_Generation, self).__init__()

        # Inputs
        self.add_trait("ASL_M0_files", traits.File(output = False, optional = False))
        self.add_trait("directory", traits.Str('/home/nietob/Data/mia/projects/test2/data/derived_data', output = False, optional = True))
        self.add_trait("timing_units", traits.Str('scans', output = False, optional = True))
        self.add_trait("exclude_fname", traits.File(output = False, optional = True)) #optional??
        self.add_trait("phys_reg_file", traits.File(output = False, optional = True)) #optional??
        self.add_trait("build_BOLD_phys", traits.Bool(False, output = False, optional = True)) #optional??
        self.add_trait("HPF_cutoff", traits.Float(64.0, output = False, optional = True))
        self.add_trait("CBF_regressor_name", traits.Str('baseline_CBF', output = False, optional = True))
        self.add_trait("explicit_mask", traits.File('', output = False, optional = True))

        # Outputs
        self.add_trait("matrix_object", traits.File(output=True))
        self.add_trait("ASL_M0", traits.File(output=True)) #still return ASL and M0 images?? or take them later from AL_out_det?

    def list_outputs(self): #TODO change this in order to return an object
        super(ASL_Design_Matrix_Generation, self).list_outputs()
        if not self.ASL_M0_files or self.ASL_M0_files in ["<undefined>", traits.Undefined]:
            return {}, {}
        path, filename = os.path.split(self.ASL_M0_files)
        filename = "rotated_" + filename
        self.matrix_object = os.path.join(path, filename)
        # return {"matrix_object": self.matrix_object}, {}
        return {}, {}

    def run_process_mia(self):
        super(ASL_Design_Matrix_Generation, self).run_process_mia()

        # take the label values
        PathName = '/'.join(self.ASL_M0_files.split('/')[-3:])
        self.TR = self.project.session.get_value(COLLECTION_INITIAL,PathName, "RepetitionTime")

        try:
            self.excluded = self.project.session.get_value(COLLECTION_INITIAL,PathName, "Excluded")
        except TypeError:
            self.excluded = []

        im = nib.load(self.ASL_M0_files)
        params = {"ASL_M0_files": self.ASL_M0_files,
                  "directory": self.directory,
                  "timing_units": self.timing_units,
                  "exclude_fname": self.exclude_fname,
                  "phys_reg_file": self.phys_reg_file,
                  "build_BOLD_phys": self.build_BOLD_phys,
                  "HPF_cutoff": self.HPF_cutoff,
                  "CBF_regressor_name": self.CBF_regressor_name,
                  "explicit_mask": self.explicit_mask,
                  "TR": self.TR,
                  "excluded": self.excluded}

        Processes = ASL_processes()
        P = Processes.Generate_design_matrix_process(params)

class ASL_TMap_Generation(Process_Mia, ASL_processes):
    """
    Methods:
        - list_outputs: generates the outputs of the process
        - run_process_mia: runs Generate_TMap_process
    """
    def __init__(self):
        super(ASL_TMap_Generation, self).__init__()
        # inputs
        self.add_trait("SPM_file", traits.File(output=False, optional=False))
        self.add_trait("regressor_filename", traits.File(output=False, optional=True))
        self.add_trait("regressor_to_contrast", traits.Str(output=False, optional=False))
        self.add_trait("delete_existing_contrasts", traits.Bool(False, output=False, optional=True))

        # outputs
        self.add_trait("output", traits.File(output=True))

    def list_outputs(self):
        super(ASL_TMap_Generation, self).list_outputs()
        if not self.SPM_file or self.SPM_file in ["<undefined>", traits.Undefined]:
            return {}, {}
        path, filename = os.path.split(self.SPM_file)
        filename = "Tmap_" + filename
        self.output = os.path.join(path, filename)
        return {"output": self.output}, {}

    def run_process_mia(self):
        super(ASL_TMap_Generation, self).run_process_mia()

        params = {"SPM_file": self.SPM_file,
                  "regressor_filename": self.regressor_filename,
                  "regressor_to_contrast": self.regressor_to_contrast,
                  "delete_existing_contrasts": self.delete_existing_contrasts}

        Processes = ASL_processes()
        P = Processes.ASL_TMap_Generation(params)

class ASL_CBFmap_Generation(Process_Mia, ASL_processes):
    """
    Methods:
        - list_outputs: generates the outputs of the process
        - run_process_mia: runs Calculate_CBFmap_process
    """
    def __init__(self):
        return None
    def list_outputs(self):
        return None
    def run_process_mia(self):
        return None

class ASL_Mask_Generation(Process_Mia, Anat):
    """
    Calls Generate_mask_process from Anat class to create masks.
    Methods:
        - list_outputs: generates the outputs of the process
        - run_process_mia: runs the functions needed, Resample_image, Get_Def, Apply_Def
    """

    def __init__(self):
        super(ASL_Mask_Generation, self).__init__()

        # Inputs
        # self.add_trait("masks", traits.List(['FRON','PAR','TEMP','OCC','CING','INSULA'], output=False, optional = True))
        self.add_trait("masks", traits.File('/home/nietob/Data/mia/projects/anat/data/raw_data/brainmask.nii', output=False, optional = True)) #TEMPORAL
        self.add_trait("mask_directory", traits.Str('/home/nietob/Data/template', output = False, optional = True))  # this should be not optional (I guess) but its termporal
        self.add_trait("normalised_data", traits.Bool(True, output = False, optional = True))
        self.add_trait("normalisation_parameters", traits.File(output = False, optional = True))
        self.add_trait("reference_volume", traits.File('', output = False, optional = True)) #filename or file?? #TODO ask jibril if i can put a defaul value of a file as ''
        self.add_trait("tissue_map", traits.File(output = False, optional = False))
        self.add_trait("tissue_map_threshold", traits.Float(0.9, output = False, optional = True))

        # Outputs
        self.add_trait("output_files", traits.File(output = True))

    def list_outputs(self):
        super(ASL_Mask_Generation, self).list_outputs()
        return {}, {}

    def run_process_mia(self):
        super(ASL_Mask_Generation, self).run_process_mia()

        # parameters handling
        if self.normalised_data == False and not self.normalisation_parameters:
            print('When the data is not normalised, the normalisation parameters are mandatory.')
            return

        else:
            tissue_m = nib.load(self.tissue_map)
            if self.reference_volume != '':
                ref_vol = nib.load(self.reference_volume)
                out_dir = os.path.dirname(self.reference_volume)
            else:
                out_dir = os.path.dirname(self.tissue_map)

        base = os.path.basename(self.tissue_map)
        tiss_name = os.path.splitext(base)[0]

        if 'c1' in tiss_name:
            tiss_name = '_GM'
        elif 'c2' in tiss_name:
            tiss_name = '_WM'
        elif 'c3' in tiss_name:
            tiss_name = '_CSF'
        else:
            tiss_name = ''

        out_name = tiss_name + '.nii'
        tissue_map_data = tissue_m.get_fdata()
        print(self.masks)
        # create params
        params = {"masks": self.masks,
                  "mask_directory": self.mask_directory,
                  "normalised_data": self.normalised_data,
                  "normalisation_parameters": self.normalisation_parameters,
                  "reference_volume": self.reference_volume,
                  "tissue_map_threshold": self.tissue_map_threshold,
                  "tissue_map_data": tissue_map_data}

        # normalise and create mask
        Processes = Anat()
        P = Processes.Generate_mask_process(params)


        return None





#IN DATA_BROWSER LINE 848
# self.project.session.get_fields_names(COLLECTION_CURRENT)
# tags.remove(TAG_CHECKSUM)
# tags.remove(TAG_FILENAME)
# tags = sorted(tags)

# TO SHOW ALL THE DATA IN ANN ARRAY
#         import sys
#         import numpy
#         numpy.set_printoptions(threshold=sys.maxsize)
#
# import scipy.io as spio
# mat = spio.loadmat('data.mat', squeeze_me=True)

# CONVERT .IMG TO .NII

# import nibabel as nb
# fname = '<file name>.img'
# img = nb.load(fname)
# nb.save(img, fname.replace('.img', '.nii'))



# TO CREATE TAGS
# self.project.session.add_field(COLLECTION_CURRENT, new_tag_name,
#                                tag_type, new_tag_description, True,
#                                TAG_ORIGIN_USER, new_tag_unit,
#                                new_default_value)
# self.project.session.add_field(COLLECTION_INITIAL, new_tag_name,
#                                tag_type, new_tag_description, True,
#                                TAG_ORIGIN_USER, new_tag_unit,
#                                new_default_value)