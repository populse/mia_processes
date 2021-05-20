"""The other preprocess library of the mia_processe package.

The purpose of this module is to provide bricks generally necessary for the 
pre-processing steps, which are not found in nipype.

:Contains:
    :Class:
        - Resample
        - Threshold
    :Function:
        - threshold

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
import nibabel.processing as nibp

# nipype import
from nipype.interfaces.base import (OutputMultiPath, InputMultiPath, File,
                                    traits, TraitListObject, Undefined)
from nipype.interfaces.spm.base import ImageFileSPM

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# soma-base imports
from soma.qt_gui.qt_backend.Qt import QMessageBox

# Other import
import os


class Resample(ProcessMIA):
    """
    *Resamples an image to the resolution of a reference image*

    Please, see the complete documention for the `Resample brick in the populse.mia_processes web site
    https://populse.github.io/mia_processes/documentation/bricks/preprocess/other/Resample.html

   """
    
    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Resample, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        reference_image_desc = ('A 3D or 4D image used as reference to '
                                'resample the files_to_resample images (a '
                                'pathlike object or string representing a file '
                                'with extension in [.img, .nii, .hdr]).')
        
        files_to_resample_desc = ('The 3D images that will be resampled (a '
                                  'list of pathlike object or string '
                                  'representing a file or a list of items '
                                  'which are a pathlike object or string '
                                  'representing a file, valid extensions: '
                                  '[.img, .nii, .hdr]).')
        suffix_to_delete_desc = ('The suffix to delete from the '
                                 'files_to_resample, when creating the '
                                 'out_files (a string).')
        suffix_desc = 'The suffix for the out_files image (a string).'
        prefix_desc = 'The prefix for the  out_files image (a string).'
        interp_desc = ('The order of the spline interpolation (an integer '
                       'between 0 and 5; trilinear == 3).')
        
        # Outputs description
        out_files_desc = ('The resulting image after resampling (a pathlike '
                         'object or string representing a file, or a list of '
                         'pathlike objects or strings representing a file).')
        
        # Inputs traits

        self.add_trait("reference_image",
                       ImageFileSPM(output=False,
                                    desc=reference_image_desc))

        self.add_trait("files_to_resample",
                       InputMultiPath(ImageFileSPM(),
                                      output=False,
                                      desc=files_to_resample_desc))
        
        self.add_trait("suffix_to_delete",
                       traits.String("_002",
                                     output=False,
                                     optional=True,
                                     desc=suffix_to_delete_desc))

        self.add_trait("suffix",
                       traits.String("_003",
                                     output=False,
                                     optional=True,
                                     desc=suffix_desc))
        
        self.add_trait("prefix",
                       traits.String(" ",
                                     output=False,
                                     optional=True,
                                     desc=prefix_desc))
        
        self.add_trait("interp",
                       traits.Range(value=3,
                                    low=0,
                                    high=5,
                                    output=False,
                                    optional=True,
                                    desc=interp_desc))

        # Outputs traits
        self.add_trait("out_files",
                       OutputMultiPath(File(),
                                       output=True,
                                       desc=out_files_desc))

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
        super(Resample, self).list_outputs()

        # Outputs definition
        if ((self.reference_image) and
                (self.files_to_resample) and
                (self.interp)):
            files = []
            files_name = self.files_to_resample

            try:
                refName = nib.load(self.reference_image)
                
            except (nib.filebasedimages.ImageFileError,
                    FileNotFoundError, TypeError) as e:
                print("\nError with reference_image, during initialisation: ",
                      e)
                refName = None
            
            if ((not self.suffix_to_delete) or
                    (self.suffix_to_delete.isspace()) or
                    (self.suffix_to_delete in [Undefined, "<undefined>"])):
                self.suffix_to_delete = " "

            if ((not self.suffix) or
                    (self.suffix.isspace()) or
                    (self.suffix in [Undefined, "<undefined>"])):
                self.suffix = " "

            if ((not self.prefix) or
                    (self.prefix.isspace()) or
                    (self.prefix in [Undefined, "<undefined>"])):
                self.prefix = " "
           
            for file_name in files_name:

                try:
                    fileName = nib.load(file_name)

                except (nib.filebasedimages.ImageFileError,
                        FileNotFoundError, TypeError) as e:
                    print("\nError with files_to_resample, during "
                           "initialisation: ", e)
                    fileName = None

                if (refName) and (fileName):

                    if ((len(fileName.shape) != 3) or
                                               (3 > len(refName.shape) > 4)):
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Warning)
                        msg.setWindowTitle("mia_processes - "
                                           "Resample brick Warning!")
                        msg.setText("Currently, the Resample brick only allows "
                                    "to resample 3D images from 3D or 4D "
                                    "reference images.\n\n However the '{0}' "
                                    "reference image is a {1}D and the '{2}' "
                                    "image is a {3}D ...\n\nPlease, modify "
                                    "your input and initialise again this "
                                    "brick.".format(self.reference_image,
                                                    len(refName.shape),
                                                    file_name,
                                                    len(fileName.shape)))
                        msg.setStandardButtons(QMessageBox.Close)
                        msg.buttonClicked.connect(msg.close)
                        msg.exec()
                        print("\nResample brick: Initialisation failed ... "
                              "Please check the input parameters!")
                        files_name = None
                        break

                else:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Warning)
                    msg.setWindowTitle("mia_processes - "
                                       "Resample brick Warning!")
                    msg.setText("files_to_resample '{0}' or/and "
                                "reference_image '{1}' is (are) empty ... Do "
                                "you want to continue? \n\n"
                                "- To correct the input parameters "
                                "click 'Abort'.\n\n"
                                "- If the Resample brick is located in a "
                                "pipeline, it may be usual that this(these) "
                                "parameter(s) is(are) empty at the "
                                "initialisation time. In this case, if the "
                                "files_to_resample and the reference_image "
                                "parameters correspond to a 3D image and a 3D "
                                "or a 4D, respectively, click 'Yes', "
                                "otherwise click 'Abort' and change the input "
                                "parameters.".format(file_name,
                                                     self.reference_image))
                    msg.setStandardButtons(QMessageBox.Yes | QMessageBox.Abort)
                    retval = msg.exec_()
                    
                    if retval == QMessageBox.Yes:
                        print("\nResample brick Warning!: Empty input file(s). "
                              "No check of the dimensions of the input images "
                              "is done ... (files_to_resample must be a 3D and "
                              "reference_image must be a 3D or a 4D) ...")

                    else:
                         print("\nResample brick: Initialisation failed ... "
                               "Please check the input parameters!")
                         files_name = None
                         break
                
                path, file_name = os.path.split(file_name)
                file_name_no_ext, file_extension = os.path.splitext(file_name)

                if ((self.suffix_to_delete != " ") and
                        (file_name_no_ext[-len(self.suffix_to_delete):] ==
                                                        self.suffix_to_delete)):
                    file_name_no_ext = (
                                 file_name_no_ext[:-len(self.suffix_to_delete)])
                    
                files.append(os.path.join(self.output_directory,
                                          (self.prefix.strip() +
                                                 file_name_no_ext +
                                                 self.suffix.strip() +
                                                 file_extension)))
                
            # May overwrite an input image with an output image
            if (self.suffix == " ") and (self.prefix == " ") and (files_name):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle("mia_processes - "
                                   "Resample brick Warning!")
                msg.setText("Suffix and prefix input parameters are not "
                            "defined or consist only of one or more white "
                            "spaces.\n\nThe files_to_resample input parameter "
                            "(if suffix_to_delete input parameter is not set) "
                            "or other files (if suffix_to_delete input "
                            "parameter is set) could be overwritten ... this "
                            "concerns the following images:\n{}\n\nDo you "
                            "agree to use these input "
                            "parameters?".format(set(files)))
                msg.setStandardButtons(QMessageBox.Yes | QMessageBox.Abort)
                retval = msg.exec_()

                if retval == QMessageBox.Yes:
                    print('\nResample brick warning: suffix and prefix input '
                          'parameters are not defined, the following files '
                          'could be overwrite!:\n{0} ...\n'.format(set(files)))

                else:
                    files_name = None

            # Two (at least) input images give the same output image
            if (len(files) != len(set(files))) and (files_name):
                dupes = set()
                seen = set()

                for x in files:

                    if x in seen:
                        dupes.add(x)

                    seen.add(x)

                msg = QMessageBox()
                msg.setIcon(QMessageBox.Warning)
                msg.setWindowTitle("mia_processes - "
                                   "Resample brick Warning!")
                msg.setText("The suffix, prefix and suffix_to_delete input "
                            "parameters combination seems to create for at "
                            "least two different input images the same output "
                            "image:\n{}\n\nPlease change the suffix, prefix "
                            "and suffix_to_delete input parameters combination "
                            "and initialise again this brick.".format(dupes))
                msg.setStandardButtons(QMessageBox.Close)
                msg.buttonClicked.connect(msg.close)
                msg.exec()
                files_name = None

            if files_name:
                self.outputs['out_files'] = files

        #Tags inheritance (optional)
        if self.outputs:

            for key, val in self.outputs.items():

                if key == "out_files":

                    for in_val, out_val in zip(files_name, val):
                        _, fileOval = os.path.split(out_val)
                        fileOval_no_ext, _ = os.path.splitext(fileOval)

                        if ((self.suffix != " ") and
                                (fileOval_no_ext[-len(self.suffix):] ==
                                                                  self.suffix)):
                            fileOval_no_ext = fileOval_no_ext[:-len(self.
                                                                        suffix)]

                        if ((self.prefix != " ") and
                                (fileOval_no_ext[0:len(self.prefix)] ==
                                                                  self.prefix)):
                            fileOval_no_ext = (fileOval_no_ext[len(self.
                                                                      prefix):])

                        _, fileIval = os.path.split(in_val)
                        fileIval_no_ext, _ = os.path.splitext(fileIval)

                        if ((self.suffix_to_delete != " ") and
                            (fileIval_no_ext[-len(self.suffix_to_delete):] ==
                                                        self.suffix_to_delete)):
                            fileOval_no_ext = (fileOval_no_ext +
                                               self.suffix_to_delete)

                        if fileOval_no_ext == fileIval_no_ext:
                            self.inheritance_dict[out_val] = in_val

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def _check_interp(self):
        """Checks the order of the splines interpolation.

        :returns: a message according to the selected order
        """

        # skimage.transform.resize (does not seem to be the best solution for
        # resampling):
        # - B-splines of the order 0 and 1 correspond to nearest neighbour and
        #   linear interpolation, respectively.
        # - B-splines of a higher order can be defined by a repetitive
        #   convolution of the zeroth-order spline (the box function) with
        #   itself.
        #
        # nibabel.processing.resample_from_to (seems to be a good solution, but
        # can give non-zero intensity for the background -> needs thresolding ?)
        #
        # nilearn.image.resample_to_img (seems to be a good solution, no
        # non-zero intensity for the background observed)
        
        if self.interp == 0:
            # skimage and nibabel: Nearest neighbour
            return "spline interpolation of order 0 (Nearest-neighbour)"

        elif self.interp == 1:
            # skimage: It seems that with dimensions=3 and B-spline order=1,
            #          the interpolation would be equivalent to "trilinear"
            #return ("Bi-linear (with dimensions=3, Bi-linear is equivalent "
            #        "to 'trilinear')")

            # nibabel: spline interpolation of order 1
            return "spline interpolation of order 1"

        elif self.interp == 2:
            # skimage:
            # return "Bi-quadratic"

            # nibabel: spline interpolation of order 2
            return "spline interpolation of order 2"
        
        elif self.interp == 3:
            # skimage
            # return "Bi-cubic"

            # nibabel: spline interpolation of order 3: trilinear
            return "spline interpolation of order 3 (trilinear)"
        
        elif self.interp == 4:
            # skimage
            # return "Bi-quartic"

            # nibabel: spline interpolation of order 4
            return "spline interpolation of order 4"
         
        elif self.interp == 5:
            # skimage
            # return "Bi-quintic"

            # nibabel: spline interpolation of order 5
            return "spline interpolation of order 5"

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Resample, self).run_process_mia()
        files_name = self.files_to_resample
        refName = nib.load(self.reference_image)
        print('\nResample process from mia_processes: \n',
              self._check_interp(), ' \n')

        for file_name in files_name:
            fileName = nib.load(file_name)
            #mask_data = mask.get_fdata()  # skimage
              
            # Currently, resampling is only done for 3D images using 3D
            # or 4D images for reference. It would be interesting to widen the
            # possibilities to more cases (at least a 4D using a 3D, and why
            # not, to larger dimensions?). Currently we use
            # nibabel.processing.resample_from_to(), (see
            # https://mail.python.org/pipermail/neuroimaging/2019-January/001902.html).
            # It seems this method produces a little noise (can give non-zero
            # intensity for the background). In addition to the resolution
            # settings it seems the size of the reference_image are also
            # transferred to the resampled image but not the orientations and
            # positions (output space comes from the affine of
            # files_to_resample).

            if len(fileName.shape) == len(refName.shape) == 3:
                # nibabel:
                fileFinal = nibp.resample_from_to(fileName,
                                                  refName,
                                                  order=self.interp)

            if len(fileName.shape) == 3 and len(refName.shape) == 4:
                # skimage
                #ref_size = ref_data.shape[:3]
                #resized_mask_data = resize(mask_data,
                #                           ref_size,
                #                           order=self.interp,
                #                           mode='reflect')
                # TODO: Taking info of mask's or ref's header?
                #mask_final = nib.Nifti1Image(resized_mask_data,
                #                             ref.affine,
                #                             ref.header)

                # nibabel:
                fileFinal = nibp.resample_from_to(fileName,
                                                  refName.slicer[:,:,:,0],
                                                  order=self.interp)

            # Image save
            path, file_name = os.path.split(file_name)
            file_name_no_ext, file_extension = os.path.splitext(file_name)

            if ((self.suffix_to_delete != " ") and
                        (file_name_no_ext[-len(self.suffix_to_delete):] ==
                                                        self.suffix_to_delete)):
                file_name_no_ext = (
                                 file_name_no_ext[:-len(self.suffix_to_delete)])
                    
            file_out= (os.path.join(self.output_directory,
                                    (self.prefix.strip() +
                                     file_name_no_ext +
                                     self.suffix.strip() +
                                     file_extension)))
            nib.save(fileFinal, file_out)


class Threshold(ProcessMIA):
    """
    *Image thresholding*

    Please, see the complete documention for the `Threshold brick in the populse.mia_processes web site
    https://populse.github.io/mia_processes/documentation/bricks/preprocess/other/Threshold.html

    """
    def __init__(self):
        """Dedicated to the attributes initialisation / instanciation.
        
        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(Threshold, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = []

        # Inputs description
        in_files_desc = ('A list of items with string elements corresponding '
                         'to existing path files.')
        threshold_desc = ('Value for the applied threshold (a float between 0 '
                          'and 1).')
        GM_filtering_desc = ('Filtering (based on the keyword c1) to keep only '
                             'grey matter images (a boolean)')
        suffix_desc = 'Suffix of the output image (a string).'
        prefix_desc = 'Prefix of the output image (a string).'
        
        # Outputs description
        out_files_desc = ('Path of the scan after application of the threshold '
                         '(a pathlike object or string representing a file, or '
                         'a list of pathlike objects or strings representing a '
                         'file).')
        
        # Inputs traits
        self.add_trait("in_files",
                       InputMultiPath(traits.Either(ImageFileSPM(),
                                                   traits.List(ImageFileSPM())),
                                      output=False,
                                      desc=in_files_desc))
        
        self.add_trait("threshold",
                       traits.Range(value=0.3,
                                    low=0.0,
                                    high=1.0,
                                    output=False,
                                    optional=True,
                                    desc=threshold_desc))

        self.add_trait("GM_filtering",
                       traits.Bool(default_value=True,
                                   output=False,
                                   optional=True,
                                   desc=GM_filtering_desc))
        
        self.add_trait("suffix",
                       traits.String("_002",
                                     output=False,
                                     optional=True,
                                     desc=suffix_desc))
        
        self.add_trait("prefix",
                       traits.String("",
                                     output=False,
                                     optional=True,
                                     desc=prefix_desc))

        # Outputs traits
        self.add_trait("out_files",
                       OutputMultiPath(File(),
                                       output=True,
                                       desc=out_files_desc))

        self.init_default_traits()

    def _namesFilter(self):
        """Filtering of in_files input parameter in order to keep only GM.

        :returns: GM images (containing "c1" in the first 5 characters)
        """
        files = []

        for file_name in self.in_files:

            if isinstance(file_name, (list, TraitListObject)):
                file_name = file_name[0]

            # Take the first 5 characters in the case if the GM was processed
            # ex. normalisation = wc1, normalisation then smooth = swc1 ...
            # This is not the cleaner way ... in case of a file name
            # (PatientName) starting with c1 !
            if 'c1' in os.path.basename(os.path.normpath(file_name))[:5]:
                files.append(file_name)

        return files

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
        super(Threshold, self).list_outputs()

        # Outputs definition
        if ((self.in_files) and (self.threshold)):

            if self.GM_filtering is True:
                files_name = self._namesFilter()

            else:
                files_name = self.in_files

            if ((not self.suffix) or (self.suffix.isspace()) or
                                     self.suffix in [Undefined, "<undefined>"]):
                self.suffix = " "

            if ((not self.prefix) or (self.prefix.isspace()) or
                                   (self.prefix in [Undefined, "<undefined>"])):
                self.prefix = " "

            files = []
            flag = True # If False, suf/pref check will not be performed later
            retval = None

            for file_name1 in files_name:
                path, file_name = os.path.split(file_name1)

                if (self.suffix == " " and
                        self.prefix == " " and
                            path == self.output_directory and flag is True ):
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Warning)
                    msg.setWindowTitle("mia_processes - "
                                       "Threshold brick Warning!")
                    msg.setText("Suffix and prefix input parameters are not "
                                "defined or consist only of one or more white "
                                "spaces.\nThe {0} input parameter will be "
                                "overwritten ...\n Yes or "
                                "Abort?".format(file_name1))
                    msg.setStandardButtons(QMessageBox.Yes|QMessageBox.YesToAll|
                                           QMessageBox.Abort)
                    retval = msg.exec_()

                    if retval !=  QMessageBox.Abort:
                        (file_name_no_ext,
                         file_extension) = os.path.splitext(file_name)
                        files.append(os.path.join(self.output_directory,
                                                  (self.prefix.strip() +
                                                   file_name_no_ext +
                                                   self.suffix.strip() +
                                                   file_extension)))
                        print('\nThreshold brick warning: the out_files output '
                              'parameter is the same as the in_files input '
                              'parameter (suffix and prefix are not defined):'
                              '\n{0} will be overwrited ...'.format(file_name1))

                        if retval == QMessageBox.YesToAll:
                            flag = False
                            print('\nYesToAll selected: end of overwrite '
                                  'checks on input images ...')
                    else:
                        files_name = []
                        print('\nAborted. Please check your input parameters ...')
                        break

                else:
                    (file_name_no_ext,
                         file_extension) = os.path.splitext(file_name)
                    files.append(os.path.join(self.output_directory,
                                              (self.prefix.strip() +
                                               file_name_no_ext +
                                               self.suffix.strip() +
                                               file_extension)))

            if files_name:
                self.outputs['out_files'] = files

            else:
                print('- There was no output file deducted during '
                      'initialisation. Please check the input parameters...!')

        # tags inheritance (optional)
        if self.outputs:

            for key, val in self.outputs.items():

                if (key == "out_files"):

                    for in_val, out_val in zip(files_name, val):
                        _, fileOval = os.path.split(out_val)
                        fileOval_no_ext, _ = os.path.splitext(fileOval)
                        _, fileIval = os.path.split(in_val)
                        fileIval_no_ext, _ = os.path.splitext(fileIval)

                        if ((self.prefix) and
                                     (fileOval_no_ext.startswith(self.prefix))):
                            fileOval_no_ext = fileOval_no_ext[len(self.prefix):]

                        if ((self.suffix) and
                                       (fileOval_no_ext.endswith(self.suffix))):
                            fileOval_no_ext = fileOval_no_ext[:-len(
                                                                   self.suffix)]
                        
                        if fileOval_no_ext == fileIval_no_ext:
                            self.inheritance_dict[out_val] = in_val

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Threshold, self).run_process_mia()
        
        if self.GM_filtering is True:
            files_name = self._namesFilter()

        else:
            files_name = self.in_files

        for file_name in files_name:
            # Image processing
            img_final = threshold(file_name, self.threshold)

            # Image save
            _, file_name = os.path.split(file_name)
            file_name_no_ext, file_extension = os.path.splitext(file_name)
            out_file = os.path.join(self.output_directory,
                                    (self.prefix.strip() +
                                     file_name_no_ext +
                                     self.suffix.strip() +
                                     file_extension))
            nib.save(img_final, out_file)

def threshold(file_name, thresh):
    """
    Basic method for image thresholding

    :param file_name: Image to be thresholded
    :param thresh: Threshold value (a float between 0 and 1)
    :returns: Image after thresholding
    """
    img = nib.load(file_name)
    img_data = img.get_fdata()
    _max =  img_data.max()
    img_thresh = (img_data > (_max*thresh)).astype(float)
    img_final = nib.Nifti1Image(img_thresh, img.affine, img.header)
    return img_final
