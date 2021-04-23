"""The other preprocess library of the mia_processe package.

The purpose of this module is to provide bricks generally necessary for the 
pre-processing steps, which are not found in nipype.

:Contains:
    :Class:
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
                    msg.buttonClicked.connect(msg.close)
                    retval = msg.exec()

                    if retval == QMessageBox.YesToAll: flag = False

                if retval !=  QMessageBox.Abort:
                    (file_name_no_ext,
                     file_extension) = os.path.splitext(file_name)
                    files.append(os.path.join(self.output_directory,
                                              (self.prefix.strip() +
                                               file_name_no_ext +
                                               self.suffix.strip() +
                                               file_extension)))

                    if (path == self.output_directory and
                                   self.suffix == " " and self.prefix == " "):
                        print('\nThreshold brick warning: The out_files output '
                              'parameter is the same as the in_files input '
                              'parameter (suffix and prefix are not defined):\n'
                              '{0} will be overwrited ...'.format(file_name1))

                else:
                    files_name = []
                    break

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
