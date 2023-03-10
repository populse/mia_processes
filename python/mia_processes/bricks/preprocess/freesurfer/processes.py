"""The freesurfer preprocess library of the mia_processes package.

The purpose of this module is to customise the main freesurfer preprocessing bricks
provided by nipype and to correct some things that do not work directly in
populse_mia.

:Contains:
    :Class:
        - Binarize
        - SynthStrip

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# populse_mia import
from populse_mia.user_interface.pipeline_manager.process_mia import ProcessMIA

# Other import
import os
from traits.api import Bool, Either, Enum, File, Float, Int, List, Undefined, String 
import capsul
from capsul.in_context import freesurfer

EXT = {'NIFTI_GZ': 'nii.gz',
       'NIFTI': 'nii',
       'MGZ': 'mgz'}


class Binarize(ProcessMIA):
    """
    * Use FreeSurfer mri_binarize to binarize a volume
    (or volume-encoded surface file).
    Can also be used to merge with other binarizations.
    Binarization can be done based on threshold or on matched values. *

    Please, see the complete documentation for the `Binarize' brick in the populse.mia_processes website
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/freesurfer/Binarize.html>`

    """

    def __init__(self):
        '''Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        '''
        # Initialisation of the objects needed for the launch of the brick
        super(Binarize, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['freesurfer', 'nipype']

        # Inputs description
        in_file_desc = ('Input file (a pathlike object or string '
                        'representing a file).')
        min_desc = 'Minimum voxel threshold(float).'
        max_desc = 'Maximum voxel threshold(float).'
        rmin_desc = 'Compute min based on rmin*globalmean.'
        rmax_desc = 'Compute max based on rmax*globalmean.'
        match_desc = 'Match instead of threshold'
        wm_desc = 'set match vals to 2 and 41 (aseg for cerebral WM)'
        ventricles_desc = ('set match vals those for aseg '
                           'ventricles+choroid (not 4th)')
        wm_ven_csf_desc = ('WM and ventricular CSF,'
                           'including choroid (not 4th)')
        count_file_desc = ('save number of hits in ascii file'
                           '(hits, ntotvox, pct)')
        bin_val_desc = 'set vox outside range to val (default is 0)'
        bin_val_not_desc = 'set vox outside range to val (default is 0)'
        invert_desc = 'set binval=0, binvalnot=1'
        frame_no_desc = 'use 0-based frame of input (default is 0)'
        merge_file_desc = 'merge with mergevol'
        mask_file_desc = 'must be within mask'
        mask_thresh_desc = 'set thresh for mask'
        abs_desc = 'take abs of invol first (ie, make unsigned)'
        bin_col_num_desc = 'set binarized voxel value to its column number'
        zero_edges_desc = 'zero the edge voxels'
        zero_slice_edge_desc = 'zero the edge slice voxels'
        dilate_desc = 'niters: dilate binarization in 3D'
        erode_desc = ('nerode: erode binarization in 3D '
                      '(after any dilation)')
        erode2d_desc = ('nerode2d: erode binarization in 2D '
                        '(after any 3D erosion)')

        output_type_desc = ('Typecodes of the output image formats (one '
                            'of NIFTI, MGZ, NIFTI_GZ).')
        out_suffix_desc = 'Suffix of the output image (a string).'

        # Outputs description
        out_file_desc = ('The binanized file (a pathlike object or a '
                         'string representing a file).')

        # Inputs traits
        self.add_trait('in_file',
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait('min',
                       Either(Undefined,
                              Float(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=min_desc))

        self.add_trait('max',
                       Either(Undefined,
                              Float(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=max_desc))
        self.add_trait('rmin',
                       Either(Undefined,
                              Float(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=rmin_desc))

        self.add_trait('rmax',
                       Either(Undefined,
                              Float(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=rmax_desc))

        self.add_trait('match',
                       Either(Undefined,
                              List(Int()),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=match_desc))

        self.add_trait('wm',
                       Bool(default=False,
                            output=False,
                            optional=True,
                            desc=wm_desc))

        self.add_trait('ventricles',
                       Bool(default=False,
                            output=False,
                            optional=True,
                            desc=ventricles_desc))

        self.add_trait('wm_ven_csf',
                       Bool(default=False,
                            output=False,
                            optional=True,
                            desc=wm_ven_csf_desc))

        self.add_trait('count_file',
                       Either(Undefined,
                              File(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=count_file_desc))

        self.add_trait('output_type',
                       Enum('NIFTI',
                            'NIFTI_GZ',
                            'MGZ',
                            output=False,
                            optional=True,
                            desc=output_type_desc))

        self.add_trait('out_suffix',
                       String('_thresh',
                              output=False,
                              optional=True,
                              desc=out_suffix_desc))

        self.add_trait('bin_val',
                       Either(Undefined,
                              Int(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=bin_val_desc))

        self.add_trait('bin_val_not',
                       Either(Undefined,
                              Int(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=bin_val_not_desc))

        self.add_trait('invert',
                       Bool(default=False,
                            output=False,
                            optional=True,
                            desc=invert_desc))

        self.add_trait('frame_no',
                       Either(Undefined,
                              Int(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=frame_no_desc))

        self.add_trait('merge_file',
                       Either(Undefined,
                              File(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=merge_file_desc))

        self.add_trait('mask_file',
                       Either(Undefined,
                              File(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=mask_file_desc))

        self.add_trait('mask_thresh',
                       Either(Undefined,
                              Float(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=mask_thresh_desc))

        self.add_trait('abs',
                       Bool(default=False,
                            output=False,
                            optional=True,
                            desc=abs_desc))

        self.add_trait('bin_col_num',
                       Bool(default=False,
                            output=False,
                            optional=True,
                            desc=bin_col_num_desc))

        self.add_trait('zero_edges',
                       Bool(default=False,
                            output=False,
                            optional=True,
                            desc=zero_edges_desc))

        self.add_trait('zero_slice_edge',
                       Bool(default=False,
                            output=False,
                            optional=True,
                            desc=zero_slice_edge_desc))

        self.add_trait('dilate',
                       Either(Undefined,
                              Int(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=dilate_desc))

        self.add_trait('erode',
                       Either(Undefined,
                              Int(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=erode_desc))

        self.add_trait('erode2d',
                       Either(Undefined,
                              Int(),
                              default=Undefined,
                              output=False,
                              optional=True,
                              desc=erode2d_desc))
        # Outputs traits
        self.add_trait('out_file',
                       File(output=True,
                            desc=out_file_desc))

        self.init_default_traits()
        self.init_process('nipype.interfaces.freesurfer.Binarize')

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
        super(Binarize, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if (self.min != Undefined or
                self.max != Undefined) and self.match != Undefined:
            print('\nInitialisation failed. "match" parameter can not be used '
                  'with "min" and/or "max" parameters'
                  ' Please, define only "min"/"max" paremeters or "match"'
                  'parameters (set the other as Undefined) ...!')
            return
        if (self.min != Undefined or
                self.max != Undefined) and self.wm_ven_csf:
            print('\nInitialisation failed. "wm_ven_csf" parameter can not'
                  'be used  with "min" and/or "max" parameters'
                  ' Please, define only "min"/"max" paremeters or "wm_ven_csf"'
                  'parameters (set the other as Undefined) ...!')
            return
        if self.in_file:
            if self.output_directory:
                output_type = self.output_type

                try:
                    ifile = os.path.split(self.in_file)[-1]
                    file_name, in_ext = ifile.rsplit('.', 1)
                    if in_ext == 'gz':
                        (file_name_2, in_ext_2) = file_name.rsplit('.', 1)
                        if in_ext_2 == 'nii':
                            in_ext = 'nii.gz'
                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                self.outputs['out_file'] = os.path.join(
                    self.output_directory,
                    os.path.split(self.in_file)[1].replace(
                        '.' + in_ext,
                        self.out_suffix + '.' + EXT[output_type]))

                self.inheritance_dict[self.outputs[
                    'out_file']] = self.in_file
            else:
                print('No output_directory was found...!\n')
                return

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(Binarize, self).run_process_mia()

        # madatory inputs / outputs
        self.process.in_file = self.in_file
        self.process.binary_file = self.out_file
        self.process.out_type = EXT[self.output_type]

        # optionnal inputs
        self.process.min = self.min
        self.process.max = self.max
        self.process.rmin = self.rmin
        self.process.rmax = self.rmax
        self.process.match = self.match
        self.process.count_file = self.count_file
        self.process.bin_val = self.bin_val
        self.process.bin_val_not = self.bin_val_not
        self.process.frame_no = self.frame_no
        self.process.merge_file = self.merge_file
        self.process.mask_file = self.mask_file
        self.process.mask_thresh = self.mask_thresh
        self.process.dilate = self.dilate
        self.process.erode = self.erode
        self.process.erode2d = self.erode2d
        self.process.wm = self.wm
        self.process.ventricles = self.ventricles
        if self.wm_ven_csf:
            # only add wm_ven_csf when True because Mutually
            # exclusive with min/max
            self.process.wm_ven_csf = self.wm_ven_csf
        self.process.invert = self.invert
        self.process.abs = self.abs
        self.process.bin_col_num = self.bin_col_num
        self.process.zero_edges = self.zero_edges
        self.process.zero_slice_edge = self.zero_slice_edge

        return self.process.run(configuration_dict={})


class SynthStrip(ProcessMIA):
    """
    * Skul stripping using SynthStrip *

    SynthStrip is a skull-stripping tool that extracts brain signal
    from a landscape of image types, ranging across imaging modality,
    contrast, resolution, and subject population. It leverages a deep
    learning strategy  that synthesizes arbitrary training images
    from segmentation maps to optimize a robust model agnostic
    to acquisition specifics.

    SynthStrip: Skull-Stripping for Any Brain Image
    Andrew Hoopes, Jocelyn S. Mora, Adrian V. Dalca, Bruce Fischl*,
    Malte Hoffmann* (*equal contribution)
    NeuroImage 260, 2022, 119474
    https://doi.org/10.1016/j.neuroimage.2022.119474

    Please, see the complete documention for the `Segment' brick in the populse.mia_processes web site
    <https://populse.github.io/mia_processes/documentation/bricks/preprocess/freesurfer/SynthStrip.html>`_

    """

    def __init__(self):
        """Dedicated to the attributes initialisation/instantiation.

        The input and output plugs are defined here. The special
        'self.requirement' attribute (optional) is used to define the
        third-party products necessary for the running of the brick.
        """
        # Initialisation of the objects needed for the launch of the brick
        super(SynthStrip, self).__init__()

        # Third party softwares required for the execution of the brick
        self.requirement = ['freesurfer', 'nipype']

        # Inputs description
        in_file_desc = 'Input image to be brain extracted'
        border_mm_desc = ('Mask border threshold in mm')
        output_type_desc = ('Typecodes of the output image formats (one '
                            'of NIFTI, MGZ, NIFTI_GZ).')
        no_csf_desc = 'Exclude CSF from brain border'

        # Outputs description
        out_file_desc = 'Brain-extracted path'
        out_mask_desc = 'Brainmask path'

        # Inputs traits
        self.add_trait('in_file',
                       File(output=False,
                            optional=False,
                            desc=in_file_desc))

        self.add_trait('border_mm',
                       Int(1,
                           output=False,
                           optional=True,
                           desc=border_mm_desc))

        self.add_trait('output_type',
                       Enum('NIFTI',
                            'NIFTI_GZ',
                            'MGZ',
                            output=False,
                            optional=True,
                            desc=output_type_desc))

        self.add_trait('no_csf',
                       Bool(False,
                            output=False,
                            optional=True,
                            desc=no_csf_desc))

        # Outputs traits
        self.add_trait('out_file',
                       File(output=True,
                            optional=True,
                            desc=out_file_desc))

        self.add_trait('out_mask',
                       File(output=True,
                            desc=out_mask_desc))

        self.init_default_traits()

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
        super(SynthStrip, self).list_outputs()

        # Outputs definition and tags inheritance (optional)
        if self.in_file:
            if self.output_directory:
                output_type = self.output_type

                try:
                    ifile = os.path.split(self.in_file)[-1]
                    file_name, in_ext = ifile.rsplit('.', 1)
                    if in_ext == 'gz':
                        (file_name_2, in_ext_2) = file_name.rsplit('.', 1)
                        if in_ext_2 == 'nii':
                            in_ext = 'nii.gz'
                except ValueError:
                    print('\nThe input image format is not recognized ...!')
                    return

                self.outputs['out_file'] = os.path.join(
                    self.output_directory,
                    os.path.split(self.in_file)[1].replace(
                        '.' + in_ext,
                        '_desc-brain.' + EXT[output_type]))

                self.outputs['out_mask'] = os.path.join(
                    self.output_directory,
                    os.path.split(self.in_file)[1].replace(
                        '.' + in_ext,
                        '_desc-brain_mask.' + EXT[output_type]))

            else:
                print('No output_directory was found...!\n')
                return

        if self.outputs:
            self.inheritance_dict[self.outputs[
                'out_file']] = self.in_file
            self.inheritance_dict[self.outputs[
                'out_mask']] = self.in_file

        # Return the requirement, outputs and inheritance_dict
        return self.make_initResult()

    def run_process_mia(self):
        """Dedicated to the process launch step of the brick."""
        super(SynthStrip, self).run_process_mia()

        # default input
        fconf = capsul.engine.configurations.get(
            'capsul.engine.module.freesurfer')
        model_path = os.path.join(os.path.dirname(fconf['setup']), 'models',
                                  'synthstrip.1.pt')

        cmd = ['mri_synthstrip', '-i', self.in_file, '-o',
               self.out_file, '-m', self.out_mask, '-b', str(self.border_mm),
               '--model', model_path]

        if self.no_csf:
            cmd += ['--no-csf']

        return freesurfer.freesurfer_call(cmd)
