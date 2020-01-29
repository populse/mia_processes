# -*- coding: utf-8 -*- #

"""The pipelines library of the mia_processes package.

The objective of this module is to propose pipelines built from atomic bricks
or other pipelines provided in the mia_processes library.

:Contains:
    :Class:
        - Spatial_preprocessing_1

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

# capsul import
from capsul.api import Pipeline


class Spatial_preprocessing_1(Pipeline):
    """
 - Spatial_preprocessing_1 (mia_processes.studies.spatial_preprocessing_1.Spatial_preprocessing_1)
*** Data preprocessing for cerebrovascular reserve analyse ***
    Main pipeline modules:
        # Anatomic Images: NewSegment -> Normalize12 
        # Functional Images: Realign -> Coregister (to Anat Images) -> Normalize12 -> Smooth
    Inputs:
        # Anatomic Images (ex. 3D T1 sequence sush as T1 turbo field echo): <anat_file>
            <ex. /home/ArthurBlair/data/raw_data/Anat.nii>
        # Functional Images under hypercapnic challenge (ex. 3D T2* sequence sush as echo planar imaging): <func_files>
            <ex. ['/home/ArthurBlair/data/raw_data/Func.nii']
        # The voxel sizes (x, y & z, in mm) of the written normalised Functional Images: <voxel_sizes_func>
            <ex. [2, 2, 2]> Depends of the study; this example is valid for cevastoc, cevastoc32, etc.
    Outputs:
        # Deformation field from the NewSegment module and the Anatomic Images: <forward_deformation_field>
            <ex. /home/ArthurBlair/data/raw_data/y_Anat.nii>
        # Native space probaility maps from Anatomic Images (Grey matter, White matter, etc.): <native_class_images>
            <ex. [['/home/ArthurBlair/data/raw_data/c1Anat.nii'],
                 ['/home/ArthurBlair/data/raw_data/c2Anat.nii'],
                 ['/home/ArthurBlair/data/raw_data/c3Anat.nii'],
                 ['/home/ArthurBlair/data/raw_data/c4Anat.nii'],
                 ['/home/ArthurBlair/data/raw_data/c5Anat.nii']]>
        # Estimated translation and rotation parameters from Functional Images: <realignment_parameters>
            <ex. /home/ArthurBlair/data/raw_data/rp_Func.txt>
        # Anatomic Images normalised to the MNI space: <normalized_anat>
            <ex. /home/ArthurBlair/data/raw_data/wAnat.nii>
        # Functional Images, realigned, coregistered to Anat Images, Normalised to the MNI Space and smoothed: <smoothed_func>
            <ex. /home/ArthurBlair/data/raw_data/swrFunc.nii>
        # Mean Functional Images coregistered to the anat Images: <coregistered_source>
            <ex. /home/ArthurBlair/data/raw_data/meanFunc.nii>
    """

    def pipeline_definition(self):
        # nodes
        self.add_process("newsegment1", "mia_processes.preprocess.spm.spatial_preprocessing.NewSegment")
        self.add_process("realign1", "mia_processes.preprocess.spm.spatial_preprocessing.Realign")
        self.add_process("list_duplicate1", "mia_processes.tools.tools.List_Duplicate")
        self.add_process("normalize1", "mia_processes.preprocess.spm.spatial_preprocessing.Normalize12")
        self.add_process("normalize2", "mia_processes.preprocess.spm.spatial_preprocessing.Normalize12")
        self.add_process("smooth1", "mia_processes.preprocess.spm.spatial_preprocessing.Smooth")
        self.add_process("coregister1", "mia_processes.preprocess.spm.spatial_preprocessing.Coregister")

        # links
        self.export_parameter("list_duplicate1", "file_name", "anat_file")
        self.export_parameter("realign1", "in_files", "func_files")
        self.export_parameter("normalize2", "write_voxel_sizes", "voxel_sizes_func")
        self.export_parameter("newsegment1", "forward_deformation_field")
        self.add_link("newsegment1.forward_deformation_field->normalize1.deformation_file")
        self.add_link("newsegment1.forward_deformation_field->normalize2.deformation_file")
        self.export_parameter("newsegment1", "bias_field_images")
        self.export_parameter("newsegment1", "bias_corrected_images")
        self.export_parameter("newsegment1", "native_class_images")
        self.add_link("realign1.realigned_files->coregister1.apply_to_files")
        self.add_link("realign1.mean_image->coregister1.source")
        self.export_parameter("realign1", "realignment_parameters")
        self.add_link("list_duplicate1.out_file->coregister1.target")
        self.add_link("list_duplicate1.out_list->newsegment1.channel_files")
        self.add_link("list_duplicate1.out_list->normalize1.apply_to_files")
        self.export_parameter("normalize1", "normalized_files", "normalized_anat")
        self.add_link("normalize2.normalized_files->smooth1.in_files")
        self.export_parameter("smooth1", "smoothed_files", "smoothed_func")
        self.add_link("coregister1.coregistered_files->normalize2.apply_to_files")
        self.export_parameter("coregister1", "coregistered_source")

        # nodes positions
        self.node_position = {
            "outputs": (286.6713, 405.9444000000001),
            "inputs": (-747.6544000000002, 676.9191999999996),
            "normalize1": (-105.31994999999998, 361.1199999999999),
            "list_duplicate1": (-698.5690125000001, 322.7524000000001),
            "normalize2": (11.760050000000234, 814.1444),
            "newsegment1": (-529.162, 239.17460000000008),
            "realign1": (-572.1639500000001, 611.8702999999997),
            "coregister1": (-254.54776249999983, 717.4246000000004),
            "smooth1": (293.10992500000015, 837.1643999999998),
        }

        self.do_autoexport_nodes_parameters = False
