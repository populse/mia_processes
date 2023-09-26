# -*- coding: utf-8 -*-

"""The pipelines library of the mia_processes package.


:Contains:
    :Class:
        - Dwi_whole_brain_tractography

"""

##########################################################################
# mia_processes - Copyright (C) IRMaGe/CEA, 2018
# Distributed under the terms of the CeCILL license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html
# for details.
##########################################################################

from capsul.api import Pipeline


class Dwi_whole_brain_tractography(Pipeline):
    """
    *DWI whole brain tractography*

    Please, see the complete documentation for the
    `Dwi_whole_brain_tractography pipeline in
    the populse.mia_processes website:
    <https://populse.github.io/mia_processes/html/documentation/pipelines/DwiTractography/Dwi_whole_brain_tractography.html>`_
    """

    def pipeline_definition(self):
        """Building the pipeline"""
        # nodes
        self.add_process(
            "dwi_preprocessing",
            "mia_processes.pipelines.preprocess."
            "dwi_preprocessing.Dwi_preprocessing",
        )
        self.nodes["dwi_preprocessing"].process.nodes_activation = {
            "Convert_dwi": True,
            "Convert_b0": True,
            "Denoising": True,
            "Unringing": True,
            "MotionDistortionCorrection": True,
            "Extractb0": True,
            "Meanb0": True,
            "Meanb0Corr": True,
            "FilesToList": True,
            "Concatenateb0": True,
            "BiasCorrection": True,
            "BrainMask": True,
        }
        self.add_process(
            "dwi_tissue_boundaries",
            "mia_processes.pipelines.preprocess."
            "dwi_tissue_boundaries.Dwi_tissue_boundaries",
        )
        self.nodes["dwi_tissue_boundaries"].process.nodes_activation = {
            "ConvertT1w": True,
            "Generate5tt": True,
            "Extractb0": True,
            "Meanb0": True,
            "CoregistrationFlirt": True,
            "Transformfslconvert": True,
            "ConvertMeanb0NIfTI": True,
            "Transform_1": True,
            "Transform_2": True,
            "Generate5ttseed": True,
        }
        self.add_process(
            "dwi_fod",
            "mia_processes.pipelines.preprocess."
            "dwi_fod_msmt_csd.Dwi_fod_msmt_csd",
        )
        self.nodes["dwi_fod"].process.nodes_activation = {
            "ResponseFunctionEstimation": True,
            "FODEstimation": True,
            "IntensityNormalization": True,
            "FilesToList": True,
            "ListToFile_1": True,
            "ListToFile_2": True,
            "ListToFile_3": True,
        }
        self.add_process(
            "tractography",
            "mia_processes.bricks.preprocess." "mrtrix.processes.Tractography",
        )
        self.nodes["tractography"].process.select = 10000000
        self.add_process(
            "editingtrack_1",
            "mia_processes.bricks.preprocess.mrtrix.processes.EditingTrack",
        )
        self.nodes["editingtrack_1"].process.suffix = "smaller_200k"
        self.nodes["editingtrack_1"].process.number = 200000
        self.add_process(
            "filteringtrack",
            "mia_processes.bricks.preprocess.mrtrix.processes.FilteringTrack",
        )
        self.add_process(
            "editingtrack_2",
            "mia_processes.bricks.preprocess.mrtrix.processes.EditingTrack",
        )
        self.nodes["editingtrack_2"].process.suffix = "superSmall_10k"
        self.nodes["editingtrack_2"].process.number = 10000
        self.add_process(
            "editingtrack_3",
            "mia_processes.bricks.preprocess.mrtrix.processes.EditingTrack",
        )
        self.nodes["editingtrack_3"].process.suffix = "superSmall_10k"
        self.nodes["editingtrack_3"].process.number = 10000

        self.add_process(
            "spherical_harmonic_extraction",
            "mia_processes.bricks.preprocess."
            "mrtrix.processes.SphericalHarmonicExtraction",
        )

        # links
        self.export_parameter("dwi_preprocessing", "in_dwi", is_optional=False)
        self.export_parameter(
            "dwi_preprocessing", "in_bvec", is_optional=False
        )
        self.export_parameter(
            "dwi_preprocessing", "in_bval", is_optional=False
        )
        self.export_parameter(
            "dwi_preprocessing", "in_dwi_pe_dir", is_optional=True
        )
        self.export_parameter(
            "dwi_preprocessing", "in_dwi_ro_time", is_optional=True
        )
        self.export_parameter(
            "dwi_preprocessing", "in_b0_reverse", is_optional=False
        )
        self.export_parameter(
            "dwi_tissue_boundaries", "in_T1w", is_optional=False
        )
        self.add_link("dwi_preprocessing.brain_mask->dwi_fod.brain_mask")
        self.add_link("dwi_preprocessing.preproc_dwi->dwi_fod.in_dwi")
        self.add_link(
            "dwi_preprocessing.preproc_dwi->dwi_tissue_boundaries.in_dwi"
        )
        self.add_link(
            "dwi_tissue_boundaries.5tt_coreg->tractography.act_image"
        )
        self.add_link(
            "dwi_tissue_boundaries.5tt_coreg->filteringtrack.act_image"
        )
        self.add_link(
            "dwi_tissue_boundaries.gmwmSeed_coreg->tractography.seed_gmwmi"
        )
        self.add_link("dwi_fod.wm_fod_norm->tractography.in_file")
        self.add_link("dwi_fod.wm_fod_norm->filteringtrack.in_fod")
        self.add_link("tractography.out_file->editingtrack_1.in_tracks")
        self.add_link("tractography.out_file->filteringtrack.in_tracks")
        self.add_link(
            "dwi_fod.wm_fod_norm->spherical_harmonic_extraction.in_SH_coeff"
        )
        self.export_parameter(
            "tractography", "out_file", "tracks_10mio", is_optional=False
        )
        self.add_link("tractography.out_file->editingtrack_3.in_tracks")
        self.export_parameter(
            "editingtrack_1",
            "tracks_out",
            "tracks_out_200k",
            is_optional=False,
        )
        self.export_parameter(
            "filteringtrack", "tracks_out", "tracks_sift", is_optional=False
        )
        self.add_link("filteringtrack.tracks_out->editingtrack_2.in_tracks")
        self.export_parameter(
            "editingtrack_2",
            "tracks_out",
            "tracks_sift_10k",
            is_optional=False,
        )
        self.export_parameter(
            "editingtrack_3", "tracks_out", "tracks_out_10k", is_optional=False
        )
        self.export_parameter(
            "spherical_harmonic_extraction",
            "output_image",
            "sh_peaks",
            is_optional=False,
        )

        # parameters order

        self.reorder_traits(
            (
                "in_dwi",
                "in_bvec",
                "in_bval",
                "in_dwi_pe_dir",
                "in_dwi_ro_time",
                "in_b0_reverse",
                "in_T1w",
                "sh_peaks",
                "tracks_10mio",
                "tracks_sift",
                "tracks_sift_10k",
                "tracks_out_200k",
                "tracks_out_10k",
            )
        )

        # nodes positions
        self.node_position = {
            "dwi_preprocessing": (-1205.5599999999997, -55.80800000000008),
            "dwi_tissue_boundaries": (-921.1039999999998, 261.61599999999993),
            "dwi_fod": (-866.376, -162.76),
            "tractography": (-577.1279999999999, -31.440000000000055),
            "inputs": (-1408.8198216389594, 91.07199999999989),
            "editingtrack_1": (69.12, -205.63199999999995),
            "filteringtrack": (-32.83200000000001, 513.2159999999998),
            "outputs": (1326.8182808243942, 536.9664559310762),
            "spherical_harmonic_extraction": (
                543.8188943836914,
                486.26420974153626,
            ),
            "editingtrack_2": (397.4399999999997, 751.6799999999998),
            "editingtrack_3": (461.37599999999986, -188.35199999999992),
        }

        # nodes dimensions
        self.node_dimension = {
            "dwi_preprocessing": (244.46875, 250.0),
            "dwi_tissue_boundaries": (214.59375, 145.0),
            "dwi_fod": (270.28125, 215.0),
            "tractography": (412.96875, 1265.0),
            "inputs": (130.421875, 285.0),
            "editingtrack_1": (298.4375, 600.0),
            "filteringtrack": (318.765625, 600.0),
            "outputs": (142.578125, 285.0),
            "editingtrack_2": (298.4375, 600.0),
            "editingtrack_3": (301.4375, 600.0),
            "spherical_harmonic_extraction": (220.59375, 355.0),
        }

        self.do_autoexport_nodes_parameters = False
