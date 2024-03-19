:orphan:

.. toctree::

+-----------------------+----------------------------------------------------+
|`Home <../index.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------+----------------------------------------------------+

Mia_processes's documentation
=============================

To use certain mia_processes bricks (atomic processes), `third-party softwares used in these bricks must be installed <https://populse.github.io/populse_mia/html/installation/3rd-party_installations.html>`_.

- **bricks**

  - **preprocess**

    - **afni**

      - `Automask <bricks/preprocess/afni/Automask.html>`_
      - `Calc <bricks/preprocess/afni/Calc.html>`_
      - `CalcDropTRs <bricks/preprocess/afni/CalcDropTRs.html>`_
      - `Despike <bricks/preprocess/afni/Despike.html>`_
      - `FWHMx <bricks/preprocess/afni/FWHMx.html>`_
      - `GCOR <bricks/preprocess/afni/GCOR.html>`_
      - `OutliersCount <bricks/preprocess/afni/OutlierCount.html>`_
      - `QualityIndex <bricks/preprocess/afni/QualityIndex.html>`_
      - `RefitDeoblique <bricks/preprocess/afni/RefitDeoblique.html>`_
      - `SkullStrip <bricks/preprocess/afni/SkullStrip.html>`_
      - `TShift <bricks/preprocess/afni/TShift.html>`_
      - `TStatMean <bricks/preprocess/afni/TStatMean.html>`_
      - `Volreg <bricks/preprocess/afni/Volreg.html>`_

    - **ants**

      - `AffineInitializer <bricks/preprocess/ants/AffineInitializer.html>`_
      - `ApplyTransform <bricks/preprocess/ants/ApplyTransform.html>`_
      - `N4BiasFieldCorrection <bricks/preprocess/ants/N4BiasFieldCorrection.html>`_
      - `Registration <bricks/preprocess/ants/Registration.html>`_

    - **dipy**

      - `ComputeDKI <bricks/preprocess/dipy/ComputeDKI.html>`_
      - `Denoise <bricks/preprocess/dipy/Denoise.html>`_

    - **freesurfer**

      - `Binarize <bricks/preprocess/freesurfer/Binarize.html>`__
      - `SynthStrip <bricks/preprocess/freesurfer/SynthStrip.html>`_
      - `SynthStripMriqc <bricks/preprocess/freesurfer/SynthStripMriqc.html>`_

    - **fsl**

      - `BetSurfacesExtraction <bricks/preprocess/fsl/BetSurfacesExtraction.html>`_
      - `ConvertXFM <bricks/preprocess/fsl/ConvertXFM.html>`_
      - `EpiReg <bricks/preprocess/fsl/EpiReg.html>`_
      - `ExtractROI <bricks/preprocess/fsl/ExtractROI.html>`_
      - `FastSegment <bricks/preprocess/fsl/FastSegment.html>`_
      - `Flirt <bricks/preprocess/fsl/Flirt.html>`_
      - `Smooth <bricks/preprocess/fsl/Smooth.html>`__

    - **matlab_wrap**

      - `ComputeBrainVolume <bricks/preprocess/matlab_wrap/ComputeBrainVolume.html>`__

    - **mrtrix**

      - `ConstrainedSphericalDeconvolution <bricks/preprocess/mrtrix/ConstrainedSphericalDeconvolution.html>`_
      - `DWIBiasCorrect <bricks/preprocess/mrtrix/DWIBiasCorrect.html>`_
      - `DWIBrainMask <bricks/preprocess/mrtrix/DWIBrainMask.html>`_
      - `DWIDenoise <bricks/preprocess/mrtrix/DWIDenoise.html>`_
      - `DWIExtract <bricks/preprocess/mrtrix/DWIExtract.html>`_
      - `DWIPreproc <bricks/preprocess/mrtrix/DWIPreproc.html>`_
      - `EditingTrack <bricks/preprocess/mrtrix/EditingTrack.html>`_
      - `FilteringTrack <bricks/preprocess/mrtrix/FilteringTrack.html>`_
      - `FitTensor <bricks/preprocess/mrtrix/FitTensor.html>`_
      - `Generate5tt2gmwmi <bricks/preprocess/mrtrix/Generate5tt2gmwmi.html>`_
      - `Generate5ttfsl <bricks/preprocess/mrtrix/Generate5ttfsl.html>`_
      - `MRCat <bricks/preprocess/mrtrix/MRCat.html>`_
      - `MRConvert <bricks/preprocess/mrtrix/MRConvert.html>`_
      - `MRDeGibbs <bricks/preprocess/mrtrix/MRDeGibbs.html>`_
      - `MRMath <bricks/preprocess/mrtrix/MRMath.html>`_
      - `MRTransform <bricks/preprocess/mrtrix/MRTransform.html>`_
      - `MTNormalise <bricks/preprocess/mrtrix/MTNormalise.html>`_
      - `ResponseSDDhollander <bricks/preprocess/mrtrix/ResponseSDDhollander.html>`_
      - `SphericalHarmonicExtraction <bricks/preprocess/mrtrix/SphericalHarmonicExtraction.html>`_
      - `TensorMetrics <bricks/preprocess/mrtrix/TensorMetrics.html>`_
      - `Tractography <bricks/preprocess/mrtrix/Tractography.html>`_
      - `TransformFSLConvert <bricks/preprocess/mrtrix/TransformFSLConvert.html>`_

    - **others**

      - `ApplyBiasCorrection <bricks/preprocess/others/ApplyBiasCorrection.html>`_
      - `ArtifactMask <bricks/preprocess/others/ArtifactMask.html>`_
      - `Binarize <bricks/preprocess/others/Binarize.html>`__
      - `ConformImage <bricks/preprocess/others/ConformImage.html>`_
      - `ConvROI <bricks/preprocess/others/ConvROI.html>`_
      - `Enhance <bricks/preprocess/others/Enhance.html>`_
      - `EstimateSNR <bricks/preprocess/others/EstimateSNR.html>`_
      - `ExtractROIbyLabel <bricks/preprocess/others/ExtractROIbyLabel.html>`_
      - `ExtractSignalROI <bricks/preprocess/others/ExtractSignalROI.html>`_
      - `GradientThreshold <bricks/preprocess/others/GradientThreshold.html>`_
      - `Harmonize <bricks/preprocess/others/Harmonize.html>`_
      - `IntensityClip <bricks/preprocess/others/IntensityClip.html>`_
      - `Mask <bricks/preprocess/others/Mask.html>`_
      - `NonSteadyDetector <bricks/preprocess/others/NonSteadyDetector.html>`_
      - `Resample1 <bricks/preprocess/others/Resample1.html>`_
      - `Resample2 <bricks/preprocess/others/Resample2.html>`_
      - `RotationMask <bricks/preprocess/others/RotationMask.html>`_
      - `Sanitize <bricks/preprocess/others/Sanitize.html>`_
      - `TSNR <bricks/preprocess/others/TSNR.html>`_
      - `TemplateFromTemplateFlow <bricks/preprocess/others/TemplateFromTemplateFlow.html>`_
      - `Threshold <bricks/preprocess/others/Threshold.html>`_

    - **spm**

      - `Coregister <bricks/preprocess/spm/Coregister.html>`_
      - `GM_WM_Normalize <bricks/preprocess/spm/GM_WM_Normalize.html>`_
      - `NewSegment <bricks/preprocess/spm/NewSegment.html>`_
      - `Normalize12 <bricks/preprocess/spm/Normalize12.html>`_
      - `Realign <bricks/preprocess/spm/Realign.html>`_
      - `SliceTiming <bricks/preprocess/spm/SliceTiming.html>`_
      - `Smooth <bricks/preprocess/spm/Smooth.html>`__

  - **reports**

    - `AnatIQMs <bricks/reports/AnatIQMs.html>`_
    - `BoldIQMs <bricks/reports/BoldIQMs.html>`_
    - `BoldIQMsPlot <bricks/reports/BoldIQMsPlot.html>`_
    - `CarpetParcellation <bricks/reports/CarpetParcellation.html>`_
    - `ComputeDVARS <bricks/reports/ComputeDVARS.html>`_
    - `FramewiseDisplacement <bricks/reports/FramewiseDisplacement.html>`_
    - `Mean_stdDev_calc <bricks/reports/Mean_stdDev_calc.html>`_
    - `PlotSignalROI <bricks/reports/PlotSignalROI.html>`_
    - `ReportAnatMriqc <bricks/reports/ReportAnatMriqc.html>`_
    - `ReportCO2inhalCvr <bricks/reports/ReportCO2inhalCvr.html>`_
    - `ReportFuncMriqc <bricks/reports/ReportFuncMriqc.html>`_
    - `ReportGroupMriqc <bricks/reports/ReportGroupMriqc.html>`_
    - `Result_collector <bricks/reports/Result_collector.html>`_
    - `Spikes <bricks/reports/Spikes.html>`_

  - **stats**

    - **spm**

      - `EstimateContrast <bricks/stats/spm/EstimateContrast.html>`_
      - `EstimateModel <bricks/stats/spm/EstimateModel.html>`_
      - `Level1Design <bricks/stats/spm/Level1Design.html>`_
      - `MultipleRegressionDesign <bricks/stats/spm/MultipleRegressionDesign.html>`_
      - `OneSampleTTestDesign <bricks/stats/spm/OneSampleTTestDesign.html>`_
      - `PairedTTestDesign <bricks/stats/spm/PairedTTestDesign.html>`_
      - `TwoSampleTTestDesign <bricks/stats/spm/TwoSampleTTestDesign.html>`_

  - **tools**

    - `Concat_to_list <bricks/tools/Concat_to_list.html>`_
    - `Concat_to_list_of_list <bricks/tools/Concat_to_list_of_list.html>`_
    - `Delete_Data <bricks/tools/Delete_data.html>`_
    - `Files_To_List <bricks/tools/Files_To_List.html>`_
    - `Filter_Files_List <bricks/tools/Filter_Files_List.html>`_
    - `Find_In_List <bricks/tools/Find_In_List.html>`_
    - `Get_Conditions_From_csv <bricks/tools/Get_Conditions_From_csv.html>`_
    - `Get_Patient_Name <bricks/tools/Get_Patient_Name.html>`_
    - `Import_Data <bricks/tools/Import_Data.html>`_
    - `Input_Filter <bricks/tools/Input_Filter.html>`_
    - `List_Duplicate <bricks/tools/List_Duplicate.html>`_
    - `List_To_File <bricks/tools/List_To_File.html>`_
    - `Make_A_List <bricks/tools/Make_A_List.html>`_
    - `Make_CVR_reg_physio <bricks/tools/Make_CVR_reg_physio.html>`_

- **pipelines**

  - **CerebVascularReact**

    - `CO2_inhalation <pipelines/CerebVascularReact/CO2_inhalation.html>`_

  - **DWITractography**

    - `DWI_whole_brain_tractograpy <pipelines/DWITractography/DWI_whole_brain_tractograpy.html>`_

  - **preprocess**

    - `Anat_airmask <pipelines/preprocess/Anat_airmask.html>`_
    - `Anat_headmask <pipelines/preprocess/Anat_headmask.html>`_
    - `Anat_mni_tpms <pipelines/preprocess/Anat_mni_tpms.html>`_
    - `Anat_skullstrip <pipelines/preprocess/Anat_skullstrip.html>`_
    - `Anat_skullstrip_synthstrip <pipelines/preprocess/Anat_skullstrip_synthstrip.html>`_
    - `Anat_spatial_norm <pipelines/preprocess/Anat_spatial_norm.html>`_
    - `Bold_hmc <pipelines/preprocess/Bold_hmc.html>`_
    - `Bold_mni_align <pipelines/preprocess/Bold_mni_align.html>`_
    - `Bold_spatial_preprocessing1 <pipelines/preprocess/Bold_spatial_preprocessing1.html>`_
    - `Bold_spatial_preprocessing2 <pipelines/preprocess/Bold_spatial_preprocessing2.html>`_
    - `Dwi_fod_msmt_csd <pipelines/preprocess/Dwi_fod_msmt_csd.html>`_
    - `Dwi_preprocessing <pipelines/preprocess/Dwi_preprocessing.html>`_
    - `Dwi_tissue_boundaries <pipelines/preprocess/Dwi_tissue_boundaries.html>`_
    - `Spatial_mask <pipelines/preprocess/Spatial_mask.html>`_

  - **qualityControl**

    - `Anat_mriqc <pipelines/qualityControl/Anat_mriqc.html>`_
    - `Bold_mri_qc <pipelines/qualityControl/Bold_mriqc.html>`_

  - **reports**

    - `Bold_iqms <pipelines/reports/Bold_iqms.html>`_
    - `Extract_roi_param <pipelines/reports/Extract_roi_param.html>`_

  - **stat**

    - `Bold_stat_cvr <pipelines/stat/Bold_stat_cvr.html>`_
