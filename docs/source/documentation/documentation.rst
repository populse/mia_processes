:orphan:

.. toctree::

+-----------------------+----------------------------------------------------+
|`Home <../index.html>`_|`GitHub <https://github.com/populse/mia_processes>`_|
+-----------------------+----------------------------------------------------+

Mia_processes's documentation
=============================

- **bricks**

  - **preprocess**

    - **afni**

      - `Automask <bricks/preprocess/afni/Automask.html>`_
      - `Calc <bricks/preprocess/afni/Calc.html>`_
      - `CalcDropTRs <bricks/preprocess/afni/CalcDropTRs.html>`_
      - `Despike <bricks/preprocess/afni/Despike.html>`_
      - `FWHMx <bricks/reports/FWHMx.html>`_
      - `GCOR correlation <bricks/reports/GCOR.html>`_
      - `OutliersCount <bricks/reports/OutlierCount.html>`_
      - `QualityIndex <bricks/reports/QualityIndex.html>`_
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

      - `Denoise <bricks/preprocess/dipy/Denoise.html>`_

    - **freesurfer**

      - `Binarize <bricks/preprocess/freesurfer/Binarize.html>`_
      - `SynthStrip <bricks/preprocess/freesurfer/SynthStrip.html>`_

    - **fsl**

      - `BetSurfacesExtraction <bricks/preprocess/fsl/BetSurfacesExtraction.html>`_
      - `FastSegment <bricks/preprocess/fsl/FastSegment.html>`_
      - `Smooth <bricks/preprocess/fsl/Smooth.html>`_
    
    - **others**

      - `ApplyBiasCorrection <bricks/preprocess/others/ApplyBiasCorrection.html>`_
      - `ArtifactMask <bricks/preprocess/others/ArtifactMask.html>`_
      - `Binarize <bricks/preprocess/others/Binarize.html>`_
      - `ConformImage <bricks/preprocess/others/ConformImage.html>`_
      - `ConvROI <bricks/preprocess/others/ConvROI.html>`_
      - `Enhance <bricks/preprocess/others/Enhance.html>`_
      - `EstimateSNR <bricks/preprocess/others/EstimateSNR.html>`_
      - `GradientThreshold <bricks/preprocess/others/GradientThreshold.html>`_
      - `Harmonize <bricks/preprocess/others/Harmonize.html>`_
      - `IntensityClip <bricks/preprocess/others/IntensityClip.html>`_
      - `Mask <bricks/preprocess/others/Mask.html>`_
      - `NonSteadyDetector <bricks/preprocess/others/NonSteadyDetector.html>`_
      - `Resample <bricks/preprocess/others/Resample.html>`_
      - `RotationMask <bricks/preprocess/others/RotationMask.html>`_
      - `Sanitize <bricks/preprocess/others/Sanitize.html>`_
      - `TemplateFromTemplateFlow <bricks/preprocess/others/TemplateFromTemplateFlow.html>`_
      - `Threshold <bricks/preprocess/others/Threshold.html>`_
      - `TSNR <bricks/preprocess/others/TSNR.html>`_


    - **spm**

      - `Coregister <bricks/preprocess/spm/Coregister.html>`_
      - `GM_WM_Normalize <bricks/preprocess/spm/GM_WM_Normalize.html>`_
      - `NewSegment <bricks/preprocess/spm/NewSegment.html>`_
      - `Normalize12 <bricks/preprocess/spm/Normalize12.html>`_
      - `Realign <bricks/preprocess/spm/Realign.html>`_
      - `SliceTiming <bricks/preprocess/spm/SliceTiming.html>`_
      - `Smooth <bricks/preprocess/spm/Smooth.html>`_
 
  - **reports**

    - `AnatIQMS <bricks/reports/AnatIQMS.html>`_
    - `BoldIQMS <bricks/reports/BoldIQMS.html>`_
    - `CarpetParcellation <bricks/reports/CarpetParcellation.html>`_
    - `ComputeDVARS <bricks/reports/ComputeDVARS.html>`_
    - `FramewiseDisplacement <bricks/reports/FramewiseDisplacement.html>`_
    - `Spikes <bricks/reports/Spikes.html>`_

  - **tools**

    - `Files_To_List <bricks/tools/Files_To_List.html>`_
    - `Filter_Files_List <bricks/tools/Filter_Files_List.html>`_
    - `Input_Filter <bricks/tools/Input_Filter.html>`_
    - `List_Duplicate <bricks/tools/List_Duplicate.html>`_
    - `List_To_File <bricks/tools/List_To_File.html>`_  

- **pipelines**

  - **preprocess**

    - `Anat_airmask <pipelines/preprocess/Anat_airmask.html>`_
    - `Anat_headmask <pipelines/preprocess/Anat_headmask.html>`_
    - `Anat_mni_tpms <pipelines/preprocess/Anat_mni_tpms.html>`_
    - `Anat_skullstrip_synthstrip <pipelines/preprocess/Anat_skullstrip_synthstrip.html>`_
    - `Anat_skullstrip <pipelines/preprocess/Anat_skullstrip.html>`_
    - `Anat_spatial_norm <pipelines/preprocess/Anat_spatial_norm.html>`_
    - `Bold_hmc <pipelines/preprocess/Bold_hmc.html>`_
    - `Bold_mni_align <pipelines/preprocess/Bold_mni_align.html>`_
    - `Spatial_preprocessing_1 <pipelines/preprocess/Spatial_preprocessing_1.html>`_

  - **qualityControl**

    - `Anat_mriqc <pipelines/qualityControl/Anat_mriqc.html>`_
    - `Bold_mri_qc <pipelines/qualityControl/Bold_mriqc.html>`_
  
  - **reports**

    - `Bold_iqms <pipelines/reports/Bold_iqms.html>`_


