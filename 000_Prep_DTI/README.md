## 000_Prep_DTI

This folder contains all the code used to get from the our raw DWI (see example particiant) to the individual connectivity matrices.

DTI_dcm2nii.m
-  requires MATLAB tooolbox Dicm2Nii  Xiangrui Li (2024). xiangruili/dicm2nii (https://github.com/xiangruili/dicm2nii/releases/tag/v2023.02.23), GitHub. Abgerufen 9. Juli 2024. 
- creates folder structre and transforms all necessary DICOM files to NIFTI and moves them in preparation to teh DT/GQI track Analysis
- also includes a SPM12 Reorientation step, that would have to be done manually at the first iteration for each subject

DTI_pocessing.m
- Merges the necessary files (using prepareDWI_for2107_FSL6.sh)
- moves files to prepare them for Freesurfer 
- runs Freesurfers recon-all 
- runs CATO structural processing pipeline (Structural_configuration_merged_topup_eddy_noSubcortical.conf) including:
	- strucutral preprocessng (topup, eddy, movement correction)
	- parcellation into templates (aparc, lausanne120, lausanne250, lausanne500, economo, BB50human)
	- Region properties 
	- reconstruction of diffusion (using DTi and combined DTI with GQI)
	- reconstruction of fibres
	- Determination of fibre properties
	- reconstruction of network

Outlier_detection_10kin1day.m
- based on: Van Den Heuvel, M. P., Scholtens, L. H., Van Der Burgh, H. K., Agosta, F., Alloza, C., Arango, C., ... & Lange, S. C. D. (2019). 10Kin1day: a bottom-up neuroimaging initiative. Frontiers in neurology, 10, 425. https://doi.org/10.3389/fneur.2019.00425 
- determines outliers in connectvity matrices in dataset of subjects

### QC (Quality Control
Using the Enigma protocols: https://enigma.ini.usc.edu/protocols/imaging-protocols/
- requires the code from the protocals 
- scripts to run the internal, external and outlier analyses
- QC needs to be done by visual inspection



