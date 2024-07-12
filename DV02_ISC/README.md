# Inter-subject correlation analysis
This folder contains all the code to get from the task fMRI nifti files to the finished model comparision.

- *Step1_preprocessing_fMRI_ISC.m*: The preprocessing followed the tuturial given in Nastase et al. (2019) for the preprocessing for inter-subject correlation using SPM12. It includes the follwoing steps:
	- PREPARATION
		- creating all directories and moving the necessary files to those directories
	- PREPROCESSING
		- slice timing correction to the temporally and spatially middle slice 33 (66 slices in total)
		- realignment to the mean Image
		- coregistration to the anatimical Image
		- segmentation of anatomical Image
		- normalisation of data to MNI space
		- smoothing using a 6 mm smoothing kernel
	- REGRESSION OF NUISANCE
		- extracting the white matter and cerebral spinal fluid Signal
		- model specification using WM and CSF Signal, and button presses in the task as regressors
		- model estimation (classical)
	- AFTER PROCESSING (preparation for ISC Analysis)
		- creation of a Grey Matter mask
		- organizing residuals into their conditions (Single, Triplet, Nonet, Complete) 
		- z-scaling the residuals
		- creating one 4D-nifti for the individual and synchronized Task conditions
- *Step2_isc_loo.m*: Adaption of published code (Sam Nastase, https://github.com/snastase/isc-tutorial) based on code by Christian Keysers (doesn't work with a single participant, of course)
	- Loads the 4D-nifti of the residuals per condition of a subject and correlates it with the sum of the Signal of all other subjects but the one (leave one out).
	- Then saves 2 nifti-maps containing the voxelwise correlations in raw and Fisher *z*-transformed format.
- *Step3_ExtractISCfromROI.m*: To get one ISC value per parcel for the analysis, the average value is extracted per parcel. Requires MNI space parcellation of each subject (shared), WFU pick atlas Toolbox (in SPM12) for dialation of ROI and Nifti.
	- extracts all ROIs/parcels per subject, dialates them and creates binary coregistered Maps
	- extracts the ISC values from that parcel and averages them
	- saves all ISC values per subject in seperate files per condition
- *Step4_PrepareAllData.R*: takes all the data from the organizational schemes and the extracted averages ISC values for all conditions of each participant and parcel and combines them into one data set for analysis. (requires all the info about participants, so the final dataset will be shared but not the particpant information)
- *iscR2data_Final_2024_03_13.RData*: final R data set used for the Analysis.
- *Step5_1_ISCAnalysis.R*: Loads the final data (extreme values excluded and z-scaled)
Than calculates the models for all organizational schemes and does the model comparision.
Loo model comparision requires the loo() function from brms (model_loo = loo(model)). This step was not done in this script but seperately for each model on a super computer, since the compting power of a desktop computer (32 GB RAM) was not enough.
- *Step5_2_Results_ISCanalysis.R*: On basis of the estimated models, calculates and saves all the results of the hypothesis test within each model. Also the assummption check is included in this script.
- *Step5_3_TS_Plots.R*: Creates and saves all the result plots depicted in the paper.

#### References
- Nastase, S. A., Gazzola, V., Hasson, U., & Keysers, C. (2019). Measuring shared responses across subjects using intersubject correlation. Social Cognitive and Affective Neuroscience, 14(6), 669–687. 
	https://doi.org/10.1093/scan/nsz037