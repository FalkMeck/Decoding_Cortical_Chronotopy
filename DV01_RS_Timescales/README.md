# Restin-state timescale analysis
This folder contains all the code to get from the resting-state fMRI nifti files to the finished model comparision.

- Step1_Preprocessing: The preprocessing of the data is based on the Methods presented in Ito et al. (2020). The used code also was also partially taken from https://github.com/ColeLab/hierarchy2020
	- *rest_prep_cato.m*: prepares the resting state data fro later steps using CATO
	this scripts calls on *Functional_configuration_prep.conf* which includes:
		- functional preprocessing: registration to anatomy, Segmentation, slice timing, movement correction 		  (*preprocess_default*)
		- parcellation (taken from structural CATO steps)
		- Collection of Region properties
		- computation of Motion metrics
		- CATO also could nuisance regress the data but we opted for the more thorough method (Ito et al., 2020)
	- *prep_nuisance.m*: to run the nuisance regression we need specific masks that can be extracted from the parcelled structrual images and the movement parameters must be in a specific text file. Both steps are completed here.
	- *nuisanceRegressionPipeline_FM.py*: all credit to Takuya Ito, we used the 24pXaCompCorXVolterra-method with 64 regressors + spike Regression
	- *run_nuisanceRegressionPipeline_FM.py*: just runs the defined nuisance Regression
- *Step2_reconstructrion_timeSeries_network.m*: This scripts contains the remaining steps of the CATO functional processing pipeline but extracted to MATLAB. The process continues after the nuisance regession and contains bandpass filtering and the extraction of timeseries data which will be used for the following steps. Technically also functional connectivity Matrices are calculated but not further used.
- *Step3_spectral_template.m*: based on the model free calcualtion of the timescales by the autocorrealtion function in Raut et al. (2020) also requires the shared code by Ryan Raut (https://github.com/ryraut/lag-code). Uses the extracted timeseries to determine timescales by the precise abscissa of the autocorrelation function at half its height.
- Step4_CombineAllData.R: takes all the data from the organizational schemes and the resting-state timescales of each participant and parcel and combines them into one data set for Analysis.
(requires all the info about participants, so the final dataset will be shared but not the particpant information)
- *R2data_lasuanne250_2024_09_21.RData*: R dataset that was used for the analysis (Step5_1 - Step5_3)
- *Step5_1_TS_analysis.R*: Loads the final data (outlier extrated, z-scaling and logarithmic transforamtion of timescales). Than calculates the models for all organizational schemes, checks their assumptions and does the model comparision.
*Step5_2_TS_BayesResults.R*: On Basis of the estimated models, calculates and saves all the results of the hypothesis test within each model
- *Step5_3_TS_Plots.R*: Creates and saves all the result plots depicted in the paper.


#### References
- Ito, T., Hearne, L. J., & Cole, M. W. (2020). A cortical hierarchy of localized and distributed processes revealed 
	via dissociation of task activations, connectivity changes, and intrinsic timescales. NeuroImage, 221, 117141. https://doi.org/10.1016/j.neuroimage.2020.117141
- Raut, R. V., Snyder, A. Z., & Raichle, M. E. (2020). Hierarchical dynamics as a macroscopic organizing principle of the human brain. Proceedings of the National Academy of Sciences, 117(34), 20890-20897. 
	https://doi.org/10.1073/pnas.2003383117
