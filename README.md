## Decoding_Cortical_Chronotopy
This repository contains all the code used for the analyses in the paper "Decoding Cortical Chronotopy - Comparing the Influence of Different Cortical Organizational Schemes".

Additionally, you will fin one examplatory participant, to run test the code on as well as the final datasets the models for the rest-state timescales and the inter-subject 
correlation analyses were ran on. 


The requirements for running this code include:
- MATLAB, R, FSL, FreeSurfer
- SPM12
- Tools for NifTI and ANALYZE Image Jimmy Shen (2024). Tools for NIfTI and ANALYZE image (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image), MATLAB Central File Exchange. Abgerufen 9. Juli 2024. 
- fdr_bh.m  David Groppe (2024). fdr_bh (https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh), MATLAB Central File Exchange. Abgerufen 9. Juli 2024. 
- matlab_FS (the matlab repository from your Freesurfer Installation, comes with the installation)
- the git repository https://github.com/ryraut/lag-code by Ryan Raut
- Brain Connectivity Toolbox (BCT; https://sites.google.com/site/bctnet)
- Connectivity Analysis TOolbox (CATO; http://dutchconnectomelab.nl/CATO/)
- WFU pickatlas Toolbox for SPM12 (https://www.nitrc.org/projects/wfu_pickatlas/)
- MarsBaR Toolbox for SPM12 (https://marsbar-toolbox.github.io)

# ----------------------------------------------------------------------------------------------------------#

### Preparation

000_Experiment:
Contains all the inforamtion necessary to create the digit sequences and the code to run the experiment for the fMRI-Task Day (Day 1) and the DWI-Rest Day (Day 2). 
Instructions are mainly in German though an English Version of the experiment/task is included.

000_Prep_DTI:
Contains all code to reconstruct the structural white matter network from the DWI data.
This also includes the Quality Control of the parcellation and the reconstructed networks.

# ----------------------------------------------------------------------------------------------------------#

### Organizational Schemes

01_RC (Rich Club architecture)
The scripts in this Folder are used to idenfy the set of RC hub nodes in individual (binarized) connectivity Matrices.
Run the two scripts in order (Step1 the Step2).
variation_of_inforamtion.m is necessary for finding coherent module structure


02_UM (Uni-to-multimodal gradient)
Contains all code for both methods in itentifying the uni-to-multimodal Gradient
./01_SFC
	Folder contains al necessary scripts for the step-wise fucntional connectivity analysis 
	("ExampleCentralityandStepwiseConnectivity", by Jorge Sepulcre)
	Additionally, we provide the seed masks we used for ("Size_8_mm_adjusted")
	missings.mat is necessary to make the code work, to handle voxels with missing information
	
	Step1_SFC_voxelwise.m: contains all steps from preprocessed (nuisance regressed) resting-state data to SFC 
		Maps used for the Analysis
	- preparing data
	- segementation
	- normalization
	- voxelwise connectivity Matrix
	- ROIs of seeds
	- calculation of the step-wise functional connectivty
	- comparision of Maps
	- Export indivdual nodal SFC based on parcellation scheme
./02_Ji_atlas
	Contains all information and templates of the Glasser Atlas in MNI space 
	("Glasser_plus.csv", "MNI_2009c_asym") as well as the parcellation of the same MNI template accoring to the
	lausanne250 parcellation ("msbp_anat") using the multiscale brain parcellator 
	(https://multiscalebrainparcellator.readthedocs.io/en/latest/)
	Tourbier S, Aleman-Gomez Y, Griffa A, Bach Cuadra M, Hagmann P (2019, October 22) 
		sebastientourbier/multiscalebrainparcellator: Multi-Scale Brain Parcellator (Version v1.1.0). Zenodo. http://doi.org/10.5281/zenodo.2536778
	Step1_Ji_Glasser_lausanne_Overlap.m then determines the parcelwise percentatg of "multimodality"/overlap between a lausanne paracel with a Glasser parcel lavelled as multimodal
	--> in MNI space, used fro entire Group

03_PL_AM (Posterior/lateral-to-anterior/medial gradient)
Step1_annot2labels.sh: extracts the label files of each parcel of each subject, based on the atlas specific Annotation file
Step2_label2coord.R: extracts the x-, y-, and z-coordintes from the label files for each individual parcel, den coordiantes are in MNI305/fsaverage space
--> (absolute) x and y coordinates can be used to model the lateral-to-medial and anterior-to-postrior trends respectively

04_DC (Diverse Club Architecture)
Step1_DCAnalysis_SingleSub.m: Determines indivudal DC hubs (same process like for RC hubs, but in a single script)
requirements: 
all_club_bu.m: Adaption of rc_club_bu.m that determines the clubness in a binary undirected graph based on any nodal grpah measure not just Degree
consensual_partition.m: accumulates all code to determine the best resultion of modularity in one function that results in the best gamma and network modularity

05_GM (Graph measures of centrality)
Degree, betweenness centrality and closeness centrlity are extracted during RC identification
Step2_ExtractDC_PartCoef_WithMod.m: extracts the participation coefficient and within module degree z-score based on the individually optimized module structure from the DC identification

06_CT (Cortical thickness (as an approximation of cytoarchtectural differences))
extractStats.R: uses .../Freesurfer/stats to extract the individual stats of each parcel of all used parcealltion schemes

# ----------------------------------------------------------------------------------------------------------#

### Dependend Variables and Analyses
DV00_Spintest
01_PepData4SpinTest.R: Aggregate and Export the organizational scheme data, resting-state timescales data and ISC values across participants and Exports them in the correct form for the spin test.
02_SpinTest.py: Reads in all the data sets and computes spin Tests between all different organizational schemes and Output variables (Timescale and ISC) averaged across the Group. It threby checks if the varible is numeric or a label and uses Pearson correaltion, Kruskal-Wallis or χ²-test, respectively. Set up from 100,000 permutations. Saves the null distribution and all the correlation/effect size matrices.
Step3_correctResultMatrixEF.m: loads results matrix from spin test and uses different possible correction Methods to correct for multiple compatision. 
04_MeanCorrMat_plotting.R: Creates the correaltion Figure (Figure 3) from the Paper from the corrected correlation/effect size Matrix.

DV01_RS_Timescales
Contains all code from the "raw" resting-state fMRI nifti files over the preprocessing up to the analysis, including the calcualtion of all hypothesis tests and creation of graphics (consult separate README for more information).


DV02_ISC
Contains all code from the "raw" task fMRI nifti files over the preprocessing and determining ISC values up to the analysis, including the calcualtion of all hypothesis tests and creation of graphics (consult separate README for more information).





