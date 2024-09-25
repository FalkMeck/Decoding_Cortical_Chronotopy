# Decoding Cortical Chronotopy
This repository contains all the code used for the analyses in the paper _Decoding Cortical Chronotopy - Comparing the Influence of Different Cortical Organizational Schemes_ (Mecklenbrauck, Sepulcre, Fehring, & Schubotz, in preparation). This also includes the final datasets the model estimation and comparision was run on.

Additionally, you can donload one examplary, minimally processed participant, to run test the code [here](https://uni-muenster.sciebo.de/s/zCNXmi0KUONf1Pz).

The requirements for running this code include:
- MATLAB, R, FSL, FreeSurfer
- SPM12
- [Tools for NifTI and ANALYZE Image](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) by Jimmy Shen
- [fdr_bh.m](https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh) by David Groppe
- matlab_FS (the matlab repository from your Freesurfer Installation, comes with the installation)
- [Lag-code](https://github.com/ryraut/lag-code) by Ryan Raut
- Brain Connectivity Toolbox ([BCT](https://sites.google.com/site/bctnet))
- Connectivity Analysis TOolbox ([CATO](http://dutchconnectomelab.nl/CATO/))
- [WFU pickatlas Toolbox](https://www.nitrc.org/projects/wfu_pickatlas/) for SPM12
- [MarsBaR Toolbox](https://marsbar-toolbox.github.io) for SPM12 

## Preparation

### 000_Experiment:
This folder contains all the inforamtion necessary to create the digit sequences and the code to run the experiment for the fMRI-Task Day (Day 1) and the DWI-Rest Day (Day 2). 
Instructions are mainly in German though an English Version of the experiment/task is included.

### 000_Prep_DTI:
Contains all code to reconstruct the structural white matter network from the DWI data.
This also includes the Quality Control of the parcellation and the reconstructed networks.



## Organizational Schemes

### 01_RC (Rich Club architecture)
- Requirements 
	- *variation_of_inforamtion.m* and *consensual_partition.m* are necessary for finding coherent module structure
- The scripts in this Folder are used to idenfy the set of RC hub nodes in individual (binarized) connectivity Matrices.
- Run the two scripts in order (Step1 the Step2).



### 02_UM (Uni-to-multimodal gradient)
Contains all code for both methods in itentifying the uni-to-multimodal Gradient
#### ./01_SFC
Folder contains al necessary scripts for the step-wise functional connectivity (SFC) analysis
- Requirements
	- **./ExampleCentralityandStepwiseConnectivity** (by Jorge Sepulcre): needs to be added to MATLAB paths
	- **./Size_8_mm_adjusted**: provides the masks used as seeds for SFC analysis
	- *missings.mat* is necessary to handle voxels with missing information
- *Step1_SFC_voxelwise.m*: contains all steps from preprocessed (nuisance regressed) resting-state data to SFC Maps used for the Analysis
	- preparing data
	- segementation
	- normalization
	- voxelwise connectivity Matrix
	- ROIs of seeds
	- calculation of the step-wise functional connectivty
	- comparision of Maps
	- Export indivdual nodal SFC based on parcellation scheme

#### ./02_Ji_atlas
Contains all information and templates of the Glasser Atlas in MNI space (*Glasser_plus.csv*, **./MNI_2009c_asym**) as well as the parcellation of the same MNI template accoring to the
	lausanne250 parcellation (**./msbp_anat**) using the [Multi-Scale Brain Parcellator](https://multiscalebrainparcellator.readthedocs.io/en/latest/)
- *Step1_Ji_Glasser_lausanne_Overlap.m* determines the parcelwise percentage of "multimodality"/overlap between a lausanne paracel with a Glasser parcel labelled as multimodal
	-  in MNI space, used fro entire Group

#### 03_PL_AM (Posterior/lateral-to-anterior/medial gradient)
- *Step1_annot2labels.sh*: extracts the label files of each parcel of each subject, based on the atlas specific Annotation file
- *Step2_label2coord.R*: extracts the x-, y-, and z-coordintes from the label files for each individual parcel, den coordiantes are in MNI305/fsaverage space
	- (absolute) x and y coordinates can be used to model the lateral-to-medial and anterior-to-postrior trends respectively

#### 04_DC (Diverse Club Architecture)
- Requirements
	- *variation_of_inforamtion.m* and *consensual_partition.m* are necessary for finding coherent module structure
	- *all_club_bu.m*: Adaption of *rc_club_bu.m* (in BCT) that determines the clubness in a binary undirected graph based on any nodal graph measure not just degree
	- *consensual_partition.m*: accumulates all code to determine the best resultion of modularity in one function that results in the best gamma and network modularity
- *Step1_DCAnalysis_SingleSub.m*: Determines indivudal DC hubs (same process like for RC hubs, but in a single script)


#### 05_GM (Graph measures of centrality)
- degree, betweenness centrality and closeness centrlity are extracted during RC identification
- *Step2_ExtractDC_PartCoef_WithMod.m*: extracts the participation coefficient and within module degree z-score based on the individually optimized module structure from the DC identification

#### 06_CT (Cortical thickness (as an approximation of cytoarchtectural differences))
- *extractStats.R*: uses **.../Freesurfer/stats** to extract the individual stats of each parcel of all used parcealltion schemes


## Dependend Variables and Analyses

#### DV00_Spintest
- *01_PepData4SpinTest.R*: Aggregates and exports the organizational scheme data, resting-state timescales data and ISC values across participants and exports them in the correct form for the spin test.
- *02_SpinTest.py* (adapted from Jana Fehring): Reads in all the data sets and computes spin Tests between all different organizational schemes and Output variables (Timescale and ISC) averaged across the Group. It threby checks if the varible is numeric or a label and uses Pearson correaltion, Kruskal-Wallis or χ²-test, respectively. Set up from 100,000 permutations. Saves the null distribution and all the correlation/effect size matrices.
- *Step3_correctResultMatrixEF.m*: loads results matrix from spin test and uses different possible correction Methods to correct for multiple compatision. 
- *04_MeanCorrMat_plotting.R*: Creates the correaltion Figure (Figure 3) from the Paper from the corrected correlation/effect size Matrix.

#### DV01_RS_Timescales
Contains all code from the "raw" resting-state fMRI nifti files over the preprocessing up to the analysis, including the calcualtion of all hypothesis tests and creation of graphics (consult separate **README** for more information).


#### DV02_ISC
Contains all code from the "raw" task fMRI nifti files over the preprocessing and determining ISC values up to the analysis, including the calcualtion of all hypothesis tests and creation of graphics (consult separate **README** for more information).


#### Supplementary_Analysis
Contains code that compares the two uni-to-multimodal gradient approaches from the manuscript (see folder 02_UM) with the established gradients from [Margulies et al. (2016)]https://doi.org/10.1073/pnas.1608282113 and [Sydnor et al. (2021)]https://doi.org/10.1016/j.neuron.2021.06.016
This analysis validates the methods applied in the investigation of cortical organizational schemes.