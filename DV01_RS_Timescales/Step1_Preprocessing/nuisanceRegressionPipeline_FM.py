# Takuya Ito
# 09/11/2018

# Post-processing nuisance regression using Ciric et al. 2017 inspired best-practices
## OVERVIEW
# There are two main parts to this script/set of functions
# 1. "step1_createNuisanceRegressors" 
#   Generates a variety of nuisance regressors, such as motionSpikes, aCompCor regressors, etc. that are essential to a subset of Ciric-style models, with the addition of some new combinations (e.g., aCompCor + spikeReg + movement parameters) 
#   This is actually the bulk of the script, and takes quite a while to compute, largely due to the fact that we need to load in 4D time series from the raw fMRI data (in order to compute regressors such as global signal)
# 2. "step2_nuisanceRegression"
#   This is the function that actually performs the nuisance regression, using regressors obtained from step1. There are a variety of models to choose from, including:
#       The best model from Ciric et al. (2017) (e.g., 36p + spikeReg)
#       What I call the "legacy Cole Lab models", which are the traditional 6 motion parameters, gsr, wm and ventricle time series and all their derivatives (e.g., 18p)
#       There is also 16pNoGSR, which is the above, but without gsr and its derivative.
#   Ultimately, read below for other combinations; what I would consider the best option that does NOT include GSR is the default, called "24pXaCompCorXVolterra" - read below for what it entails...

# IMPORTANT: In general, only functions step1, step2 and the parameters preceding that will need to be edited. There are many helper functions below, but in theory, they should not be edited.
# Currently, this script is defaulted to create the nuisance regressors in your current working directory (in a sub directory), and the glm output in your current working directory
# The default is set to use data from the HCP352 QC'd data set, so will need to be updated accordingly. 
# For now, this only includes extensive nuisance regression. Any task regression will need to be performed independently after this.

## EXAMPLE USAGE:
# import nuisanceRegressionPipeline as nrp
# nrp.step1_createNuisanceRegressors(nproc=8)
# nrp.step2_nuisanceRegression(nproc=5, model='24pXaCompCorXVolterra',spikeReg=False,zscore=False)

## DISCLAIMER: This is a first draft, so... keep that in mind.

import numpy as np
import os
#import glob
#from nipy.modalities.fmri.hemodynamic_models import spm_hrf
#import multiprocessing as mp
import statsmodels.api as sm
import h5py
import scipy.stats as stats
from scipy import signal
import nibabel as nib
import scipy
#import time
import warnings
import gzip
warnings.simplefilter('ignore', np.ComplexWarning)

## Define GLOBAL variables (variables accessible to all functions
# Define base data directory
#datadir = '.../RICHIE2_Analysis'
# Define number of frames to skip
framesToSkip = 0
# Define all runs you want to preprocess
allRuns = ['_fmri']
outputdir = '.../Analysis/' 
# Define the *output* directory for nuisance regressors
nuis_reg_dir = outputdir +'nuisanceRegressors/'
# Create directory if it doesn't exist
if not os.path.exists(nuis_reg_dir): os.makedirs(nuis_reg_dir)
# Define the *output* directory for preprocessed data

#
# Define subjects list
subjNums = ['HFG_121']

def step1_createNuisanceRegressors():
    """
    Function to generate subject-wise nuisance parameters in parallel
    This function first defines a local function (a function within this function) to generate each subject's nuisance regressors
    Then we use the multiprocessing module to generate regressors for multiple subjects at a time

    **Note: Parameters in this function may need to be edited for project-specific purposes. Sections in which editing should NOT be done are noted
    """
    # Make below function global, so it is accessible to the parallel process (don't change this)
    global _createNuisanceRegressorsSubject
    def _createNuisanceRegressorsSubject(subj):
        ## Potentially will need to be edited, according to project
        # Directory for all the masks
        analysisdir = '.../Analysis/'
        maskdir = analysisdir + subj + '/masks/'
        # Path and file name for the whole-brain mask
        globalmask = maskdir + subj + '_wholebrainmask.nii.gz'
        # Path and file name for the white matter mask
        wmmask = maskdir + subj + '_wmMask.nii.gz'
        # Path and file name for the ventricle mask
        ventriclesmask = maskdir + subj + '_ventricles.nii.gz'
        # This is the path and filename for the output regressors
        nuisance_reg_filename = nuis_reg_dir + subj + '_nuisanceRegressors.h5'
        # Define the directory containing the raw preprocessed data
        datadir = analysisdir + subj + '/minimally_Processed/resting_state_fMRI/'
        # Number of principal components to extract out of WM and ventricle signals
        compCorComponents = 5
        # Spike regression threshold, using relative root-mean-square displacement (in mm)
        spikeReg = .25
        ####

        for run in allRuns:
        
            print('creating nuisance regressors for subject', subj, 'run:', run)
            # This is the fMRI 4d file (volumetric) to obtain the noise signals -- done for each run
            inputname = datadir + subj + run + '.nii.gz' 

            #### Obtain movement parameters -- this will differ across preprocessing pipelines (e.g., HCP vs. typical)
            # For all 12 movement parameters (6 regressors + derivatives)
            movementRegressors = np.loadtxt(analysisdir + subj + '/Movement_Regressors.txt')
            # Separate the two parameters out for clarity
            # x, y, z + 3 rotational movements
            motionParams = movementRegressors[:,:6]
            # The derivatives of the above movements (backwards differentiated)
            motionParams_deriv = movementRegressors[:,6:] # HCP automatically computes derivative of motion parameters
            
            ####
            # DO NOT CHANGE THIS SECTION, IT IS NECESSARY FOR THE SCRIPT TO RUN
            h5f = h5py.File(nuisance_reg_filename,'a')
            try:
                h5f.create_dataset(run + '/motionParams',data=motionParams)
                h5f.create_dataset(run + '/motionParams_deriv',data=motionParams_deriv)
            except:
                del h5f[run + '/motionParams'], h5f[run + '/motionParams_deriv']
                h5f.create_dataset(run + '/motionParams',data=motionParams)
                h5f.create_dataset(run + '/motionParams_deriv',data=motionParams_deriv)
            h5f.close()
            # END OF DO NOT CHANGE
            ####
    
            #### Obtain relative root-mean-square displacement -- this will differ across preprocessing pipelines
            # A simple alternative is to compute the np.sqrt(x**2 + y**2 + z**2), where x, y, and z are motion displacement parameters
            # e.g., x = x[t] - x[t-1]; y = y[t] - y[t-1]; z = z[t] - z[t-1]
            # For HCP data, just load in the relative RMS
            relativeRMS = np.loadtxt(analysisdir + subj + '/Movement_RelativeRMS.txt')
            # Calculate motion spike regressors using helper functions defined below
            _createMotionSpikeRegressors(relativeRMS, subj, run, spikeReg=spikeReg)
            # Extract physiological noise signals using helper functions defined below
            _createPhysiologicalNuisanceRegressors(inputname, subj, run, globalmask, wmmask, ventriclesmask, aCompCor=compCorComponents)

    # Construct parallel processes to run the local function in parallel (subject-wise parallelization)
    # Outputs will be found in "nuis_reg_dir" parameter
    for subj in subjNums:
        _createNuisanceRegressorsSubject(subj)
            
            #pool = mp.Pool(processes=nproc)
            #pool.map_async(_createNuisanceRegressorsSubject,subjNums).get()
            #pool.close()
            #pool.join()
    

def step2_nuisanceRegression(model='24pXaCompCorXVolterra',spikeReg=True,zscore=False):
    """
    Function to perform nuisance regression on each run separately
    This uses parallel processing, but parallelization occurs within each subject
    Each subject runs regression on each region/voxel in parallel, thus iterating subjects and runs serially
    
    Input parameters:
        subj    : subject number as a string
        run     : task run
        outputdir: Directory for GLM output, as an h5 file (each run will be contained within each h5)
        model   : model choices for linear regression. Models include:
                    1. 24pXaCompCorXVolterra [default]
                        Variant from Ciric et al. 2017. 
                        Includes (64 regressors total):
                            - Movement parameters (6 directions; x, y, z displacement, and 3 rotations) and their derivatives, and their quadratics (24 regressors)
                            - aCompCor (5 white matter and 5 ventricle components) and their derivatives, and their quadratics (40 regressors)
                    2. 18p (the lab's legacy default)
                        Includes (18 regressors total):
                            - Movement parameters (6 directions) and their derivatives (12 regressors)
                            - Global signal and its derivative (2 regressors)
                            - White matter signal and its derivative (2 regressors)
                            - Ventricles signal and its derivative (2 regressors)
                    3. 16pNoGSR (the legacy default, without GSR)
                        Includes (16 regressors total):
                            - Movement parameters (6 directions) and their derivatives (12 regressors)
                            - White matter signal and its derivative (2 regressors)
                            - Ventricles signal and its derivative (2 regressors)
                    4. 12pXaCompCor (Typical motion regression, but using CompCor (noGSR))
                        Includes (32 regressors total):
                            - Movement parameters (6 directions) and their derivatives (12 regressors)
                            - aCompCor (5 white matter and 5 ventricle components) and their derivatives (no quadratics; 20 regressors)
                    5. 36p (State-of-the-art, according to Ciric et al. 2017)
                        Includes (36 regressors total - same as legacy, but with quadratics):
                            - Movement parameters (6 directions) and their derivatives and quadratics (24 regressors)
                            - Global signal and its derivative and both quadratics (4 regressors)
                            - White matter signal and its derivative and both quadratics (4 regressors)
                            - Ventricles signal and its derivative (4 regressors)
        spikeReg : spike regression (Satterthwaite et al. 2013) [True/False]
                        Note, inclusion of this will add additional set of regressors, which is custom for each subject/run
        zscore   : Normalize data (across time) prior to fitting regression
        nproc = number of processes to use via multiprocessing
    """
    # Iterate through each subject
    for subj in subjNums:
        # Iterate through each run
        for run in allRuns:
            print('Running regression on subject', subj, '| run', run)
            print('\tModel:', model, 'with spikeReg:', spikeReg, '| zscore:=', zscore)
            ## Load in data to be preprocessed - This needs to be a space x time 2d array
            analysisdir = '.../Analysis/'
            datadir = analysisdir + subj + '/minimally_Processed/resting_state_fMRI/'
            
            niiInput = datadir + subj + run + '.nii.gz' 
            fMRI4d = np.asarray(nib.load(niiInput).get_fdata())
               
            # convert data from 4D to 2D array
            data = fMRI4d.reshape((fMRI4d.shape[0]*fMRI4d.shape[1]*fMRI4d.shape[2]),fMRI4d.shape[3])
            
            #inputfile = datadir + '/hcpPreprocessedData/' + subj + '_GlasserParcellated_' + run + '_LR.csv'
            # Load data
            # data = np.loadtxt(inputfile,delimiter=',')
            # Run nuisance regression for this subject's run, using a helper function defined below
            # Data will be output in 'outputdir', defined above
            _nuisanceRegression(subj, run, data, outputdir, model=model,spikeReg=spikeReg,zscore=zscore)
            
            


#########################################
# Functions that probably don't need to be edited

def _nuisanceRegression(subj, run, inputdata, outputdir,  model='24pXaCompCorXVolterra', spikeReg=True, zscore=False):
    """
    This function runs nuisance regression on the Glasser Parcels (360) on a single subjects run
    Will only regress out noise parameters given the model choice (see below for model options)

    Input parameters:
        subj    : subject number as a string
        run     : task run
        outputdir: Directory for GLM output, as an h5 file (each run will be contained within each h5)
        model   : model choices for linear regression. Models include:
                    1. 24pXaCompCorXVolterra [default]
                        Variant from Ciric et al. 2017. 
                        Includes (64 regressors total):
                            - Movement parameters (6 directions; x, y, z displacement, and 3 rotations) and their derivatives, and their quadratics (24 regressors)
                            - aCompCor (5 white matter and 5 ventricle components) and their derivatives, and their quadratics (40 regressors)
                    2. 18p (the legacy default)
                        Includes (18 regressors total):
                            - Movement parameters (6 directions) and their derivatives (12 regressors)
                            - Global signal and its derivative (2 regressors)
                            - White matter signal and its derivative (2 regressors)
                            - Ventricles signal and its derivative (2 regressors)
                    3. 16pNoGSR (the legacy default, without GSR)
                        Includes (16 regressors total):
                            - Movement parameters (6 directions) and their derivatives (12 regressors)
                            - White matter signal and its derivative (2 regressors)
                            - Ventricles signal and its derivative (2 regressors)
                    4. 12pXaCompCor (Typical motion regression, but using CompCor (noGSR))
                        Includes (32 regressors total):
                            - Movement parameters (6 directions) and their derivatives (12 regressors)
                            - aCompCor (5 white matter and 5 ventricle components) and their derivatives (no quadratics; 20 regressors)
                    5. 36p (State-of-the-art, according to Ciric et al. 2017)
                        Includes (36 regressors total - same as legacy, but with quadratics):
                            - Movement parameters (6 directions) and their derivatives and quadratics (24 regressors)
                            - Global signal and its derivative and both quadratics (4 regressors)
                            - White matter signal and its derivative and both quadratics (4 regressors)
                            - Ventricles signal and its derivative (4 regressors)
        spikeReg : spike regression (Satterthwaite et al. 2013) [True/False]
                        Note, inclusion of this will add additional set of regressors, which is custom for each subject/run
        zscore   : Normalize data (across time) prior to fitting regression
        nproc = number of processes to use via multiprocessing
    """

    data = inputdata

    tMask = np.ones((data.shape[1],))
    tMask[:framesToSkip] = 0

    # Skip frames
    data = data[:,framesToSkip:]
    
    # Demean each run
    data = signal.detrend(data,axis=1,type='constant')
    # Detrend each run
    data = signal.detrend(data,axis=1,type='linear')
    tMask = np.asarray(tMask,dtype=bool)
    
    if zscore:
        data = stats.zscore(data,axis=1)

    nROIs = data.shape[0]

    # Load nuisance regressors for this data
    h5f = h5py.File(nuis_reg_dir + subj + '_nuisanceRegressors.h5','r') 
    if model=='24pXaCompCorXVolterra':
        # Motion parameters + derivatives
        motion_parameters = h5f[run]['motionParams'][:].copy()
        motion_parameters_deriv = h5f[run]['motionParams_deriv'][:].copy()
        # WM aCompCor + derivatives
        aCompCor_WM = h5f[run]['aCompCor_WM'][:].copy()
        aCompCor_WM_deriv = h5f[run]['aCompCor_WM_deriv'][:].copy()
        # Ventricles aCompCor + derivatives
        aCompCor_ventricles = h5f[run]['aCompCor_ventricles'][:].copy()
        aCompCor_ventricles_deriv = h5f[run]['aCompCor_ventricles_deriv'][:].copy()
        # Create nuisance regressors design matrix
        nuisanceRegressors = np.hstack((motion_parameters, motion_parameters_deriv, aCompCor_WM, aCompCor_WM_deriv, aCompCor_ventricles, aCompCor_ventricles_deriv))
        quadraticRegressors = nuisanceRegressors**2
        nuisanceRegressors = np.hstack((nuisanceRegressors,quadraticRegressors))
    
    elif model=='18p':
        # Motion parameters + derivatives
        motion_parameters = h5f[run]['motionParams'][:].copy()
        motion_parameters_deriv = h5f[run]['motionParams_deriv'][:].copy()
        # Global signal + derivatives
        global_signal = h5f[run]['global_signal'][:].copy()
        global_signal_deriv = h5f[run]['global_signal_deriv'][:].copy()
        # white matter signal + derivatives
        wm_signal = h5f[run]['wm_signal'][:].copy()
        wm_signal_deriv = h5f[run]['wm_signal_deriv'][:].copy()
        # ventricle signal + derivatives
        ventricle_signal = h5f[run]['ventricle_signal'][:].copy()
        ventricle_signal_deriv = h5f[run]['ventricle_signal_deriv'][:].copy()
        # Create nuisance regressors design matrix
        tmp = np.vstack((global_signal,global_signal_deriv,wm_signal,wm_signal_deriv,ventricle_signal,ventricle_signal_deriv)).T # Need to vstack, since these are 1d arrays
        nuisanceRegressors = np.hstack((motion_parameters, motion_parameters_deriv, tmp))

    elif model=='16pNoGSR':
        # Motion parameters + derivatives
        motion_parameters = h5f[run]['motionParams'][:].copy()
        motion_parameters_deriv = h5f[run]['motionParams_deriv'][:].copy()
        # white matter signal + derivatives
        wm_signal = h5f[run]['wm_signal'][:].copy()
        wm_signal_deriv = h5f[run]['wm_signal_deriv'][:].copy()
        # ventricle signal + derivatives
        ventricle_signal = h5f[run]['ventricle_signal'][:].copy()
        ventricle_signal_deriv = h5f[run]['ventricle_signal_deriv'][:].copy()
        # Create nuisance regressors design matrix
        tmp = np.vstack((wm_signal,wm_signal_deriv,ventricle_signal,ventricle_signal_deriv)).T # Need to vstack, since these are 1d arrays
        nuisanceRegressors = np.hstack((motion_parameters, motion_parameters_deriv, tmp))
    
    elif model=='12pXaCompCor':
        # Motion parameters + derivatives
        motion_parameters = h5f[run]['motionParams'][:].copy()
        motion_parameters_deriv = h5f[run]['motionParams_deriv'][:].copy()
        # WM aCompCor + derivatives
        aCompCor_WM = h5f[run]['aCompCor_WM'][:].copy()
        aCompCor_WM_deriv = h5f[run]['aCompCor_WM_deriv'][:].copy()
        # Ventricles aCompCor + derivatives
        aCompCor_ventricles = h5f[run]['aCompCor_ventricles'][:].copy()
        aCompCor_ventricles_deriv = h5f[run]['aCompCor_ventricles_deriv'][:].copy()
        # Create nuisance regressors design matrix
        nuisanceRegressors = np.hstack((motion_parameters, motion_parameters_deriv, aCompCor_WM, aCompCor_WM_deriv, aCompCor_ventricles, aCompCor_ventricles_deriv))
    
    elif model=='36p':
        # Motion parameters + derivatives
        motion_parameters = h5f[run]['motionParams'][:].copy()
        motion_parameters_deriv = h5f[run]['motionParams_deriv'][:].copy()
        # Global signal + derivatives
        global_signal = h5f[run]['global_signal'][:].copy()
        global_signal_deriv = h5f[run]['global_signal_deriv'][:].copy()
        # white matter signal + derivatives
        wm_signal = h5f[run]['wm_signal'][:].copy()
        wm_signal_deriv = h5f[run]['wm_signal_deriv'][:].copy()
        # ventricle signal + derivatives
        ventricle_signal = h5f[run]['ventricle_signal'][:].copy()
        ventricle_signal_deriv = h5f[run]['ventricle_signal_deriv'][:].copy()
        # Create nuisance regressors design matrix
        tmp = np.vstack((global_signal,global_signal_deriv,wm_signal,wm_signal_deriv,ventricle_signal,ventricle_signal_deriv)).T # Need to vstack, since these are 1d arrays
        nuisanceRegressors = np.hstack((motion_parameters, motion_parameters_deriv, tmp))
        quadraticRegressors = nuisanceRegressors**2
        nuisanceRegressors = np.hstack((nuisanceRegressors,quadraticRegressors))


    if spikeReg:
        # Obtain motion spikes
        try:
            motion_spikes = h5f[run]['motionSpikes'][:].copy()
            nuisanceRegressors = np.hstack((nuisanceRegressors,motion_spikes))
        except:
            print('Spike regression option was chosen... but no motion spikes for subj', subj, '| run', run, '!')
        # Update the model name - to keep track of different model types for output naming
        model = model + '_spikeReg' 

    if zscore:
        model = model + '_zscore'

    h5f.close()
    # Skip first 5 frames of nuisanceRegressors, too # not necessarym excluded them before frame_to_skip set to 0
    nuisanceRegressors = nuisanceRegressors[framesToSkip:,:].copy()
    
    inputs = []
    for roi in range(nROIs):
        inputs.append((data[roi,:], nuisanceRegressors))
    
    residual_ts = np.zeros((nROIs,np.sum(tMask)))
    betas = np.zeros((nROIs,nuisanceRegressors.shape[1]+1)) # All regressors, + constant regressors
    
    # calculate regrssion for each voxel
    for roi4reg in range(nROIs):
    #for roi4reg in range(10):
        print('Running regression for voxel: ', roi4reg)
        result = _regression2(inputs[roi4reg])
        residual_ts[roi4reg,:] = result[1]
        betas[roi4reg,:] = result[0]
        
    #reshape data table back to fMRI shape and save
    residfMRI = residual_ts.reshape(94,94,66,np.sum(tMask))
    nifti_file = nib.Nifti1Image(residfMRI, affine=np.eye(4))
    outputNii = outputdir + subj +'/' + subj + run + '_'+ model +'.nii'
    nib.save(nifti_file, outputNii)

    # save resid & betas
    residOut = gzip.GzipFile(outputdir + subj +'/' + subj + run + '_resid_'+ model +'.npy.gz', 'w')
    np.save(file=residOut, arr=residual_ts)
    residOut.close()
    betasOut = gzip.GzipFile(outputdir + subj +'/' + subj + run + '_betas_'+ model +'.npy.gz', 'w')
    np.save(file=betasOut, arr=betas)
    betasOut.close()

   
#    pool = mp.Pool(processes=nproc)
#    results = pool.map_async(_regression2,inputs).get()
#    pool.close()
#    pool.join()
#    
#    roi = 0
#    for result in results:
#        betas[roi,:] = result[0]
#        residual_ts[roi,:] = result[1]
#        roi += 1
    
#    outname1 = run + '/nuisanceReg_resid_' + model
#    outname2 = run + '/nuisanceReg_betas_' + model

#    outputfilename = outputdir + subj + '_glmOutput_data.h5'
#    h5f = h5py.File(outputfilename,'a')
#    try:
#        h5f.create_dataset(outname1,data=residual_ts)
#        h5f.create_dataset(outname2,data=betas)
#    except:
#        del h5f[outname1], h5f[outname2]
#        h5f.create_dataset(outname1,data=residual_ts)
#        h5f.create_dataset(outname2,data=betas)
#    h5f.close()

def _createMotionSpikeRegressors(relative_rms, subj, run,  spikeReg=.25):
    """
    relative_rms-  time x 1 array (for HCP data, can be obtained from the txt file 'Movement_RelativeRMS.txt'; otherwise see Van Dijk et al. (2011) Neuroimage for approximate calculation
    run         -   Indicate which run this is
    spikeReg    -   generate spike time regressors for motion spikes, using a default threshold of .25mm FD threshold
    """

    nTRs = relative_rms.shape[0]

    motionSpikes = np.where(relative_rms>spikeReg)[0]
    if len(motionSpikes)>0:
        spikeRegressorsArray = np.zeros((nTRs,len(motionSpikes)))

        for spike in range(len(motionSpikes)):
            spike_time = motionSpikes[spike]
            spikeRegressorsArray[spike_time,spike] = 1.0

        spikeRegressorsArray = np.asarray(spikeRegressorsArray,dtype=bool)

        # Create h5py output
        h5f = h5py.File(nuis_reg_dir + subj + '_nuisanceRegressors.h5','a')
        try:
            h5f.create_dataset(run + '/motionSpikes',data=spikeRegressorsArray)
        except:
            del h5f[run + '/motionSpikes']
            h5f.create_dataset(run + '/motionSpikes',data=spikeRegressorsArray)

        h5f.close()

def _createPhysiologicalNuisanceRegressors(inputname, subj, run, globalmask, wmmask, ventriclesmask, aCompCor=5):
    """
    inputname   -   4D input time series to obtain nuisance regressors
    run      -   fMRI run
    globalmask  -   whole brain mask to extract global time series
    wmmask      -   white matter mask (functional) to extract white matter time series
    ventriclesmask- ventricles mask (functional) to extract ventricle time series
    aCompCor    -   Create PC component time series of white matter and ventricle time series, using first n PCs
    """
    
    # Load raw fMRI data (in volume space)
    print('Loading raw fMRI data')
    fMRI4d = nib.load(inputname).get_fdata()

    ##########################################################
    ## Nuisance time series (Global signal, WM, and Ventricles)
    print('Obtaining standard global, wm, and ventricle signals and their derivatives')
    # Global signal
    globalMask = nib.load(globalmask).get_fdata()
    globalMask = np.asarray(globalMask,dtype=bool)[:,:,:,0]
    globaldata = fMRI4d[globalMask].copy()
    globaldata = signal.detrend(globaldata,axis=1,type='constant')
    globaldata = signal.detrend(globaldata,axis=1,type='linear')
    global_signal1d = np.mean(globaldata,axis=0)
    # White matter signal
    wmMask = nib.load(wmmask).get_fdata()
    wmMask = np.asarray(wmMask,dtype=bool)[:,:,:,0]
    wmdata = fMRI4d[wmMask].copy()
    wmdata = signal.detrend(wmdata,axis=1,type='constant')
    wmdata = signal.detrend(wmdata,axis=1,type='linear')
    wm_signal1d = np.mean(wmdata,axis=0)
    # Ventricle signal
    ventricleMask = nib.load(ventriclesmask).get_fdata()
    ventricleMask = np.asarray(ventricleMask,dtype=bool)[:,:,:,0]
    ventricledata = fMRI4d[ventricleMask].copy()
    ventricledata = signal.detrend(ventricledata,axis=1,type='constant')
    ventricledata = signal.detrend(ventricledata,axis=1,type='linear')
    ventricle_signal1d = np.mean(ventricledata,axis=0)

    del fMRI4d

    ## Create derivative time series (with backward differentiation, consistent with 1d_tool.py -derivative option)
    # Global signal derivative
    global_signal1d_deriv = np.zeros(global_signal1d.shape)
    global_signal1d_deriv[1:] = global_signal1d[1:] - global_signal1d[:-1]
    # White matter signal derivative
    wm_signal1d_deriv = np.zeros(wm_signal1d.shape)
    wm_signal1d_deriv[1:] = wm_signal1d[1:] - wm_signal1d[:-1]
    # Ventricle signal derivative
    ventricle_signal1d_deriv = np.zeros(ventricle_signal1d.shape)
    ventricle_signal1d_deriv[1:] = ventricle_signal1d[1:] - ventricle_signal1d[:-1]

    ## Write to h5py
    # Create h5py output
    h5f = h5py.File(nuis_reg_dir + subj + '_nuisanceRegressors.h5','a')
    try:
        h5f.create_dataset(run + '/global_signal',data=global_signal1d)
        h5f.create_dataset(run + '/global_signal_deriv',data=global_signal1d_deriv)
        h5f.create_dataset(run + '/wm_signal',data=wm_signal1d)
        h5f.create_dataset(run + '/wm_signal_deriv',data=wm_signal1d_deriv)
        h5f.create_dataset(run + '/ventricle_signal',data=ventricle_signal1d)
        h5f.create_dataset(run + '/ventricle_signal_deriv',data=ventricle_signal1d_deriv)
    except:
        del h5f[run + '/global_signal'], h5f[run + '/global_signal_deriv'], h5f[run + '/wm_signal'], h5f[run + '/wm_signal_deriv'], h5f[run + '/ventricle_signal'], h5f[run + '/ventricle_signal_deriv']
        h5f.create_dataset(run + '/global_signal',data=global_signal1d)
        h5f.create_dataset(run + '/global_signal_deriv',data=global_signal1d_deriv)
        h5f.create_dataset(run + '/wm_signal',data=wm_signal1d)
        h5f.create_dataset(run + '/wm_signal_deriv',data=wm_signal1d_deriv)
        h5f.create_dataset(run + '/ventricle_signal',data=ventricle_signal1d)
        h5f.create_dataset(run + '/ventricle_signal_deriv',data=ventricle_signal1d_deriv)

    
    ##########################################################
    ## Obtain aCompCor regressors using first 5 components of WM and Ventricles (No GSR!)
    ncomponents = 5
    nTRs = len(global_signal1d)
    print('Obtaining aCompCor regressors and their derivatives')
    # WM time series
    #wmstart = time.time()
    # Obtain covariance matrix, and obtain first 5 PCs of WM time series
    tmpcov = np.corrcoef(wmdata.T)
    eigenvalues, topPCs = scipy.sparse.linalg.eigs(tmpcov,k=ncomponents,which='LM')
    # Now using the top n PCs 
    aCompCor_WM = topPCs
#    wmend = time.time() - wmstart
#    print 'WM aCompCor took', wmend, 'seconds'
    
    # Ventricle time series
    #ventstart = time.time()
    # Obtain covariance matrix, and obtain first 5 PCs of ventricle time series
    tmpcov = np.corrcoef(ventricledata.T)
    eigenvalues, topPCs = scipy.sparse.linalg.eigs(tmpcov,k=ncomponents,which='LM')
    # Now using the top n PCs
    aCompCor_ventricles = topPCs
#    ventricletime = time.time() - ventstart 
#    print 'Ventricle aCompCor took', ventricletime, 'seconds' 
    
    # White matter signal derivative using backwards differentiation
    aCompCor_WM_deriv = np.zeros(aCompCor_WM.shape)
    aCompCor_WM_deriv[1:,:] = np.real(aCompCor_WM[1:,:]) - np.real(aCompCor_WM[:-1,:])
    # Ventricle signal derivative
    aCompCor_ventricles_deriv = np.zeros(aCompCor_ventricles.shape)
    aCompCor_ventricles_deriv[1:,:] = np.real(aCompCor_ventricles[1:,:]) - np.real(aCompCor_ventricles[:-1,:])

    ## Write to h5py
    try:
        h5f.create_dataset(run + '/aCompCor_WM',data=aCompCor_WM)
        h5f.create_dataset(run + '/aCompCor_WM_deriv',data=aCompCor_WM_deriv)
        h5f.create_dataset(run + '/aCompCor_ventricles',data=aCompCor_ventricles)
        h5f.create_dataset(run + '/aCompCor_ventricles_deriv',data=aCompCor_ventricles_deriv)
    except:
        del h5f[run + '/aCompCor_WM'], h5f[run + '/aCompCor_WM_deriv'], h5f[run + '/aCompCor_ventricles'], h5f[run + '/aCompCor_ventricles_deriv']
        h5f.create_dataset(run + '/aCompCor_WM',data=aCompCor_WM)
        h5f.create_dataset(run + '/aCompCor_WM_deriv',data=aCompCor_WM_deriv)
        h5f.create_dataset(run + '/aCompCor_ventricles',data=aCompCor_ventricles)
        h5f.create_dataset(run + '/aCompCor_ventricles_deriv',data=aCompCor_ventricles_deriv)


    ##########################################################
    ## Load motion parameters, and calculate motion spike regressors

    h5f.close()

def _regression2(tupleInput):
    """
    Hand coded OLS regression using closed form equation: betas = (X'X)^(-1) X'y
    """
    
    Y,regressors = tupleInput
    # Add 'constant' regressor
    regressors = sm.add_constant(regressors)
    X = regressors.copy()
    try:
#        #C_ss_inv = np.linalg.inv(np.dot(X.T,X))
        C_ss_inv = np.linalg.pinv(np.dot(X.T,X))
    except np.linalg.LinAlgError as err:
        C_ss_inv = np.linalg.pinv(np.cov(X.T))
    betas = np.dot(C_ss_inv,np.dot(X.T,Y.T))
    resid = Y - (betas[0] + np.dot(X[:,1:],betas[1:])).T
    return betas, resid

