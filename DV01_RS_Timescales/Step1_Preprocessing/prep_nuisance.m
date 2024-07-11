%% Preparation for Nuisance regression (À la Ito et al. (2020))
function [] = prep_nuisance()
% Needed:
% 1. Maskdir: .../subject_dir(i,:)/masks/
% 1.1 globalmask = Maskdir, subject_dir(i,:), '_wholebrainmask_func_dil1vox.nii.gz'
ROIs{1,1} = [2,4,5,7,8,10,11,12,13,14,15,16,17,18,24,26,28,30,31,41,43,44,46,47,49,50,51,52,53,54,58,60,63,77,85,251,252,253,254,255,1001,1002,1003,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,2001,2002,2003,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030,2031,2032,2033,2034,2035];
% 1.2 wmmask = Maskdir, subject_dir(i,:), '_wmMask_func_eroded.nii.gz'
ROIs{2,1} = [2,41,251,252,253,254,255];
% 1.3 ventriclemask = Maskdir, subject_dir(i,:), '_ventricles_func_eroded.nii.gz'
ROIs{3,1} = [4,5,14,15,24,43,44,72,75,76,213,221];
% 2. Movement regrssors /Movement_Regressors.txt
% - Textfile of 6 rigid movement regressors
% 2.1 Columns 7:12 derivitives of movement regressors
% 2.2 relativeRMS: sqrt(x² + y² + z²), where x, y, and z are motion displacement parameters
%      e.g., x = x[t] - x[t-1]; y = y[t] - y[t-1]; z = z[t] - z[t-1]

% 3. make Nii.gz a nii, load 4D nifti and make it a 2D csv file
%  -  TODO: Test if loading csv and can be remaped to 4D nifti

%% Prep

analysisDir = '.../';

subject_dir = ['HFG_121'];
number_subjects = size(subject_dir,1);

tic
  make_masks()
  movement_regressors()
toc


%% Step 1 creating the masks
function [] = make_masks()
addpath('.../NIfTI_20140122'); % Nifti tools



for i = 1:number_subjects
    % Creating maskDir
    maskDir = fullfile(analysisDir, subject_dir(i,:), 'masks');
    if ~exist(maskDir, 'dir')
        mkdir(maskDir);
    end
       
    
    % Extract parcellation as nii
    gunzip([analysisDir, subject_dir(i,:), 'minimally_Processed/resting_state_fMRI/fMRI_processed/', subject_dir(i,:),'_aparc+aseg_ref.nii.gz'], maskDir);
    
    % load parcellation file
    aparcIm = load_untouch_nii([maskDir,'/', subject_dir(i,:), '_aparc+aseg_ref.nii']);
    
    % whole brain
    outputIm.img = uint32(ismember(aparcIm.img, ROIs{1,1}));
    outputImFn = [maskDir,'/', subject_dir(i,:),'_wholebrainmask.nii'];
    outputIm.hdr = aparcIm.hdr;
    save_nii(outputIm,outputImFn);
    gzip(outputImFn);
    
    % White Matter
    outputIm.img = uint32(ismember(aparcIm.img, ROIs{2,1}));
    outputImFn = [maskDir,'/', subject_dir(i,:),'_wmMask.nii'];
    outputIm.hdr = aparcIm.hdr;
    save_nii(outputIm,outputImFn);
    gzip(outputImFn);
    
    % Ventricle
    outputIm.img = uint32(ismember(aparcIm.img, ROIs{3,1}));
    outputImFn = [maskDir,'/', subject_dir(i,:),'_ventricles.nii'];
    outputIm.hdr = aparcIm.hdr;
    save_nii(outputIm,outputImFn);
    gzip(outputImFn);
    
end
end

%% Step 2 creating movement regressors text files
function [] = movement_regressors()

for i = 1:number_subjects
    
    % Movement Regressor with ridig movement and derivatives
    motionPar = dlmread([analysisDir,subject_dir(i,:),'/fMRI_processed/', subject_dir(i,:), '_motion_parameters.par']); 
    motionMet = load([analysisDir,subject_dir(i,:),'/fMRI_processed/', subject_dir(i,:), '_motion_metrics.mat']); 

    % TO DO check derivitive of movement, but I think that is how CATO does
    % it
    deriv_motionPar = zeros(size(motionPar)); 
    for m = 1:size(motionPar,2)
        deriv_motionPar(2:end,m) = diff(motionPar(:,m));
    end
    
    moveReg = [motionPar, deriv_motionPar];
    dlmwrite([analysisDir,subject_dir(i,:),'/Movement_Regressors.txt'],moveReg, 'delimiter', '\t');
    
    % relativeRMS file, is used in the same way as Framewise Displacement
    % is used in CATO so we will use FD (FD is first row in MotionMetrics)
    dlmwrite([analysisDir,subject_dir(i,:),'/Movement_RelativeRMS.txt'],motionMet.motionMetrics(1,:)');
end
end

end
