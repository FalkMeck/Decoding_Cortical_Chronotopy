% This MATLAB script takes all the files in inDir that satisfy a certain
% searchstring then it calculates the average time course and all the
% leave-one-out correlations and saves them in outDir. The script uses the
% NIfTI toolbox, so make sure to install it in your path via:
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
addpath('...\NIfTI_20140122');

% Make sure all the 4D files are in the same image space and have the same
% number of time-points or you'll get an error.

% Script written by Christian Keysers
% Adjusted for RICHIE2 (FM)

subject_dir = {'HFG_121'};


study_dir ='...\ISC_Analysis\'; %set this to where all your 4D files are
outDir = fullfile(study_dir, '_Sums4LOO');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
conditions = {'Single', 'Triplet', 'Nonet', 'Complete'}; % create all paths

n = numel(subject_dir);

%% Calculates the overall sum.
% Because most PCs have limited memory, we do not load all the 4D files,
% but build up the mean by simply opening a file, adding it to a running sum,
% and then opening the next file into the same variable. Because we later use
% correlations, where the mean or the sum has the same effect, we use the sum.
for c = 1:numel(conditions)
    % This will setup sum to contain the first image only
    disp(['Calculating ISC maps (LOO) for condition: ', conditions{c}]);
    % ISC files are provided in .../minmimally_Processed/task_fMRI/ISC_files/
    firstFile = [study_dir, subject_dir{1}, filesep, 'ISC_files', filesep,conditions{c}, '_z_Res_RICHIE.nii'];
    sum=load_nii(firstFile);
    
    for i=2:n
        fileName = [study_dir, subject_dir{i}, filesep, 'ISC_files', filesep,conditions{c}, '_z_Res_RICHIE.nii'];
        tmp=load_nii(fileName);
        sum.img=sum.img+tmp.img;
    end
    %% Saves the sum for future use
    save_nii(sum,[outDir,filesep,'sum_',conditions{c},'.nii']);
    sd=std(sum.img,[],4);
    
    %% Gets the sizes of the 4D files to step through all the voxels
    [x,y,z,t]=size(sum.img);
    %% Calculate for each image, the correlation with others and the
    % Fisher-z-transformed maps and saves into outDir.
    % This script uses a simple trick: the average of other participants is not
    % recalculated using all other participants, but simply as the sum minus this
    % particular participant. This allows the process to only need the sum and
    % one particular participant's 4D file to be open simultaneously.
    % To make sure that the correlation images have the same header as the input
    % images, I load the first volume of an image into rho, then zero the values in
    % the .img.
    rho=load_nii(firstFile,1);
    f = waitbar(0,'Please wait...');
    for subj=1:n
        waitbar(subj/n,f,strcat(['starting to process subject ',subject_dir{subj}]));
        fileName = [study_dir, subject_dir{subj}, filesep, 'ISC_files', filesep,conditions{c}, '_z_Res_RICHIE.nii'];
        tmp=load_nii(fileName);
        rho.img=NaN(x,y,z); % maybe need to change to zeros, if ther is a mistake
        
        for xi=1:x
            for yi=1:y
                for zi=1:z
                    if not(sd(xi,yi,zi)==0)  % excludes voxels with zero variance
                        rho.img(xi,yi,zi)=corr(squeeze(tmp.img(xi,yi,zi,:)),squeeze(sum.img(xi,yi,zi,:))-squeeze(tmp.img(xi,yi,zi,:)));
                    end
                end
            end
        end
        
        subOut = fullfile(study_dir, subject_dir{subj}, 'ISC_LOO_maps');
        if ~exist(subOut, 'dir')
            mkdir(subOut);
        end
        save_nii(rho,strcat([subOut,filesep,'ISC_LOO_corr_',conditions{c},'.nii']));
        fisherz=rho;
        fisherz.img=atanh(rho.img);
        save_nii(fisherz,strcat([subOut,filesep,'ISC_LOO_FishZ_',conditions{c},'.nii']));
        disp(['DONE: ', subject_dir{subj}]); 
    end 
end