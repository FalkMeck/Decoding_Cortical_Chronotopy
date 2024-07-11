%% SCRIPT HAS TO BE RUN FROM A UNIX-PC %%
function [] = rest_prep_cato()
% DTI Data preprocessing steps with FSL and CATO toolbox
cd('...'); 

addpath(genpath('.../freesurfer/matlab')); % is part of freesurfer installation
addpath(genpath('.../CATO_3.1.2/CATO-main/src')); % is part of CATO installation

analysisDir = '...';

subject_dir = ['HFG_121']; 

number_subjects = size(subject_dir,1);

outdir = '.../';

%% Run function
tic
   move_files()
   CATO_fun_pipe()
toc
  

%% Move/copy files to prepare for 
    function [] = move_files()
        
        rs_dataDir = '.../Analysis/restingState/';
        
        
        for i = 1:number_subjects
            disp(['Moving files for:', subject_dir(i,:)]); 
            
            % Subject Dir
            if ~exist(fullfile(outdir, subject_dir(i,:)), 'dir')
                mkdir(fullfile(outdir, subject_dir(i,:)));
            endx
            
            % move epis
            if ~exist(fullfile(outdir, subject_dir(i,:), 'fMRI'), 'dir')
                mkdir(fullfile(outdir, subject_dir(i,:), 'fMRI')); 
            end
            epiSrc = [analysisDir, '/minimally_Processed/resting_state_fMRI/', subject_dir(i,:),'_rsfmri.nii.gz'];
            epiDes = fullfile(outdir, subject_dir(i,:), 'fMRI');
            copyfile(epiSrc, epiDes); 
            
            % move Freesurfer
            fsSrc = fullfile(analysisDir, 'minimally_Processed', 'FreeSurfer');
            fsDes = fullfile(outdir, subject_dir(i,:),'FreeSurfer'); 
            copyfile(fsSrc, fsDes); 
            
            % copy the sliceTiming file to each subject folder, might be
            % better for the sliceTimer Config
            sliceTiming = [analysisDir, '/minimally_Processed/resting_state_fMRI/sliceTimesMB.txt']; 
            subDir = fullfile(outdir, subject_dir(i,:)); 
            copyfile(sliceTiming, subDir);
            
            disp('Finished moving');  
       end
    end



%% CATO
    function [] = CATO_fun_pipe()


%         % Functional_pipleline command was adjusted to FREESURFER v7

          configFile = [outdir, '_Code/Functional_configuration_prep.conf']; 
          % only necessary if structurall_pipeline runs via MATLAB and not
          % via Executables/Shell-script
                   
           for i = 1:number_subjects
             
               % create Subject Directory
               subjectDirectory = fullfile(outdir,subject_dir(i,:));
              
              % RUN CATO pipeline via MATLAB (for each participant) 
              functional_pipeline(subjectDirectory,'configurationFile', configFile, 'runType', 'overwrite');

           end
          
    end        
end
