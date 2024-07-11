function[] = preprocessing_fMRI_ISC()

% batch script for preprocessing of fMRI data using spm12
% modified from spm8 script (mind and brain, charit?, c stelzel) by IT

spm_dir = '...\spm12\'; % complete spm12 folder with backslash at the end
addpath(genpath(spm_dir));
spm('defaults', 'FMRI');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% General Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For not-Windows user:
% check all slashes, path names might be diffrent
% du Replace all slashes use Command Key ⌘ + F (Apple) or Ctrl + F
% (Linux (an also Windows))

study_dir = '...\ISC_Analysis\';
% Directory where subject folders with subfolders should be saved

sliceTimes = load('...\minimally_Processed\task_fMRI\slice_times.mat'); 

data_dir = '...\';

% enter subjects only subjects that had the new beamer
subject_dir = ['HFG_121'];

number_subjects = size(subject_dir, 1);

TR = 1;               % enter Repetition Time
number_slices = 66;   % enter number of slices
fwhm = 6;
HPfilter = 160; %use 160 like in RICHIE I
GM_masking_final_images = true;

%for multiband slice_time and reference will be calculated in the slice
%timing step

runs = 4;

conditions = {'Single_1','Single_2','Single_3','Single_4','Single_5',...
    'Triplet_1','Triplet_2','Triplet_3','Triplet_4','Triplet_5',...
    'Nonet_1','Nonet_2','Nonet_3','Nonet_4','Nonet_5',...
    'Correct_1','Correct_2','Correct_3','Correct_4','Correct_5'};
levels = {'Single','Triplet','Nonet','Correct';...
    'Single','Triplet','Nonet','Complete'};

mat_dir = '\Design\';    % place for design file
analysis_dir = '\Model\'; % place for model directory


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Preprocessing Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use the tic toc function to switch on and switch off function with %-sign
% firstly the one before Reorientation
% the reorient manually
% the switch on all teh rest

tic
% PREPARTION
%  make_folders() %folders is provided
%  move_allNii_files() % files are provided 

% PREPROCESSING % Start here for the example participant
  slicetiming()
  realignment()
  coregistration()
  segmentation() %% also get normalized masks
  normalisation_epis()
  normalisation_anat()
  smoothing()

% REGRESSION
  extract_WM_CSF_signals() % and button presses
  model_specification()
  model_estimation()

% AFTER PROCESSING
  make_GM_mask()
  organize_residuals()
  z_scale_epochs() % also applies an GM mask if you want (GM_masking_final_images = true)
  condtions_4D() % two options: for comparison, needs exakt same number of trials in all blocks across conditions
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create folders

    function [] = make_folders()
        
        for i=1:number_subjects
            if ~exist(fullfile(study_dir,subject_dir(i,:)),'dir')
                fprintf('Creating new subject directory %s.\n',subject_dir(i,:));
                success = mkdir(fullfile(study_dir,subject_dir(i,:)));
                if ~success, error('Cannot create new subject directory.\n'); end
            end
            if ~exist(fullfile(study_dir,subject_dir(i,:),'Anatomy'),'dir')
                mkdir(fullfile(study_dir,subject_dir(i,:)),'Anatomy');
            end
            if ~exist(fullfile(study_dir,subject_dir(i,:),'Epis'),'dir')
                mkdir(fullfile(study_dir,subject_dir(i,:)),'Epis');
            end
            
            if ~exist(fullfile(study_dir,subject_dir(i,:),'Design'),'dir')
                mkdir(fullfile(study_dir,subject_dir(i,:)),'Design');
            end
            if ~exist(fullfile(study_dir,subject_dir(i,:),'Model'),'dir')
                mkdir(fullfile(study_dir,subject_dir(i,:)),'Model');
            end
            for run = 1:runs
                if ~exist(fullfile(study_dir,subject_dir(i,:),'Epis',['Run_',char(string(run))]),'dir')
                    mkdir(fullfile(study_dir,subject_dir(i,:),'Epis'),['Run_',char(string(run))]);
                end
            end
            
            % Anatomy File Move anat file here since move all nifit doesn't
            % work here
            anatFile = fullfile(data_dir, 'minimally_Processed', 'Anatomy', ...
                [subject_dir(i,:), '_D1.nii']);
            destFile = fullfile(study_dir,subject_dir(i,:),'Anatomy','_anat.nii'); 
            copyfile(anatFile, destFile); 
        end
        fprintf('Creating Directories:\n Done.\n')
    end

%% Moving All Nii Files

    function [] = move_allNii_files()
        for i = 1:number_subjects
            
            disp('Moving Files For: ')
            disp(subject_dir(i,:))
            
            for r = 1:runs
                disp(['Run: ', num2str(r)]);
                
                curRun = ['Run_',num2str(r)];
                
                epiRunDir = fullfile(data_dir,subject_dir(i,:),'Epis',curRun);
                [epiRunFiles, ~] = spm_select('List',epiRunDir, '^_.*\.nii');
                files_for_moving = strcat(epiRunDir,'/',epiRunFiles(6:end,:));
                
                for f = 1:size(files_for_moving,1)
                    destDir = fullfile(study_dir,subject_dir(i,:),'Epis',curRun);
                    copyfile(files_for_moving(f,:),destDir);
                end
            end
           
            
            disp(['DONE ', num2str(i),'/', num2str(number_subjects)]);
        end
    end


%% Slice Timing

    function [] = slicetiming()
        
        for i = 1:number_subjects
            
            disp('Slice Timing for: ')
            disp(subject_dir(i,:))
            
            % define slice_times/slice_order and reference slice from dicom
            slice_times = sliceTimes.slice_times;
            [~,slice_order] = sort(slice_times);
            
            ref_slice = ceil(median(slice_order));
            ref_time = slice_times(ref_slice);
            
            files_for_slicetime = cell(1,4);
            for run = 1:runs
                % reading in all Epis
                epi_dir = [study_dir, subject_dir(i,:), '\Epis\Run_',char(string(run)),'\'];
                [func_filenames, ~] = spm_select('List',epi_dir,'nii');
                files_for_slicetime_run = strcat(epi_dir,func_filenames);
                files_for_slicetime{1,run} = cellstr(files_for_slicetime_run);
            end
            
            % Batch commands
            jobs{1}.spm.temporal.st.scans = files_for_slicetime;
            jobs{1}.spm.temporal.st.nslices = number_slices;
            jobs{1}.spm.temporal.st.tr = TR;
            jobs{1}.spm.temporal.st.ta = TR - (TR/number_slices); % enter TA (TR-(TR/slice_number)), sometimes problematic
            % TA is not relevant in slice timing for multiband, since it is
            % already implicitly included in the slice times and thus does
            % not need to be specified
            jobs{1}.spm.temporal.st.so = slice_times;
            jobs{1}.spm.temporal.st.refslice = ref_time;
            jobs{1}.spm.temporal.st.prefix = 'a';
            % RUN
            spm_jobman('run', jobs);
            clear jobs
            disp(['DONE ', num2str(i),'/', num2str(number_subjects)]);
        end
        fprintf('Slice Timing:\n Done.\n')
    end


%% Realignment

    function [] = realignment()
        
        for i = 1:number_subjects
            
            disp('Realignment for: ')
            disp(subject_dir(i,:))
                      
            files_for_realign = cell(1,4);
            for run = 1:runs
                % reading in all Epis
                epi_dir = [study_dir, subject_dir(i,:), '\Epis\Run_',char(string(run)),'\'];
                [func_filenames, ~] = spm_select('List',epi_dir,'^a.*\.nii');
                files_for_realign_run = strcat(epi_dir,func_filenames);
                files_for_realign{1,run} = cellstr(files_for_realign_run);
            end
            
            % Batch commands
            jobs{1}.spm.spatial.realign.estwrite.data = files_for_realign;
            jobs{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;      % 0.9 is default
            jobs{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            jobs{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            jobs{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1; % register to first or mean: 1 means mean, 0 means first
            %initially SPM realigns each session to each other by aligning the first volume
            %from each session to the first volume of the first session and subsequently
            %all volumes within each session are aligned to the first volume of that session
            %after that the volumes from the first realignment step are used to create
            %a mean image and than all volumes are aligned to that mean image
            jobs{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            jobs{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            jobs{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
            jobs{1}.spm.spatial.realign.estwrite.roptions.which = [0 1]; % mean image only
            jobs{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
            jobs{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            jobs{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
            jobs{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
            % RUN
            spm_jobman('run', jobs);
            clear jobs
            
        end
        fprintf('Realignment:\n Done.\n')
    end

%% Coregistration

    function [] = coregistration()
        
        for i = 1:number_subjects
            
            disp('Coregistration for: ')
            disp(subject_dir(i,:))
            
            % read in Mean image from Realignment and Anatomy
            %%%%%      % mean image is only in the Run_1 Folder ???
            meanepi_dir = [study_dir, subject_dir(i,:), '/Epis/Run_1'];
            find_meanepi = dir(fullfile(meanepi_dir,'mean*.nii'));
            meanepi_file = find_meanepi.name;
            
            highres_dir = [study_dir, subject_dir(i,:), '/Anatomy'];
            find_highres = dir(fullfile(highres_dir, '_anat*.nii'));
            highres_file = find_highres.name;
            
            % define Dar as in Batch
            jobs{1}.spm.spatial.coreg.estwrite.ref = {[meanepi_dir filesep meanepi_file]}; % reference image: mean epi (where to coreg to)
            jobs{1}.spm.spatial.coreg.estwrite.source = {[highres_dir filesep highres_file]}; % source image: anatomy (what to coreg)
            %Anatomy diffent session? DTI if possible???
            jobs{1}.spm.spatial.coreg.estwrite.other = {''};
            jobs{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
            jobs{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
            jobs{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            jobs{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
            jobs{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
            jobs{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
            jobs{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
            jobs{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
            % RUN
            spm_jobman('run', jobs);
            clear jobs
            
        end
        fprintf('Coregistration:\n Done.\n')
    end

%% Segmentation

    function [] = segmentation()
        
        for i = 1:number_subjects
            
            disp('Segmentation for: ')
            disp(subject_dir(i,:))
            
            % load anatomy
            highres_dir = [study_dir, subject_dir(i,:), '/Anatomy'];
            find_highres = dir(fullfile(highres_dir, '_anat*.nii'));
            highres_file = find_highres.name;
            
            % Batch commands (lots of it is just the parameters for the
            % different materials
            jobs{1}.spm.spatial.preproc.channel.vols = {[highres_dir filesep highres_file]};
            jobs{1}.spm.spatial.preproc.channel.biasreg = 0.001;
            jobs{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            jobs{1}.spm.spatial.preproc.channel.write = [0 1];
            jobs{1}.spm.spatial.preproc.tissue(1).tpm = {[spm_dir,'tpm',filesep,'TPM.nii,1']}; % grey matter
            jobs{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
            jobs{1}.spm.spatial.preproc.tissue(1).native = [1 1]; % for dartel just in case
            jobs{1}.spm.spatial.preproc.tissue(1).warped = [1 1];
            jobs{1}.spm.spatial.preproc.tissue(2).tpm = {[spm_dir,'tpm',filesep,'TPM.nii,2']}; % white matter
            jobs{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
            jobs{1}.spm.spatial.preproc.tissue(2).native = [1 1]; % for dartel just in case
            jobs{1}.spm.spatial.preproc.tissue(2).warped = [1 1];
            jobs{1}.spm.spatial.preproc.tissue(3).tpm = {[spm_dir,'tpm',filesep,'TPM.nii,3']}; % CSF
            jobs{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
            jobs{1}.spm.spatial.preproc.tissue(3).native = [1 0];
            jobs{1}.spm.spatial.preproc.tissue(3).warped = [1 1]; % save the noramlized
            jobs{1}.spm.spatial.preproc.tissue(4).tpm = {[spm_dir,'tpm',filesep,'TPM.nii,4']}; % bone
            jobs{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
            jobs{1}.spm.spatial.preproc.tissue(4).native = [1 0];
            jobs{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
            jobs{1}.spm.spatial.preproc.tissue(5).tpm = {[spm_dir,'tpm',filesep,'TPM.nii,5']}; % soft tissue
            jobs{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
            jobs{1}.spm.spatial.preproc.tissue(5).native = [1 0];
            jobs{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
            jobs{1}.spm.spatial.preproc.tissue(6).tpm = {[spm_dir,'tpm',filesep,'TPM.nii,6']}; % air/background
            jobs{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
            jobs{1}.spm.spatial.preproc.tissue(6).native = [0 0];
            jobs{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
            jobs{1}.spm.spatial.preproc.warp.mrf = 1;
            jobs{1}.spm.spatial.preproc.warp.cleanup = 1;
            jobs{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
            jobs{1}.spm.spatial.preproc.warp.affreg = 'mni'; % Standard brain
            jobs{1}.spm.spatial.preproc.warp.fwhm = 0;
            jobs{1}.spm.spatial.preproc.warp.samp = 3;
            jobs{1}.spm.spatial.preproc.warp.write = [1 1]; % foward
            % RUN
            spm_jobman('run', jobs);
            clear jobs
        end
        fprintf('Segmentation:\n Done.\n')
    end

%% Normalisation Epis

    function [] = normalisation_epis()
        
        for i = 1:number_subjects
            
            disp('Writing the normalised Epis for: ')
            disp(subject_dir(i,:))
            
            highres_dir = [study_dir, subject_dir(i,:), '\Anatomy'];
            find_deform = dir(fullfile(highres_dir, 'y_*'));
            deform_file = find_deform.name;
            
            for run = 1:runs
                epi_dir = [study_dir, subject_dir(i,:), '\Epis\Run_',char(string(run)),'\'];
                [func_filenames_o, ~] = spm_select('List',epi_dir,'^a.*\.nii');
                files_for_norm = strcat(epi_dir,func_filenames_o);
                jobs{1}.spm.spatial.normalise.write.subj.def = {[highres_dir filesep deform_file]};
                jobs{1}.spm.spatial.normalise.write.subj.resample = cellstr(files_for_norm);
                jobs{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                    78 76 85];
                % Voxel size/dimensions
                jobs{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2]; % enter voxel size of normalised epis
                jobs{1}.spm.spatial.normalise.write.woptions.interp = 4;
                jobs{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
                % RUN
                spm_jobman('run', jobs);
                clear jobs
                disp(['DONE: RUN ', num2str(run)]);
            end
        end
        fprintf('Normalisation of Epis:\n Done.\n')
    end

%% Normalisation Anatomy

    function [] = normalisation_anat()
        
        for i = 1:number_subjects
            
            disp('Writing the normalised Anatomy for: ')
            disp(subject_dir(i,:))
            
            % Same procedure for Anatomy
            highres_dir = [study_dir, subject_dir(i,:), '\Anatomy'];
            find_deform = dir(fullfile(highres_dir, 'y_*'));
            deform_file = find_deform.name;
            find_highres = dir(fullfile(highres_dir, 'm_anat*.nii'));
            highres_file = find_highres.name;
            
            jobs{1}.spm.spatial.normalise.write.subj.def = {[highres_dir filesep deform_file]};
            jobs{1}.spm.spatial.normalise.write.subj.resample = {[highres_dir filesep highres_file]};
            jobs{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                78 76 85];
            % now smaller Voxel size
            jobs{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1]; % enter voxel size of normalised anatomy
            jobs{1}.spm.spatial.normalise.write.woptions.interp = 4;
            jobs{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
            
            spm_jobman('run', jobs);
            clear jobs
        end
        fprintf('Normalisation of Anatomy:\n Done.\n')
    end


%% Smoothing

    function [] = smoothing()
        
        for i = 1:number_subjects
            
            disp('Smoothing for: ')
            disp(subject_dir(i,:))
            
            %load Epis
            for run = 1:runs
                epi_dir = [study_dir, subject_dir(i,:), '\Epis\Run_',char(string(run)),'\'];
                [func_filenames_sm, ~] = spm_select('List',epi_dir,'^wa.*\.nii');
                files_for_smooth = strcat(epi_dir,func_filenames_sm);
                
                % Batch commands
                jobs{1}.spm.spatial.smooth.data = cellstr(files_for_smooth);
                jobs{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];  % enter FWHM
                jobs{1}.spm.spatial.smooth.dtype = 0;
                jobs{1}.spm.spatial.smooth.im = 0;
                jobs{1}.spm.spatial.smooth.prefix = ['s',num2str(fwhm)]; % Renaming here makes Renaming smooth unecessary
                %RUN
                spm_jobman('run', jobs);
                clear jobs
                disp(['DONE: RUN ', num2str(run)]);
            end
        end
        fprintf('Smoothing:\n Done.\n')
    end

%% Grey Matter Mask
    function[] = make_GM_mask()
        normGM = cell(2,1); %cell(number_subjects,1);
        expression_string = '(';
        for i = 1:number_subjects
            normGM{i,1} = fullfile(study_dir, subject_dir(i,:), 'Anatomy', 'wc1_anat.nii');
            expression_string = [expression_string,'i', num2str(i), '+'];
        end
        expression_string = [expression_string(1:(end-1)), ')/',num2str(number_subjects)];
        clear jobs;
        
        jobs{1}.spm.util.imcalc.input = normGM;
        jobs{1}.spm.util.imcalc.output = 'norm_GM_mask.nii';
        jobs{1}.spm.util.imcalc.outdir = {study_dir};
        jobs{1}.spm.util.imcalc.expression = expression_string;
        jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        jobs{1}.spm.util.imcalc.options.dmtx = 0;
        jobs{1}.spm.util.imcalc.options.mask = 0;
        jobs{1}.spm.util.imcalc.options.interp = 1;
        jobs{1}.spm.util.imcalc.options.dtype = 4;
        
        spm_jobman('run', jobs);
        clear jobs;
        
        % resize to later fit image dimensions
        epi_dir_run1 = fullfile(study_dir, subject_dir(1,:), 'Epis', 'Run_1');
        [func_filenames_wa, ~] = spm_select('List',epi_dir_run1,'^wa.*\.nii');
        jobs{1}.spm.spatial.coreg.write.ref = {fullfile(epi_dir_run1,func_filenames_wa(15,:))}; % 15 ist just any number
        jobs{1}.spm.spatial.coreg.write.source = {fullfile(study_dir,'norm_GM_mask.nii')};
        jobs{1}.spm.spatial.coreg.write.roptions.interp = 4;
        jobs{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
        jobs{1}.spm.spatial.coreg.write.roptions.mask = 0;
        jobs{1}.spm.spatial.coreg.write.roptions.prefix = 'resized_';
        
        spm_jobman('run', jobs);
        clear jobs;
        
        % Smoothing in same dimensions as functional kernel
        jobs{1}.spm.spatial.smooth.data = {fullfile(study_dir, 'resized_norm_GM_mask.nii')};
        jobs{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
        jobs{1}.spm.spatial.smooth.dtype = 0;
        jobs{1}.spm.spatial.smooth.im = 0;
        jobs{1}.spm.spatial.smooth.prefix = ['s',num2str(fwhm),'_'];
        
        spm_jobman('run', jobs);
        clear jobs;
        
        % threshold at 0.25
        jobs{1}.spm.util.imcalc.input = {[study_dir,'/s',num2str(fwhm),'_resized_norm_GM_mask.nii']};
        jobs{1}.spm.util.imcalc.output = 'bin_GM_mask.nii';
        jobs{1}.spm.util.imcalc.outdir = {study_dir};
        jobs{1}.spm.util.imcalc.expression = 'i1 > 0.25'; % based on Nastase et al., 2019
        jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        jobs{1}.spm.util.imcalc.options.dmtx = 0;
        jobs{1}.spm.util.imcalc.options.mask = 0;
        jobs{1}.spm.util.imcalc.options.interp = 1;
        jobs{1}.spm.util.imcalc.options.dtype = 4;
        
        spm_jobman('run', jobs);
        clear jobs;
        
        fprintf('Grey Matter Mask:\n Done.\n')
    end
%% Extrcat the WM and CSF signal for regression && DesignFile
    function [] = extract_WM_CSF_signals()
        % normalized WM and CSF masks to roi mat
        for i = 1:number_subjects
            disp('Extracting avg. WM and CF signal for: ')
            disp(subject_dir(i,:))
            
            subAnatDir = fullfile(study_dir, subject_dir(i,:), 'Anatomy');
            for m = 2:3
                imgname = [subAnatDir, filesep, 'wc', num2str(m),'_anat.nii'];
                MarsO = maroi_image(struct('vol', spm_vol(imgname), 'binarize',0,'func', 'img'));
                saveroi(MarsO, [subAnatDir, filesep, 'wc', num2str(m),'_anat_roi.mat']);
            end
            % I just want to get raw timecourses from some images; how do I do that?
            %
            % This can be done from the GUI. Select Data - Extract ROI data (full options).
            % Select the ROIs, say No to use SPM design, Other for type of images,
            % 1 for number of subjects. Select the images you want to extract data from,
            % Raw data for scaling, and 0 for grand mean. Now you can plot the data from the GUI,
            % or save in various formats using Data - Export data. The script to do this might be:
            
            
            subFile = [study_dir, subject_dir(i,:), filesep, subject_dir(i,:), '_log.txt'];
            
            Data = readtable(subFile);
            Labels = {'ID','block','trial_block','trial_cat','trial_total','number','violation',...
                'level_violation','onset','duration','Response','Hit_Miss', 'RT','pre_trial_fixcross_ISI',...
                'onset_correct'};
            Data.Properties.VariableNames = Labels(1:end-1);
            Data = rmmissing(Data);
            
            
            % get onsets/triggers
            onsetData = readmatrix(subFile, ...
                delimitedTextImportOptions('DataLines',[2,Inf]),'OutputType', 'string'); % read logfile in a form, so onsets are alse read
            onsetData = onsetData([1,2,813,814,1625,1626,2437,2438,3249],:); % get lines with onsets (for everybody the same, sinnce equal block length)
            onset_times = zeros(4,2);
            onset1 = split(onsetData(1,:)," ",2);
            onsetRest = split(onsetData(2:end,:)," ",2); % extract  trigger time information from logfile
            onset_times(1,:) = [1, str2double(onset1(1,2))];
            onsetRest_short = str2double(onsetRest(2:(end-1),[3,5])); % just get important columns
            onsetRest_short(:,3) = onsetRest_short(:,2) - [0; onsetRest_short(1:(end-1),2)]; % get distances between triggers
            onset_times(2:4,:) = onsetRest_short([2,4,6],[1,3]); % the triggers after the break give the addidtional offset between presentation and the scanner
            onset_times(:,3) = cumsum(onset_times(:,2)); % accumulated offset, necessary if time would run continously for the uns
            onset_times(:,4) = [onset1(1,2);onsetRest([3,5,7],5)]; % Starting times after break, necessary as every run starts modelling at t0
            % since runs/sessions are handled as new subjects in first level
            % modelling --> fixed effects model with multiple elements
            % we need to the modelling for each run effectively at 0
            
            Data.onset_correct = zeros(size(Data,1),1);
            
            % onsets in sekunden, ersten trigger (nach pause) subtrahieren
            for part = 1:4
                Data.onset_correct((810*(part-1)+1):(part*810)) = (Data.onset((810*(part-1)+1):(part*810))- onset_times(part,4))/1000;
            end
            Data.onset_correct_adjusted = Data.onset_correct - (5 * TR);
            Data.onset_scans = ceil(Data.onset_correct_adjusted);
            
            %----------------------------------------------------------------------------------------------------------
            for run = 1:runs % create the onset-mat-files for each run seperately
                
                runDir = fullfile(study_dir, subject_dir(i,:), 'Epis', ['Run_',num2str(run)]);
                subData = Data((810*(run-1)+1):(run*810),:); % only the data from one of the runs, as the Model-Files need to be run-specific
                
                % button presses: Onset scan button press, duration 0
                button_press = subData(subData.Response == 1,:); % get button presses
                scans_button_press = (button_press.onset_scans);
                durations_button_press = zeros(size(scans_button_press));
                
                save(fullfile(study_dir, subject_dir(i,:),mat_dir, ...
                    [subject_dir(i,:),'_run',sprintf('%02d',run),'_BP.mat']),...
                    'scans_button_press', 'durations_button_press');
                
                [roi_files, ~] = spm_select('List',subAnatDir, '.*roi.mat');
                roi_paths = cellstr(strcat(subAnatDir,'\',roi_files));
                [epi_files, ~] = spm_select('List',runDir,['^s',num2str(fwhm),'wa.*\.nii']);
                P = strcat(runDir,'\',epi_files);
                rois = maroi('load_cell', roi_paths);  % make maroi ROI objects
                mY = get_marsy(rois{:}, P, 'mean');  % extract data into marsy data object
                y = summary_data(mY); % get summary time course(s)
                y_WM = y(:,1); y_CSF = y(:,2);
                save(fullfile(study_dir, subject_dir(i,:),mat_dir, ...
                    [subject_dir(i,:),'_run',sprintf('%02d',run),'_WM_CSF.mat']),...
                    'y_WM', 'y_CSF');
                disp(['DONE: RUN ', num2str(run)]);
            end
        end
        fprintf('Extraction:\n Done.\n')
    end

%% MODEL SPECIFICATION
    function [] = model_specification()
        
        for i = 1:number_subjects
            tic
            disp('Model specification for:')
            disp(subject_dir(i,:))
            
            clear jobs;
            model_dir = [study_dir, subject_dir(i,:), analysis_dir];
            
            % define slice_times/slice_order and reference slice from dicom
            slice_times = sliceTimes.slice_times;
            [~,slice_order] = sort(slice_times);
            ref_slice = ceil(median(slice_order));
            
            
            jobs{1}.spm.stats{1}.fmri_spec.timing.units = 'scans';
            jobs{1}.spm.stats{1}.fmri_spec.timing.RT = TR;
            jobs{1}.spm.stats{1}.fmri_spec.timing.fmri_t = number_slices;
            jobs{1}.spm.stats{1}.fmri_spec.timing.fmri_t0 = find(slice_order == ref_slice);
            jobs{1}.spm.stats{1}.fmri_spec.bases.hrf.derivs = [0 0];            % Canonical HRF, ohne time und dispersion derivative: sonst [1 0] oder [1 1]
            jobs{1}.spm.stats{1}.fmri_spec.cvi = 'FAST'; % for TR < 1 ist FAST better
            jobs{1}.spm.stats{1}.fmri_spec.dir = {model_dir};
            jobs{1}.spm.stats{1}.fmri_spec.global = 'None';
            jobs{1}.spm.stats{1}.fmri_spec.mthresh = 0.8;
            jobs{1}.spm.stats{1}.fmri_spec.mask = {''};
            %jobs{1}.spm.stats{1}.fmri_spec.mask = {%%ADD GM MASK IF NECESSARY};
            jobs{1}.spm.stats{1}.fmri_spec.volt = 1;
            
            
            matfile_dir = fullfile(study_dir, subject_dir(i,:), mat_dir);
            
            for run = 1:runs
                epi_dir = fullfile(study_dir, subject_dir(i,:),'Epis' ,['Run_',num2str(run)]);
                [func_filenames, ~] = spm_select('List', epi_dir, '^s6wa.*\.nii');  % prefix der vorverarbeiteten EPIS eingeben
                model_files = strcat(epi_dir,'\',func_filenames);
                [movement_reg_file , ~]= spm_select('List', epi_dir, '^rp.*\.txt');
                movement_reg = strcat(epi_dir,'\',movement_reg_file);
                
                bp = load(fullfile(matfile_dir,[subject_dir(i,:),'_run',sprintf('%02d',run),'_BP.mat']));
                mask_signal = load(fullfile(matfile_dir,[subject_dir(i,:),'_run',sprintf('%02d',run),'_WM_CSF.mat']));
                
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).scans = cellstr(model_files);
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).cond.name = 'button_press';
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).cond.onset = bp.scans_button_press;
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).cond.duration = bp.durations_button_press;
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).cond.tmod = 0;
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).cond.orth = 0;
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).multi = {''};
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).regress(1).name = 'WM';
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).regress(1).val = mask_signal.y_WM;
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).regress(2).name = 'CSF';
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).regress(2).val = mask_signal.y_CSF;
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).multi_reg = {movement_reg};
                jobs{1}.spm.stats{1}.fmri_spec.sess(run).hpf = HPfilter;
            end
            
            spm_jobman('run',jobs)
            
            copyfile([model_dir,'SPM.mat'],...
                [model_dir,'SPM.model'])
            
            clear jobs
        end
        
        disp('Modellspezifikation: DONE ')
        toc
    end


%% Modellschaetzung
%==========================================================================
    function [] = model_estimation()
        for i = 1:number_subjects
            tic
            disp('Modellschaetzung fuer:')
            disp(subject_dir(i,:))
            
            model_dir = [study_dir, subject_dir(i,:), analysis_dir];
            jobs{1}.spm.stats{1}.fmri_est.method.Classical = 1;
            jobs{1}.spm.stats{1}.fmri_est.write_residuals = 1;
            jobs{1}.spm.stats{1}.fmri_est.spmmat = {[model_dir 'SPM.mat']};
            
            spm_jobman('run',jobs)
            copyfile([model_dir 'SPM.mat'],[model_dir 'SPM.estmodel'])
            clear jobs
            
            disp('Modellschaetzung: DONE')
            toc
        end
    end

%% Residuen Umbenennen

    function[] = organize_residuals()
        for i = 1:number_subjects
            disp('Organization for:')
            disp(subject_dir(i,:))
            model_dir = fullfile(study_dir, subject_dir(i,:),analysis_dir);
            epis_dir = fullfile(study_dir, subject_dir(i,:),'Epis');
            for run = 1:runs
                if ~exist([model_dir,'Run_',num2str(run)],'dir')
                    mkdir([model_dir,'Run_',num2str(run)]);
                end
                
                runDir = [epis_dir, '\Run_',num2str(run),'\'];
                [files_run_model,~] = spm_select('List',runDir,'^s6wa.*\.nii');
                num_runFiles= size(files_run_model,1);
                
                [files_res, ~] = spm_select('List',model_dir,'^Res_.*\.nii');
                
                for iFile = 1:num_runFiles % looping throgh all of them
                    newName = sprintf(['Res_run',char(string(run)),'_RICHIE_%05d.nii'],(iFile+5));  % %05d means numebr has 5 digits and will be completed by leading 0 if it has no 5 digits on its own                          % enter favoured names of epi files
                    if ~exist(fullfile(model_dir,['Run_',num2str(run)],newName),'file')
                        movefile(fullfile(model_dir,files_res(iFile,:)),fullfile(model_dir,['Run_',num2str(run)],newName));
                        %copyfile(fullfile(model_dir,files_res(iFile,:)),fullfile(model_dir,['Run_',num2str(run)],newName));
                        % move the old file to the ne named file, thus creating
                        % the new files
                    end
                end
            end
            
        end
        disp('Organization: DONE')
    end

%% z_scale_epochs() & also GM mask images
% When citing these tools in publications, please reference https://robjellis.net/tools.html​
    function [] = z_scale_epochs()
        if GM_masking_final_images
            maskPath = fullfile(study_dir, 'bin_GM_mask.nii');
            mask_1 = spm_vol(maskPath);
            mask  = spm_read_vols(mask_1);
        end
        for i = 1:number_subjects
            disp('z-scaling Images for:')
            disp(subject_dir(i,:));
            for run = 1:runs
                disp(['RUN: ', num2str(run)]);
                runDir = fullfile(study_dir, subject_dir(i,:), analysis_dir, ['Run_',num2str(run)]);
                [files_res, ~] = spm_select('List', runDir, '^Res_.*\.nii');
                files1 = strcat(runDir, filesep,files_res);
                newNames_files1 = strcat(runDir,filesep, 'z_', files_res);
                nImages = size(files_res,1);
                f = waitbar(0, 'Starting');
                for t = 1:nImages
                    file1 = files1(t,:);
                    v1n = spm_vol(file1);
                    v1  = spm_read_vols(v1n);  % the actual volume
                    v1(isnan(v1)) = 0;   % replace NaN if they are present
                    
                    v = v1(v1~=0); % all non-zero values (so we don't have to use nanmean or nanstd); will be a vector
                    m  = mean(v);
                    sd = std(v);
                    
                    vol = (v1 - m) / sd;
                    
                    %you could grey matter mask the images here
                    if GM_masking_final_images
                        volFinal = vol .* mask;
                    else
                        volFinal = vol;
                    end
                    
                    % check for appropriate dtype
                    dtype = dtype_check(volFinal(:)); % dtype check from the imcalc toolbox
                    v1n.('dt') = dtype;
                    v1n.fname = newNames_files1(t,:);
                    %save data
                    v1n = spm_write_vol(v1n,volFinal);
                    waitbar(t/nImages, f, sprintf('Progress: %d %%', floor(t/nImages*100)));
                    pause(0.1);
                end
                close(f);
                disp(['DONE: RUN ', num2str(run)]);
            end
        end
        disp('z-scaling Images: DONE')
    end

%% Concatenate to a single 4D image per condition
    function [] = condtions_4D()
        % thereby leave out 6-10 seconds of data per chunk
     %   allOverview = readtable(fullfile(study_dir, 'allSubj_conditionInfo_adjusted.xlsx'));
        allOverview = readtable(fullfile(study_dir, 'allSubj_conditionInfo_combined.xlsx'));
        for i = 1:number_subjects
            disp('4D transformation for:')
            disp(subject_dir(i,:));
            condFiles = cell(numel(conditions),1);
            condMat = load(fullfile(study_dir, subject_dir(i,:), [subject_dir(i,:),'_condImgInfo.mat']));
            subOverview = allOverview(contains(allOverview.Subject,subject_dir(i,:)),:);
            
            outDir = fullfile(study_dir, subject_dir(i,:), 'ISC_files');
            if ~exist(outDir, 'dir')
                mkdir(outDir);
            end
            
            for c = 1:numel(conditions)
                disp(conditions{c});
                curOnset = subOverview.Onset_Img_corrected(contains(subOverview.Condition, conditions{c}));
                curOffset = subOverview.Offset_Img_corrected(contains(subOverview.Condition, conditions{c}));
                curRun = condMat.condOverview{contains(condMat.condOverview(:,1), conditions{c}),2};
                
                curOnset_corrected = curOnset +10;
                
                number_files3D = curOffset - curOnset_corrected +1;
                
                files3D = cell(number_files3D,1);
                
                run = str2double(curRun(end));
                modelRunDir = fullfile(study_dir,subject_dir(i,:),analysis_dir,curRun);
                
                m = 1;
                for iFile = curOnset_corrected:curOffset
                    fileName = sprintf(['z_Res_run',char(string(run)),'_RICHIE_%05d.nii'],(iFile));
                    pathName = fullfile(modelRunDir, fileName);
                    if exist(pathName, 'file')
                        files3D{m,1} =  pathName;
                    else
                        disp(['Does not exist: ', fileName]);
                    end
                    m = m+1;
                end
                condFiles{c,1} = files3D;
            end
            for lev = 1:size(levels,2)
                clear jobs;
                conLev = contains(conditions,levels{1,lev})';
                files_for_4D = cat(1, condFiles{conLev,1});
                jobs{1}.spm.util.cat.vols = files_for_4D;
                jobs{1}.spm.util.cat.name = fullfile(outDir, [levels{2,lev},'_z_Res_RICHIE.nii']);
                jobs{1}.spm.util.cat.dtype = 16; %FLOAT32 (single precision float value)
                jobs{1}.spm.util.cat.RT = TR;
                
                spm_jobman('RUN', jobs);
                clear jobs;
                disp(['DONE concatenating: ', levels{2,lev}]);
            end
        end
        disp('4D transformation: DONE');
    end
end