%% Stepwise fucntional connectivity

% fMRI preprocessing
% 1. removal of first volumes % DONE
% 2. slice time correction % DONE FSL
% 3. motion correction % DONE FSL
% 8. nuisance regression: 6 ridig motion, whole brain signal, lateral
% ventricle signal, white matter signal % DONE more throughly (Ito et al., 2020)
% 4. normalization
% 5. resampled to 2mm cubic voxels
% 6. preserving MNI coordinates
% --> comparable Maps across sample
% 7. temproal filtering
%    7.1. remove offsets & linear trends
%    7.2. retainign frequencies below 0.08 Hz


% SFC
% 1. individual association matrices (pearson R voxel-wise)
%   Therefore, BOLD images transformed to N (brain) x M (time) matrix)
%   Get N x N association matrix
% 2. Binarize association matrix using 0.001 FDR correction only positive
% correlations
% 3.
function [] = SFC_voxelwise()

%% Prep
spm_dir = '...\spm12\';
addpath(genpath(spm_dir));

study_dir = '...\';
T1dir = '...\minimally_Processed\Anatomy\'; %NIFTI Data manually reoriented to MNI space

niiTools = '...\_Tools\NIfTI_20140122'; %REQUIREMENTS: NifTI Toolbox
addpath(niiTools);
addpath(genpath('...\CATO_3.1.2_unix')); %REQUIREMENTS: CATO Toolbox

primRoiDir = '...\Size_8_mm_adjusted'; %Seeds, created manually in MarsBar
resampled = true;

subjectNames= {'HFG_121'};

subject_num = [42];

nImages = 475;

normVoxel = 8;

repetitionTimeMsec = 1000;
minRepetitionTime = 100;
bandpass_filter_frequencies = [0.01,0.1];

q_fdr=0.001; % q for the fdr
num_steps= 7; % number of steps for the SFC

multi_seed = true; % Do you want to calcualte the mutli_seed variant of SFC

% Directory of MNI atlases
msbp_atlas_dir = '...\02_Ji_atlas\msbp_anat\';

% FUNCTIONS
tic
  give_head()
%data is already coregistered from CATO (using FSL)
  segmentation()
  normalisation_epis()
  normalized_brain_mask()
  create_voxelwise_matrix()
  roi_mni_space()
  SFC_calculation_JS()
  if multi_seed, SFC_calculation_multi_JS(), end
  GLM_2nd_level()
  export_node_SFC()
  normalize_anat() % and prep for Bids
  export_indiv_SFC()
toc

%% Add header information to the rs-files
    function [] = give_head()
        for i = 1:numel(subjectNames)
            disp(['Preparing: ', subjectNames{i}]);
            subjectDir = fullfile(study_dir, subjectNames{i});
            
            % make folders & move files
            sfcDir = fullfile(subjectDir, 'SFC_analysis');
            if ~exist(sfcDir, 'dir')
                mkdir(sfcDir);
                if~exist(fullfile(sfcDir, 'Anatomy'),'dir'), mkdir(fullfile(sfcDir, 'Anatomy')), end
                if~exist(fullfile(sfcDir, 'Epis'),'dir'), mkdir(fullfile(sfcDir, 'Epis')), end
            end
            
            % Functional % both in .../minimally_Processed/resting_state_fMRI
            regFile = fullfile(subjectDir, [subjectNames{i}, '_fmri_24pXaCompCorXVolterra_spikeReg.nii']);
            procFile = fullfile(subjectDir, 'fMRI_processed', [subjectNames{i},'_fmri.nii.gz']);
            
            gunzip(procFile, fullfile(sfcDir, 'Epis'));
            
            % Anatomy
            sourceFile = [T1dir, subjectNames{i},'_D2.nii'];
            destFile = [sfcDir,'\Anatomy\', subjectNames{i}, '_anat.nii'];
            copyfile(sourceFile, destFile);
            
            % create the real functional file % V01
            procNii = load_untouch_nii([fullfile(sfcDir, 'Epis'),filesep, subjectNames{i},'_fmri.nii']);
            regNii = load_untouch_nii(regFile);
            
            procNii_copy = procNii;
            
            procNii_copy.img = regNii.img;
            
            outpath = [sfcDir, '\Epis\', subjectNames{i}, '_fmriResid.nii'];
            
            save_untouch_nii(procNii_copy, outpath);
            
            disp(['DONE ', num2str(i),'/', num2str(numel(subjectNames))]);
        end
    end

%% Segmentation

    function [] = segmentation()
        
        for i = 1:numel(subjectNames)
            
            disp('Segmentation for: ')
            disp(subjectNames{i})
            
            subjectDir = fullfile(study_dir, subjectNames{i});
            
            % make folders & move files
            sfcDir = fullfile(subjectDir, 'SFC_analysis');
            
            % load anatomy
            highres_dir = fullfile(sfcDir, 'Anatomy');
            find_highres = dir(fullfile(highres_dir, '*_anat.nii'));
            highres_file = find_highres.name;
            
            % Batch commands (lots of it is just the parameters for the
            % different materials
            jobs{1}.spm.spatial.preproc.channel.vols = {[highres_dir filesep highres_file]};
            jobs{1}.spm.spatial.preproc.channel.biasreg = 0.001;
            jobs{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            jobs{1}.spm.spatial.preproc.channel.write = [0 1];
            jobs{1}.spm.spatial.preproc.tissue(1).tpm = {[spm_dir,'tpm',filesep,'TPM.nii,1']}; % grey matter
            jobs{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
            jobs{1}.spm.spatial.preproc.tissue(1).native = [1 1];
            jobs{1}.spm.spatial.preproc.tissue(1).warped = [1 1];
            jobs{1}.spm.spatial.preproc.tissue(2).tpm = {[spm_dir,'tpm',filesep,'TPM.nii,2']}; % white matter
            jobs{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
            jobs{1}.spm.spatial.preproc.tissue(2).native = [1 1];
            jobs{1}.spm.spatial.preproc.tissue(2).warped = [1 1];
            jobs{1}.spm.spatial.preproc.tissue(3).tpm = {[spm_dir,'tpm',filesep,'TPM.nii,3']}; % CSF
            jobs{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
            jobs{1}.spm.spatial.preproc.tissue(3).native = [1 0];
            jobs{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
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
            jobs{1}.spm.spatial.preproc.warp.write = [0 1]; % foward
            % RUN
            spm_jobman('run', jobs);
            clear jobs
        end
        fprintf('Segmentation:\n Done.\n')
    end
%% SPM Vanilla Normalization
%% Normalisation Epis

    function [] = normalisation_epis()
        
        for i =  1:numel(subjectNames)
            
            disp('Writing the normalised Epis for: ')
            disp(subjectNames{i})
            
            sfcDir = fullfile(study_dir,subjectNames{i} ,'SFC_analysis');
            
            epi_dir = fullfile(sfcDir,'Epis');
            highres_dir = fullfile(sfcDir,'Anatomy');
            find_deform = dir(fullfile(highres_dir,'y_*'));
            deform_file = find_deform.name;
            % Einladen der Epis mit a* (Seit Slice Timing ist nichts passiert)und deformation field ('y_*')
            [func_filenames_o, ~] = spm_select('ExtList',epi_dir,'Resid*.',1:nImages);
            files_for_norm = strcat(epi_dir,filesep, func_filenames_o(1:nImages,:));
            
            % Batch commands
            jobs{1}.spm.spatial.normalise.write.subj.def = {[highres_dir filesep deform_file]};
            jobs{1}.spm.spatial.normalise.write.subj.resample = cellstr(files_for_norm);
            jobs{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                78 76 85];
            % Voxel size/dimensions
            jobs{1}.spm.spatial.normalise.write.woptions.vox = [normVoxel normVoxel normVoxel]; % enter voxel size of normalised epis
            jobs{1}.spm.spatial.normalise.write.woptions.interp = 4;
            jobs{1}.spm.spatial.normalise.write.woptions.prefix = ['w',num2str(normVoxel),'_'];
            % RUN
            spm_jobman('run', jobs);
            clear jobs
        end
        fprintf('Normalisation of Epis:\n Done.\n')
    end

%% Create a normalized Whole brain mask...

% take the normalice wc1 and wc2 of all subject and create an "individual"
% normalized brain mask, then average that over all subject
% lastly, coregister brain mask with functional space

    function [] = normalized_brain_mask()
        
        normMasks = cell(numel(subjectNames), 1);
        expression = '((';
        % Step 1: "invidual noramlized brain masks"
        for i = 1:numel(subjectNames)
            disp(['Brain Masking for: ', subjectNames{i}]);
            anatDir = fullfile(study_dir, subjectNames{i}, 'SFC_analysis', 'Anatomy');
            grayMatNorm = fullfile(anatDir, ['wc1',subjectNames{i}, '_anat.nii']);
            whiteMatNorm = fullfile(anatDir, ['wc2',subjectNames{i}, '_anat.nii']);
            
            jobs{1}.spm.util.imcalc.input = {grayMatNorm; whiteMatNorm};
            jobs{1}.spm.util.imcalc.output = ['wBrainMask_',subjectNames{i},'.nii'];
            jobs{1}.spm.util.imcalc.outdir = {anatDir};
            jobs{1}.spm.util.imcalc.expression = '(i1 + i2) > 0.2';
            jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            jobs{1}.spm.util.imcalc.options.dmtx = 0;
            jobs{1}.spm.util.imcalc.options.mask = 0;
            jobs{1}.spm.util.imcalc.options.interp = 1;
            jobs{1}.spm.util.imcalc.options.dtype = 4;
            
            spm_jobman('run', jobs);
            clear jobs
            
            jobs{1}.spm.spatial.coreg.write.ref = {fullfile(study_dir, subjectNames{i},'SFC_analysis', 'Epis', ['w',num2str(normVoxel),'_', subjectNames{i},'_fmriResid.nii,15'])};
            jobs{1}.spm.spatial.coreg.write.source = {fullfile(anatDir, ['wBrainMask_',subjectNames{i},'.nii'])};
            jobs{1}.spm.spatial.coreg.write.roptions.interp = 4;
            jobs{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            jobs{1}.spm.spatial.coreg.write.roptions.mask = 0;
            jobs{1}.spm.spatial.coreg.write.roptions.prefix = ['r',num2str(normVoxel),'_'];
            spm_jobman('run', jobs);
            clear jobs
            
            normMasks{i,1} = fullfile(anatDir, ['r',num2str(normVoxel),'_','wBrainMask_',subjectNames{i},'.nii']);
            normMasks{i,1} = fullfile(anatDir, ['wBrainMask_',subjectNames{i},'.nii']);
            
            expression = [expression,'i', num2str(i),'+'];
            disp(['DONE ', num2str(i),'/', num2str(numel(subjectNames))]);
            
        end
        
        %Step 2: create whole brain mask for all subjects
        expression = [expression(1:end-1),')/',num2str(numel(subjectNames)),')> 0.2'];
        
        jobs{1}.spm.util.imcalc.input = normMasks;
        jobs{1}.spm.util.imcalc.output = 'wBrainMask.nii';
        jobs{1}.spm.util.imcalc.outdir = {study_dir};
        jobs{1}.spm.util.imcalc.expression = expression;
        jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        jobs{1}.spm.util.imcalc.options.dmtx = 0;
        jobs{1}.spm.util.imcalc.options.mask = 0;
        jobs{1}.spm.util.imcalc.options.interp = 1;
        jobs{1}.spm.util.imcalc.options.dtype = 512;
        disp('Brain Masking whole data set');
        spm_jobman('run', jobs);
        clear jobs
        
        % Step 3: coregistier (but test first if necessary)
        jobs{1}.spm.spatial.coreg.write.ref = {fullfile(study_dir, subjectNames{1},'SFC_analysis', 'Epis', ['w',num2str(normVoxel),'_', subjectNames{1},'_fmriResid.nii,15'])};
        jobs{1}.spm.spatial.coreg.write.source = {fullfile(study_dir, 'wBrainMask.nii')};
        jobs{1}.spm.spatial.coreg.write.roptions.interp = 4;
        jobs{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
        jobs{1}.spm.spatial.coreg.write.roptions.mask = 0;
        jobs{1}.spm.spatial.coreg.write.roptions.prefix = ['r',num2str(normVoxel),'_'];
        spm_jobman('run', jobs);
        clear jobs
        
    end

%% Partially load normalized nifits using whole brain mask
% BAsed on reconstruction of timeTeries_network
    function[] = create_voxelwise_matrix()
        
        brainmask = load_nii(fullfile(study_dir, ['r',num2str(normVoxel),'_','wBrainMask.nii']));
        brainmask = brainmask.img(:) > 0;
        missings = [];
        for i = 1:numel(subjectNames)
            disp(['Voxelwise Matrix for: ', subjectNames{i}]);
            
            epiDir = fullfile(study_dir, subjectNames{i}, 'SFC_analysis','Epis');
            niiFile = [epiDir, filesep, 'w',num2str(normVoxel),'_',subjectNames{i}, '_fmriResid.nii'];
            
            fmri = load_nifti_partially(niiFile, brainmask);
            signalIntensities = single(fmri.partialvol);
            missings = [missings;find(isnan(signalIntensities(:,1)))];
            % I don't think this works for voxelwise
            %             prevalenceMask = mean(signalIntensities ~= 0, 2) >= 0.9;
            %             signalIntensities = signalIntensities(prevalenceMask, :);
            %             brainPrevalenceMask = false(size(brainMask));
            %             brainPrevalenceMask(brainMask) = prevalenceMask;
            
            %% BANDPASS FILTER
            assert(repetitionTimeMsec >= minRepetitionTime, ...
                ['Repetition time (%g msec) reported in fmriProcessedFile', ...
                ' is smaller than the minRepetitionTime (%g msec). ', ...
                'Transform repetition time to milliseconds ', ...
                'or adjust minRepetitionTime-parameter'], repetitionTimeMsec, minRepetitionTime);
            repetitionTimeSec = repetitionTimeMsec/1000;
            
            [filter_b, filter_a] = butter(2, 2*repetitionTimeSec*bandpass_filter_frequencies);
            
            % Use for-loop to avoid memory issues from having double. (1sec difference)
            % filteredSignal = filtfilt(filter_b, filter_a, double(signalIntensities(selectedVoxels, :)));
            filteredSignal = zeros(size(signalIntensities), 'single');
            for j = 1:size(signalIntensities, 1)
                filteredSignal(j,:) = filtfilt(filter_b, filter_a, ...
                    double(signalIntensities(j, :)));
            end
            
            %% Correlaton matrix
            [connectivity, pValues] = corr(filteredSignal');
            connectivity = connectivity .* ~eye(length(connectivity));
            pValues = pValues .* ~eye(length(connectivity));
            
            outFile = fullfile(study_dir, subjectNames{i}, 'SFC_analysis',[subjectNames{i},'_norm',num2str(normVoxel),'_voxWiseCon.mat']);
            
            save(outFile,'connectivity', 'pValues','-v7.3');
            
        end
        miss = unique(missings);
        save(fullfile(study_dir, 'missings.mat'), 'miss'); % in .../02_UM/01_SFC/
    end


%% Define ROIS in MNI space/ in the matrix
% create them and them normalize them to same space then they should be the
% same position as all the rest

% Coordinates based on Sepulcre, Corrected with anatomy toolbox to fix
% double coordiantes in x direction (creates double seed)
% Visual: Left -13(-14), -78, 8; Right 9(10), -78, 8
% Auditory: Left -53(-54), -14, 8; Right 57(58), -14, 8 % but auditory is wierd in general
% Somatosensory: Left -42/-43, -29, 65; Right 38/39, -29, 65

% created with MarsBaR

    function [] = roi_mni_space()
        % Threshold
        primRois = dir(primRoiDir);
        primRois = {primRois.name};
        primRois = primRois(endsWith(primRois, 'nii'))';
        
        for ri = 1:numel(primRois)
            
            if ~resampled
                clear jobs;
                % Coregister to noramlized and resampled images
                jobs{1}.spm.spatial.coreg.write.ref = {fullfile(study_dir, subjectNames{1},'SFC_analysis', 'Epis', ['w',num2str(normVoxel),'_', subjectNames{1},'_fmriResid.nii,15'])};
                jobs{1}.spm.spatial.coreg.write.source = {fullfile(primRoiDir, primRois{ri})};
                jobs{1}.spm.spatial.coreg.write.roptions.interp = 4;
                jobs{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                jobs{1}.spm.spatial.coreg.write.roptions.mask = 0;
                jobs{1}.spm.spatial.coreg.write.roptions.prefix = ['r',num2str(normVoxel),'_'];
                spm_jobman('run', jobs);
                prefix = ['r',num2str(normVoxel),'_'];
            else
                prefix = [];
            end
            % binarize again in new reference space
            clear jobs
            
            jobs{1}.spm.util.imcalc.input = {fullfile(primRoiDir, [prefix,primRois{ri}])};
            jobs{1}.spm.util.imcalc.output = ['s','r',num2str(normVoxel),'_',primRois{ri}] ;
            jobs{1}.spm.util.imcalc.outdir = {primRoiDir};
            jobs{1}.spm.util.imcalc.expression = 'i1 > 0.2';
            jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            jobs{1}.spm.util.imcalc.options.dmtx = 0;
            jobs{1}.spm.util.imcalc.options.mask = 0;
            jobs{1}.spm.util.imcalc.options.interp = 1;
            jobs{1}.spm.util.imcalc.options.dtype = 4;
            
            spm_jobman('run', jobs);
            clear jobs;
        end
        
    end

%% Stepwise functional connectivity based on the MAil by Jorge Sepulce (2023.12.19)
    function [] = SFC_calculation_JS()
        addpath(genpath(fullfile(study_dir,'02_UM','01_SFC','ExampleCentralityandStepwiseConnectivity')));
        
        addpath('...\'); % for CATO and NIFTI Tools
        
        % load missings to find a non match
        missings = load(fullfile(study_dir, 'missings.mat'));
        
        brainmask_nii = load_nifti(fullfile(study_dir, ['r',num2str(normVoxel),'_','wBrainMask.nii']));
        brainmask_nii.scl_slope = 1;
        brainmask = brainmask_nii.vol(:) > 0;
        num_labels = 1:sum(brainmask);
        
        all_have = num_labels(~ismember(num_labels,missings.miss));
        
        %primRoiDir = DEFINITION NOW AT THE START
        primRois = dir(primRoiDir);
        primRois = {primRois.name};
        primRois = primRois(startsWith(primRois, ['sr',num2str(normVoxel)]))';
        
        for ri = 1:numel(primRois)
            primRoiFile = fullfile(primRoiDir, primRois{ri});
            roi = load_nifti_partially(primRoiFile, brainmask);
            name = strsplit(primRois{ri}(5:end),'_'); name = name{1};
            %primarySeeds.(name) = num_labels(roi.partialvol == 1);
            primarySeeds.(name) = roi.partialvol == 1;
            disp(sum(primarySeeds.(name)));
        end
        primNames = fieldnames(primarySeeds);
        
        for i = 1:numel(subjectNames)
            disp('Calculating SFC for:');
            disp(subjectNames{i});
            
            conFileName = fullfile(study_dir, subjectNames{i}, 'SFC_analysis',[subjectNames{i},'_norm',num2str(normVoxel),'_voxWiseCon.mat']);
            conMat = load(conFileName);
            
            fc = conMat.connectivity;
            if sum(isnan(fc(all_have(1),:))) > 0
                fc(isnan(fc(all_have(1),:)),:) = 0;
                fc(:,isnan(fc(all_have(1),:))) = 0;
            end
            fcFish = 0.5 * log( ( 1 + fc)./( 1 - fc) ); %fisher (optional)
            % Diez & Sepulcre, 2018: Finally, we applied a variance-stabilizing
            % transformation (Fisher's transformation) to all correlation coefficients of association connectivity matrices as a final step before our graph theory-based analysis
            fcPos = fcFish.* (fcFish > 0);
           
            % FWE Corrected
            indx=find(fcPos >0);
            indx_signif = conMat.pValues(indx) < (q_fdr/length(indx));
            fwe_matrix=zeros(size(conMat.pValues));
            fwe_matrix(indx(indx_signif))=1;
            fwe_matrix=uint8(fwe_matrix);
            adj=fcPos.*double(fwe_matrix);
            
            adjNorm=(adj - min(adj(:)))./ (max(adj(:)) - min(adj(:))); %normalize (optional)
            [m ,~]=size(adjNorm);
            adjNorm(1:(m+1):end)=0; %diagonal to zero
            
            %CENTRALITY & STEPWISE
            degree_centrality=sum(adjNorm,2);
            
            for ri = 1:numel(primRois)
                adjp=adjNorm;
                step_seed = cell(num_steps,1);
                seed_indx = find(primarySeeds.(primNames{ri}) > 0);
                step1_seed=adjp(seed_indx,:);
                step_seed{1} = step1_seed;
                
                for s=2:num_steps
                    adj_new=adjp*adjp; % for exponential increase, multiply by output; for arithmetic increase multiply always by step1
                    adj_new=(adj_new - min(adj_new(:)))./ (max(adj_new(:)) - min(adj_new(:))); %normalize (optional)
                    adj_new(1:(m+1):end)=0; %diagonal to zero
                    adjp=adj_new;
                    %                     eval(sprintf('step%s=adjp;',num2str(s)));
                    %                     eval(sprintf('step%s_seed=adjp(seed_indx,:);',num2str(s)));
                    step_seed{s} = adjp(seed_indx,:);
                end
                if length(seed_indx) > 1
                    step_seed_final = cellfun(@(x) mean(x,1), step_seed, 'UniformOutput', false);
                    step_seed_final = cell2mat(step_seed_final)';
                elseif length(seed_indx) == 1
                    step_seed_final = cell2mat(step_seed)';
                end
                
                testCorr = [];
                for ii = 1:(num_steps-1)
                    testCorr = [testCorr,corr(step_seed_final(:,ii),step_seed_final(:,ii+1))];
                end
                plot(testCorr);% testCorr
                
                outDir = fullfile(study_dir, subjectNames{i}, 'SFC_analysis',primNames{ri});
                if ~exist(outDir,'dir'), mkdir(outDir), end
                hdr0 = load_untouch_nii(fullfile(study_dir, ['r',num2str(normVoxel),'_','wBrainMask.nii']));
                hdr0.hdr.dime.scl_slope = 1;
                hdr0.hdr.dime.datatype = 16; hdr0.hdr.dime.bitpix = 32;
                img=zeros(size(hdr0.img));
                img(brainmask) = degree_centrality;
                hdr0.hdr.dime.glmax = max(degree_centrality); hdr0.hdr.dime.glmin = min(degree_centrality);
                hdr0.img=img;%./(1.5*10^-5);
                outName = [outDir,filesep, subjectNames{i}, '_degree_centrality_',primNames{ri},'.nii'];
                save_untouch_nii(hdr0,outName);
                    clear img
                for s = 1:num_steps
                    img=zeros(size(hdr0.img));
                    img(brainmask) = step_seed_final(:,s);
                    hdr0.hdr.dime.glmax = max(step_seed_final(:,s)); hdr0.hdr.dime.glmin = min(step_seed_final(:,s));
                    hdr0.img=img;
                    outName = [outDir,filesep, subjectNames{i},'_SFC_',num2str(s),'_',primNames{ri},'.nii'];
                    save_untouch_nii(hdr0,outName);
                    clear img
                end
            end
        end
    end

%% GLM 2nd Level with SPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% General Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [] = GLM_2nd_level()
        % Definition von Ordnern
        analysis_dir = '_SFC_Second_Level';
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Second Level
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        spm('defaults','FMRI');
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %% one sample t-tests
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ~exist(fullfile(study_dir,analysis_dir),'dir')
            mkdir(fullfile(study_dir,analysis_dir)); % Erstellen des 2nd level Analyse-Ordners
        end
   
        
        %% Specify 2nd-level
        
        primRois = dir(primRoiDir);
        primRois = {primRois.name};
        primRois = primRois(startsWith(primRois, ['sr',num2str(normVoxel)]))';
        primNames = cellfun(@(x) split(x,'_'),primRois,'UniformOutput', false);
        primNames = [primNames{:}];
        primNames = primNames(2,:);
        
        if multi_seed
            primNames = [primNames, {'Multi_seed'}];
        end
        
        for ri = 1:numel(primNames)
            conDir = fullfile(study_dir,analysis_dir, primNames{ri}); % hier Ordnernamen eingeben %!!!!!!!
            if ~exist(conDir, 'dir')
                mkdir(conDir);
            end
            
            clear jobs;
            
            jobs{1}.spm.stats.factorial_design.dir = {conDir};
                        
            A = cell(numel(subjectNames),1); %clear matlabbatch;
            for k=1:numel(subjectNames)
                fileName = [subjectNames{k},'_SFC_7_',primNames{ri},'.nii'];
                % matlabbatch{1}.spm.util.disp.data = {fullfile(study_dir,subjectNames{k},'SFC_analysis',primNames{ri},fileName)};
                % spm_jobman('run', matlabbatch);
                % clear matlabbatch;
                A{k} = fullfile(study_dir,subjectNames{k},'SFC_analysis',primNames{ri},fileName); % hier con_image-Nummer eingeben %!!!!!!!
            end
            jobs{1}.spm.stats.factorial_design.des.t1.scans = cellstr(strvcat(A));
            
            jobs{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); % Covariates
            jobs{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {}); % Multiple covaraites
            jobs{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;% Threshold masking
            jobs{1}.spm.stats.factorial_design.masking.im = 0; % Implicit Mask (Yes)
            jobs{1}.spm.stats.factorial_design.masking.em = {fullfile(study_dir,'r8_wBrainMask.nii') }; % Explicit Mask
            jobs{1}.spm.stats.factorial_design.globalc.g_omit = 1;% Global calculation
            jobs{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1; % Global normalisation
            jobs{1}.spm.stats.factorial_design.globalm.glonorm = 1; % Normalisation
            
            spm_jobman('run',jobs); %RUN
            clear jobs
            
            
            
            %% Estimate (2nd-level)
            
            jobs{1}.spm.stats.fmri_est.spmmat = {fullfile(conDir,'SPM.mat')};
            jobs{1}.spm.stats.fmri_est.write_residuals = 1;
            jobs{1}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',jobs) %RUN
            clear jobs;
            
            %% Contrasts (via Contrast Manager)
            jobs{1}.spm.stats.con.spmmat = {fullfile(conDir,'SPM.mat')};
            
            % CONTRASTS
            jobs{1}.spm.stats.con.consess{1}.tcon.name = [primNames{ri},'_Step7_Plus'];%NAMEN AENDERN %!!!!!!!
            jobs{1}.spm.stats.con.consess{1}.tcon.convec = 1;% Kontrastgewicht eintragen %!!!!!!!
            jobs{1}.spm.stats.con.consess{2}.tcon.name = [primNames{ri},'_Step7_Minus'];%NAMEN AENDERN %!!!!!!!
            jobs{1}.spm.stats.con.consess{2}.tcon.convec = -1;% Kontrastgewicht eintragen %!!!!!!!
            
            spm_jobman('run',jobs); %RUN
            clear jobs
            
            %% Convert t maps to cohens d
            spmT_file01 = fullfile(conDir,'spmT_0001.nii');
            % spmT_file02 = fullfile(conDir,'spmT_0002.nii');
            spmMat = load(fullfile(conDir,'SPM.mat'));
            df = spmMat.SPM.xX.erdf;
            jobs{1}.spm.util.imcalc.input = {spmT_file01};
            jobs{1}.spm.util.imcalc.output = ['spmP_',primNames{ri},'_0001.nii'];
            jobs{1}.spm.util.imcalc.outdir = {conDir};
            jobs{1}.spm.util.imcalc.expression = ['(i1*2)/(sqrt(',num2str(df),'))'];
            jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            jobs{1}.spm.util.imcalc.options.dmtx = 0;
            jobs{1}.spm.util.imcalc.options.mask = 0;
            jobs{1}.spm.util.imcalc.options.interp = 1;
            jobs{1}.spm.util.imcalc.options.dtype = 4;
            
            spm_jobman('run',jobs); %RUN
            clear jobs
            
        end
    end

%% Calcualte SFC for all Seeds (meaned)

    function [] = SFC_calculation_multi_JS()
        addpath(genpath(fullfile(study_dir,'_Tools','SFC_Sepulcre','ExampleCentralityandStepwiseConnectivity')));
        
        addpath('...\_Tools');
        
        % load missings to find a non match
        missings = load(fullfile(study_dir, 'missings.mat'));
        
        brainmask_nii = load_nifti(fullfile(study_dir, ['r',num2str(normVoxel),'_','wBrainMask.nii']));
        brainmask_nii.scl_slope = 1;
        brainmask = brainmask_nii.vol(:) > 0;
        num_labels = 1:sum(brainmask);
        
        all_have = num_labels(~ismember(num_labels,missings.miss));
        
        primRois = dir(primRoiDir);
        primRois = {primRois.name};
        primRois = primRois(startsWith(primRois, ['sr',num2str(normVoxel)]))';
        
        for ri = 1:numel(primRois)
            primRoiFile = fullfile(primRoiDir, primRois{ri});
            roi = load_nifti_partially(primRoiFile, brainmask);
            name = strsplit(primRois{ri}(5:end),'_'); name = name{1};
            %primarySeeds.(name) = num_labels(roi.partialvol == 1);
            primarySeeds.(name) = roi.partialvol == 1;
            disp(sum(primarySeeds.(name)));
        end
        primNames = fieldnames(primarySeeds);
        
        for i = 1:numel(subjectNames)
            disp('Calculating SFC for:');
            disp(subjectNames{i});
            
            conFileName = fullfile(study_dir, subjectNames{i}, 'SFC_analysis',[subjectNames{i},'_norm',num2str(normVoxel),'_voxWiseCon.mat']);
            conMat = load(conFileName);
            
            fc = conMat.connectivity;
            if sum(isnan(fc(all_have(1),:))) > 0
                fc(isnan(fc(all_have(1),:)),:) = 0;
                fc(:,isnan(fc(all_have(1),:))) = 0;
            end
            fcFish = 0.5 * log( ( 1 + fc)./( 1 - fc) ); %fisher (optional)
            % Diez & Sepulcre, 2018: Finally, we applied a variance-stabilizing
            % transformation (Fisher's transformation) to all correlation coefficients of association connectivity matrices as a final step before our graph theory-based analysis
            fcPos = fcFish.* (fcFish > 0);
            
            % FWE Corrected
            indx=find(fcPos >0);
            indx_signif = conMat.pValues(indx) < (q_fdr/length(indx));
            fwe_matrix=zeros(size(conMat.pValues));
            fwe_matrix(indx(indx_signif))=1;
            fwe_matrix=uint8(fwe_matrix);
            adj=fcPos.*double(fwe_matrix);
            
            adjNorm=(adj - min(adj(:)))./ (max(adj(:)) - min(adj(:))); %normalize (optional)
            [m ,~]=size(adjNorm);
            adjNorm(1:(m+1):end)=0; %diagonal to zero
            
            %CENTRALITY & STEPWISE
            degree_centrality=sum(adjNorm,2);
            
            seed_indx = [];
            for ri = 1:numel(primRois)
                seed_indx = [seed_indx; find(primarySeeds.(primNames{ri}) > 0)];
            end
            
            step_seed = cell(num_steps,1);
            
            adjp=adjNorm;
            step_seed{1} = adjp(seed_indx,:);
            
            for s=2:num_steps
                adj_new=adjp*adjp; % for exponential increase, multiply by output; for arithmetic increase multiply always by step1
                adj_new=(adj_new - min(adj_new(:)))./ (max(adj_new(:)) - min(adj_new(:))); %normalize (optional)
                adj_new(1:(m+1):end)=0; %diagonal to zero
                adjp=adj_new;
                % eval(sprintf('step%s=adjp;',num2str(s)));
                % eval(sprintf('step%s_seed=adjp(seed_indx,:);',num2str(s)));
                step_seed{s} = adjp(seed_indx,:);
            end
            
            step_seed_final = cellfun(@(x) mean(x,1), step_seed, 'UniformOutput', false);
            step_seed_final = cell2mat(step_seed_final)';
            
            testCorr = [];
            for ii = 1:(num_steps-1)
                testCorr = [testCorr,corr(step_seed_final(:,ii),step_seed_final(:,ii+1))];
            end
            plot(testCorr);% testCorr
            
            outDir = fullfile(study_dir, subjectNames{i}, 'SFC_analysis','Multi_seed');
            if ~exist(outDir,'dir'), mkdir(outDir), end
            hdr0 = load_untouch_nii(fullfile(study_dir, ['r',num2str(normVoxel),'_','wBrainMask.nii']));
            hdr0.hdr.dime.scl_slope = 1;
            hdr0.hdr.dime.datatype = 16; hdr0.hdr.dime.bitpix = 32;
            img=zeros(size(hdr0.img));
            img(brainmask) = degree_centrality;
            hdr0.hdr.dime.glmax = max(degree_centrality); hdr0.hdr.dime.glmin = min(degree_centrality);
            hdr0.img=img;%./(1.5*10^-5);
            outName = [outDir,filesep, subjectNames{i}, '_degree_centrality_Multi_seed.nii'];
            save_untouch_nii(hdr0,outName);
             clear img
            for s = 1:num_steps
                img=zeros(size(hdr0.img));
                img(brainmask) = step_seed_final(:,s);
                hdr0.hdr.dime.glmax = max(step_seed_final(:,s)); hdr0.hdr.dime.glmin = min(step_seed_final(:,s));
                hdr0.img=img;%./(1.5*10^-5);
                outName = [outDir,filesep, subjectNames{i},'_SFC_',num2str(s),'_Multi_seed.nii'];
                save_untouch_nii(hdr0,outName);
                clear img
            end
            
        end
    end

%% Export SFC values (per subject) on Node-Level
% Scale2 = lausanne120
% Scale3 = lausanne250

% num_steps = 7;

    function [] = export_node_SFC()
        
        % parcellations and scales
        parc = {'scale1', 'scale2', 'scale3', 'scale4';...
            'aparc', 'lausanne120', 'lausanne250', 'lausanne500'};
        
        % used Seeds and multi-seed
        seeds = {'LeftAuditory', 'LeftSomatosensory', 'LeftVisual', 'Multi_seed', ...
            'RightAuditory',  'RightSomatosensory', 'RightVisual'};
        
        colNames = [{'Subject', 'Template', 'Region'},seeds]; % column names for output table
        
        for p = 3 %1:size(parc,2)
            % load translator: what area has what number in MSBP
            translator = readtable([msbp_atlas_dir, 'transS',parc{1,p}(2:end), '.csv']);
            
            atlas_nii = [msbp_atlas_dir, 'sub-01_label-L2018_desc-',parc{1,p},'_atlas.nii'];
            atlas = load_untouch_nii(atlas_nii);
            
            for i = 1:numel(subjectNames) % 1:numel(subjectNames)
                disp(subjectNames{i});
                output = cell(size(translator,1),numel(seeds)+3); % create Output for all seeds this subject
                subject_Dir = fullfile(study_dir, subjectNames{i},'SFC_analysis\');
                
                % All step 7 files, num_steps = 7
                sfc_imgs = cell(numel(seeds),1);
                for s = 1:numel(seeds)
                    sfc_imgs{s} = [subject_Dir, seeds{s}, filesep, subjectNames{i}, '_SFC_', num2str(num_steps),'_', seeds{s},'.nii'];
                end
                
                % USE COREGISTER (RESLICE) to coregister sfc maps to atlas
                % size
                clear jobs;
                % bring individual seven Map to size of atlas...
                jobs{1}.spm.spatial.coreg.write.ref = {atlas_nii};
                jobs{1}.spm.spatial.coreg.write.source = sfc_imgs;
                jobs{1}.spm.spatial.coreg.write.roptions.interp = 4;
                jobs{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                jobs{1}.spm.spatial.coreg.write.roptions.mask = 0;
                jobs{1}.spm.spatial.coreg.write.roptions.prefix = 'atl_';
                % ALREADY GETS FLIPPED, YEY
                
                spm_jobman('run', jobs); %RUN
                clear jobs;
                
                if ~exist(fullfile(subject_Dir, '_SFC_resized'),'dir')
                    mkdir(fullfile(subject_Dir, '_SFC_resized'));
                end
                
                for s = 1:numel(seeds) % for each seed
                    disp(seeds{s});
                    % move the resized sfc_map
                    sfc_resized = [subject_Dir, seeds{s}, filesep,'atl_' ,subjectNames{i}, '_SFC_', num2str(num_steps),'_', seeds{s},'.nii'];
                    movefile(sfc_resized, fullfile(subject_Dir, '_SFC_resized'));
                    
                    % load the map
                    sfc_map = load_untouch_nii([subject_Dir, '_SFC_resized', filesep,'atl_' ,subjectNames{i}, '_SFC_', num2str(num_steps),'_', seeds{s},'.nii']);
                    
                    for r = 1:size(translator,1) % for each area
                        % disp(translator.region{r});
                        if s == 1
                            output{r,1} =  subjectNames{i};
                            output{r,2} =  parc{2,p};
                            output{r,3} =  translator.region{r};
                        end
                        % get the sfc values of that area
                        sfc_values = sfc_map.img(atlas.img == translator.msbpROI(r));
                        output{r,s+3} = mean(sfc_values, 'all', 'omitnan');
                        
                    end
                end
                
                % save out put for that subject
                outTab = cell2table(output, 'VariableNames', colNames);
                outName = [subject_Dir, subjectNames{i}, '_SFC_', num2str(num_steps),'_', parc{2,p},'.csv'];
                writetable(outTab, outName);
                
                disp(['DONE ', num2str(i),'/', num2str(numel(subjectNames))]);
            end % repeat for subjects
            
            % Second level: insteach of sfc maps use contrast maps
            disp('Second_Level');
            output = cell(size(translator,1),numel(seeds)+3); % create Output for all seeds this subject
            second_Dir = fullfile(study_dir, '_SFC_Second_Level\');
            
            % All step 7 files, num_steps = 7
            t_imgs = cell(numel(seeds),1);
            for s = 1:numel(seeds)
                t_imgs{s} = [second_Dir, seeds{s}, '\spmT_0001.nii'];
            end
            
            % USE COREGISTER (RESLICE) to coregister sfc maps to atlas
            % size
            clear jobs;
            % bring individual seven Map to size of atlas...
            jobs{1}.spm.spatial.coreg.write.ref = {atlas_nii};
            jobs{1}.spm.spatial.coreg.write.source = t_imgs;
            jobs{1}.spm.spatial.coreg.write.roptions.interp = 4;
            jobs{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            jobs{1}.spm.spatial.coreg.write.roptions.mask = 0;
            jobs{1}.spm.spatial.coreg.write.roptions.prefix = 'atl_';
            % ALREADY GETS FLIPPED, YEY
            
            spm_jobman('run', jobs); %RUN
            clear jobs;
            
            if ~exist(fullfile(second_Dir, '_SFC_resized'),'dir')
                mkdir(fullfile(second_Dir, '_SFC_resized'));
            end
            
            for s = 1:numel(seeds) % for each seed
                disp(seeds{s});
                % move the resized sfc_map
                t_resized = [second_Dir, seeds{s}, filesep, 'atl_spmT_0001.nii'];
                movefile(t_resized, [second_Dir, '_SFC_resized\',seeds{s}, '_Step_', num2str(num_steps), '_spmT_0001.nii']);
                
                % load the map
                second_map = load_untouch_nii([second_Dir, '_SFC_resized\',seeds{s}, '_Step_', num2str(num_steps), '_spmT_0001.nii']);
                
                for r = 1:size(translator,1) % for each area
                    % disp(translator.region{r});
                    if s == 1
                        output{r,1} =  'Second_Level';
                        output{r,2} =  parc{2,p};
                        output{r,3} =  translator.region{r};
                    end
                    % get the sfc values of that area
                    t_values = second_map.img(atlas.img == translator.msbpROI(r));
                    output{r,s+3} = mean(t_values, 'all', 'omitnan');
                    
                end
            end
            
            % save out put for that subject
            outTab = cell2table(output, 'VariableNames', colNames);
            outName = [second_Dir, 'Second_Level_tValues_Step_', num2str(num_steps),'_', parc{2,p},'.csv'];
            writetable(outTab, outName);
            
        end % repeat for templates
    end

%% NORMALIZE ANATOMY FOR PARCELLATION
    function [] = normalize_anat()
       
       newDir = '...\'; 
       bidsDir = [newDir, 'BidsDirSFC\'];
        
       for i =  1:numel(subjectNames)
            
            disp('Writing the normalised Anatomy for: ')
            disp(subjectNames{i})
            
            sfcDir = fullfile(study_dir,subjectNames{i} ,'SFC_analysis');
            
            highres_dir = fullfile(sfcDir,'Anatomy');
            find_deform = dir(fullfile(highres_dir,'y_*'));
            deform_file = find_deform.name;
            anat_file = fullfile(highres_dir, [subjectNames{i},'_anat.nii']); 
            
            % Batch commands
            jobs{1}.spm.spatial.normalise.write.subj.def = {[highres_dir filesep deform_file]};
            jobs{1}.spm.spatial.normalise.write.subj.resample = {anat_file};
            jobs{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                78 76 85];
            % Voxel size/dimensions
            jobs{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1]; % enter voxel size of normalised epis
            jobs{1}.spm.spatial.normalise.write.woptions.interp = 4;
            jobs{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
            % RUN
            spm_jobman('run', jobs);
            clear jobs
            
            %Make Bids direcotry
            subNum = sprintf('%02d',subject_num(i));
            subBids = [bidsDir, 'sub-', subNum, '/anat/'];
            if ~exist(subBids,'dir')
                mkdir(subBids); 
            end
            gzip(fullfile(highres_dir, ['w',subjectNames{i},'_anat.nii']));     
            movefile(fullfile(highres_dir, ['w',subjectNames{i},'_anat.nii.gz']),...
                [subBids, 'sub-',subNum, '_T1w.nii.gz']);
            copyfile(fullfile(newDir, 'sub-01_T1w.json'),...
                [subBids, 'sub-',subNum, '_T1w.json']);
            
            
        end
        fprintf('Normalisation of Anat:\n Done.\n')
    end

%% Export SFC values (per subject) on Node-Level
% Scale2 = lausanne120
% Scale3 = lausanne250

% num_steps = 7;

    function [] = export_indiv_SFC()
        
        % parcellations and scales
        parc = {'scale1', 'scale2', 'scale3', 'scale4';...
            'aparc', 'lausanne120', 'lausanne250', 'lausanne500'};
        
        % used Seeds and multi-seed
        seeds = {'LeftAuditory', 'LeftSomatosensory', 'LeftVisual', 'Multi_seed', ...
            'RightAuditory',  'RightSomatosensory', 'RightVisual'};
        
        colNames = [{'Subject', 'Template', 'Region'},seeds]; % column names for output table
        
        newDir = 'O:\07_WiSe_2023_24\'; 
        bidsDir = [newDir, 'BidsDirSFC\'];
        
        for p = 3 %1:size(parc,2)
            % load translator: what area has what number in MSBP
            translator = readtable([msbp_atlas_dir, 'transS',parc{1,p}(2:end), '.csv']);           
            for i = 39:numel(subjectNames) % 1:numel(subjectNames)
                disp(subjectNames{i});
                subNum = sprintf('%02d',subject_num(i));
                atlas_nii = fullfile(bidsDir,'derivatives', 'cmp',['sub-',subNum],'anat',['sub-',subNum,'_label-L2018_desc-',parc{1,p},'_atlas.nii.gz']);
                subject_Dir = fullfile(study_dir, subjectNames{i},'SFC_analysis\');
                gunzip(atlas_nii,subject_Dir)
                output = cell(size(translator,1),numel(seeds)+3); % create Output for all seeds this subject
                
                atlasMoved = [subject_Dir, filesep, 'sub-',subNum,'_label-L2018_desc-',parc{1,p},'_atlas.nii']; 
                atlas = load_untouch_nii(atlasMoved); 
                
                % All step 7 files, num_steps = 7
                sfc_imgs = cell(numel(seeds),1);
                for s = 1:numel(seeds)
                    sfc_imgs{s} = [subject_Dir, seeds{s}, filesep, subjectNames{i}, '_SFC_', num2str(num_steps),'_', seeds{s},'.nii'];
                end
                
                % USE COREGISTER (RESLICE) to coregister sfc maps to atlas
                % size
                clear jobs;
                % bring individual seven Map to size of atlas...
                jobs{1}.spm.spatial.coreg.write.ref = {fullfile(subject_Dir,['sub-',subNum,'_label-L2018_desc-',parc{1,p},'_atlas.nii'])};
                jobs{1}.spm.spatial.coreg.write.source = sfc_imgs;
                jobs{1}.spm.spatial.coreg.write.roptions.interp = 4;
                jobs{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
                jobs{1}.spm.spatial.coreg.write.roptions.mask = 0;
                jobs{1}.spm.spatial.coreg.write.roptions.prefix = 'inAtl_';
                % ALREADY GETS FLIPPED, YEY
                
                spm_jobman('run', jobs); %RUN
                clear jobs;
                
                if ~exist(fullfile(subject_Dir, '_SFC_resized_indi'),'dir')
                    mkdir(fullfile(subject_Dir, '_SFC_resized_indi'));
                end
                
                for s = 1:numel(seeds) % for each seed
                    disp(seeds{s});
                    % move the resized sfc_map
                    sfc_resized = [subject_Dir, seeds{s}, filesep,'inAtl_' ,subjectNames{i}, '_SFC_', num2str(num_steps),'_', seeds{s},'.nii'];
                    movefile(sfc_resized, fullfile(subject_Dir, '_SFC_resized_indi'));
                    
                    % load the map
                    sfc_map = load_untouch_nii([subject_Dir, '_SFC_resized_indi', filesep,'inAtl_' ,subjectNames{i}, '_SFC_', num2str(num_steps),'_', seeds{s},'.nii']);
                    
                    for r = 1:size(translator,1) % for each area
                        % disp(translator.region{r});
                        if s == 1
                            output{r,1} =  subjectNames{i};
                            output{r,2} =  parc{2,p};
                            output{r,3} =  translator.region{r};
                        end
                        % get the sfc values of that area
                        sfc_values = sfc_map.img(atlas.img == translator.msbpROI(r));
                        output{r,s+3} = mean(sfc_values, 'all', 'omitnan');
                        
                    end
                end
                
                % save out put for that subject
                outTab = cell2table(output, 'VariableNames', colNames);
                outName = [subject_Dir, subjectNames{i}, '_indiSFC_', num2str(num_steps),'_', parc{2,p},'.csv'];
                writetable(outTab, outName);
                
                disp(['DONE ', num2str(i),'/', num2str(numel(subjectNames))]);
            end % repeat for subjects            
        end % repeat for templates
    end

end
