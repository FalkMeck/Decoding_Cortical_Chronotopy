%% EXTRACT THE INDIVIDUAL ISC VALUES FOR EACH SUBJECT, CONDITION, and ROI

% define Variables
spmPath = '...\spm12';
addpath(genpath(spmPath));
NIfTItoolPath = '...\NIfTI_20140122';
addpath(NIfTItoolPath);

% all Subjects
subject_dir = {'HFG_121'};

% MSBP is in BIDS, so we need sub-?? names
% sub2num is "Translator"
sub2num = readtable('...\DV02_ISC\Subject2Num.csv');

% where are all Files, and were are they supposed to go
study_dir ='...\ISC_Analysis\'; %set this to where all your 4D files are

% all Conditions, for which ISC was calculated
conditions = {'Single', 'Triplet', 'Nonet', 'Complete'}; % create all paths

% not extremely nesseccary, but there is the translator between lausanne
% labels and MSBP numbers
msbp_atlas_dir = '...\02_UM\02_Ji_atlas\msbp_anat\';
TableNames = {'Subject','Condition','Region','VoxelNum','ISC'};

% How many subjects; 
n = numel(subject_dir);

% parcellations and scales
parc = {'scale1', 'scale2', 'scale3', 'scale4';...
    'aparc', 'lausanne120', 'lausanne250', 'lausanne500'};

outDir = '...\ISC_analysis\ROI_Analysis\'; 

for p = 3 %1:size(parc,2)
    translator = readtable([msbp_atlas_dir, 'transS',parc{1,p}(2:end), '.csv']); % provided in the ./02_UM/02_Ji_atlas-Folder
    for i = 1:n
        disp(subject_dir{i}); 
        subDir = fullfile(study_dir, subject_dir{i});
        subNum = sub2num.Number(contains(sub2num.Subject,subject_dir{i}));
        subNum = sprintf('%02d',subNum);

	% Multiscale brain parcellator parcelled file od scale 3 (lausanne250) 
	% provided in "...\minimally_Processed\task_fMRI", that has to be loaded here        
        atlasFile = [subDir, filesep, 'sub-',subNum,'_label-L2018_desc-',parc{1,p},'_atlas.nii'];
        

%         extract all areas as ROI-map
        outSub = fullfile(outDir, subject_dir{i});
        outSubROI = fullfile(outSub, ['ROIs_', parc{2,p}]);
        if ~exist(outSub,'dir') || ~exist(outSubROI,'dir')
            mkdir(outSub); mkdir(outSubROI); 
        end
        
        parcIm = load_untouch_nii(atlasFile); % Load parcelled MSBP file for this subject
        for r = 1:size(translator,1) % every cortical region of interest
            disp(['ROI: ', translator.region{r}]); % display region
            outputImFn = [outSubROI, filesep, subject_dir{i}, '_', translator.region{r}, '.nii']; % define output file
            currROI = translator.msbpROI(r); % what ROI
            idx_ROI = parcIm.img == currROI; % find ROI
            outputIm.img = zeros(size(parcIm.img));
            outputIm.img(idx_ROI) = 1; % create binary file
            outputIm.hdr = parcIm.hdr; % copy HDR
            %DIALTAE THE ROI % REquires WFU pickatlas toolbox in SPM12
            s = 1; % how many voxels?
            threed = 1; % 3D (1 = YES)
            [xd,yd,zd] = size(outputIm.img); % get sizes
            dcubes = zeros(xd, yd, zd); % create empty
            [x_sl,y_sl,z_sl] = ind2sub([xd yd zd], find(outputIm.img)); % find coordiantes of binary file
            
            for ip = 1:length(x_sl) % for each coordinate
                ix   = x_sl(ip);	iy  = y_sl(ip);	iz = z_sl(ip);
                mnx  = max( 1, ix-s);	mxx = min( ix+s, xd); % add and substract 1, except when at the edge of the image
                mny  = max( 1, iy-s);	mxy = min( iy+s, yd);
                if threed==1 % if §D was selected also do this for the Z coordiante
                    mnz = max( 1, iz-s);	mxz = min( iz+s, zd);
                else
                    mnz = iz;		mxz = iz;
                end
                dcubes(mnx:mxx, mny:mxy, mnz:mxz) = 1; % create larger binary
            end % End Coordinates
            outputIm.img = dcubes; % write dialted binary into file
            save_nii(outputIm,outputImFn); % save file
        end % End ROI
        
        % coregister all to ISC space
        standardImg = [study_dir, subject_dir{i},filesep,'ISC_LOO_maps',filesep,'ISC_LOO_corr_',conditions{1},'.nii'];
        
        % get all ROIs of this subject
        diaCont = dir([outSubROI,filesep,subject_dir{i},'*']);
        diaNames = {diaCont.name};
        diaPaths = strcat(outSubROI,filesep,diaNames);
        
        % coreg ROIs to ISC space
        clear jobs;
        jobs{1}.spm.spatial.coreg.write.ref = {standardImg}; % reference image
        jobs{1}.spm.spatial.coreg.write.source = diaPaths'; % all image to be cogesiterd
        jobs{1}.spm.spatial.coreg.write.roptions.interp = 4; % standard interpolation
        jobs{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0]; % standard warping
        jobs{1}.spm.spatial.coreg.write.roptions.mask = 0; % standard mask
        jobs{1}.spm.spatial.coreg.write.roptions.prefix = 'r_' ;
        % RUN
        spm_jobman('run', jobs);
        
        clear jobs
        
        % Get all coregistered ROIs
        diaCont = dir([outSubROI,filesep,'r_*']);
        diaNames = {diaCont.name};
        diaPaths = strcat(outSubROI, filesep, diaNames);
        
        % Binarize again, since coregistration makes files non-binary
        for ii = 1:numel(diaPaths)
            jobs{1}.spm.util.imcalc.input = diaPaths(ii);
            jobs{1}.spm.util.imcalc.output = ['bin_',diaNames{ii}(1:(end-4))];
            jobs{1}.spm.util.imcalc.outdir = {outSubROI};
            jobs{1}.spm.util.imcalc.expression = 'i1 > 0.2';
            jobs{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            jobs{1}.spm.util.imcalc.options.dmtx = 0;
            jobs{1}.spm.util.imcalc.options.mask = 0;
            jobs{1}.spm.util.imcalc.options.interp = 1;
            jobs{1}.spm.util.imcalc.options.dtype = 4;
            
            spm_jobman('run',jobs)
            clear jobs
        end % end ROIs
        
        % For each Condtion, extract ISC
        for c = 1:numel(conditions)
            % Empty output table
            Textract = cell2table(cell(size(translator,1),numel(TableNames)), 'VariableNames', TableNames);
            % For which file
            % I used the Fisher Z transformed ISC
            conFile = [study_dir, subject_dir{i},filesep,'ISC_LOO_maps',filesep,'ISC_LOO_FishZ_',conditions{c},'.nii'];
            
            for r = 1:size(translator,1) % For each ROI in translator
                disp(translator.region{r});
                
                % get the file
                roi = [outSubROI, filesep,'bin_r_' ,subject_dir{i}, '_', translator.region{r}, '.nii'];
                
                % load the file
                Y = spm_read_vols(spm_vol(roi),1);
                indx = find(Y > 0); % get coordiantes of included voxels
                [x,y,z] = ind2sub(size(Y),indx);
                XYZ = [x y z]';
                
                % extract data
                isc  = spm_get_data(conFile, XYZ); % isc
                iscClean = isc(~isnan(isc)); % remove NANs if present
                nVox = size(iscClean,2); % how many voxels
                
                cond = conditions(c);
                currRoiName = translator.region(r);
                
                % Write into Output-Table
                Textract{r,:} =[subject_dir(i), cond, currRoiName, num2cell(nVox), num2cell(mean(iscClean))];    
            end
            % Save Output-Table for this condition
            writetable(Textract,[subDir,filesep,subject_dir{i},'_ISC_',conditions{c},'.csv']);
        end
    end
end