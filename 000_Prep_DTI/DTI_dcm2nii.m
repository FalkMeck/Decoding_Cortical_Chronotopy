% DCM2NII
function [] = dcm2nii_prepro()

spm_dir = '...\spm12\'; % complete spm12 folder with backslash at the end
%-a Y Anonymize = YES
%-s N SPM2/Analyse = NO
%-d Y Date in filename = YES
%-e Y events in filename = YES
%-p Y protocol in filname = YES
%-i N ID in filname = NO
%-f N Sorce in filnemae = NO
%-v Y Convert every image in directory = YES
%-n Y output = .nii file = YES
%-4 N ???
%-g N gzip = NO

% DTI
% ├── DTI-1
% │   ├── 20210707_124804FORep2diso25diffmddw30RIB1000p256sls007a001.bval
% │   ├── 20210707_124804FORep2diso25diffmddw30RIB1000p256sls007a001.bvec
% │   └── 20210707_124804FORep2diso25diffmddw30RIB1000p256sls007a001.nii
% ├── DTI-2
% │   ├── 20210707_124804FORep2diso25diffmddw30RIB1000p256sls005a001.bval
% │   ├── 20210707_124804FORep2diso25diffmddw30RIB1000p256sls005a001.bvec
% │   └── 20210707_124804FORep2diso25diffmddw30RIB1000p256sls005a001.nii
% ├── DTIB0
% │   ├── 20210707_124804FORep2diso25diffmddw30RIB0p256sls009a1001.nii
% │   ├── 20210707_124804FORep2diso25diffmddw30RIB0p256sls010a1001.nii
% │   └── 20210707_124804FORep2diso25diffmddw30RIB0p256sls011a1001.nii
% ├── GreField-1
% │   ├── 20210707_124804FORgrefieldmappingsl25mms003a1001_1.nii
% │   └── 20210707_124804FORgrefieldmappingsl25mms003a1001_2.nii
% ├── GreField-2
% │   └── 20210707_124804FORgrefieldmappingsl25mms004a2001.nii
% └── ReversePhaseEncoding
%     └── 20210707_124804ReversePhaseEncodingTMPs012a1001.nii


study_dir = '...\DTI\Analysis';           % definiere Ordnerpfad zu Ordner mit DTI analysis
dicom_dir = '...\RICHIE\Day2';            % folder path data day 2

subject_dir = ['HFG_121'];

number_subjects = size(subject_dir,1); 

subFolder_names = ["DWI", "DTI_1","DTI_2","DTIB0","GreField_1","GreField_2",...
    "ReversePhaseEncoding"]; % names the Subfolders in DTI will get
folder_id = ["MT", "007","005","B0","mm_e","mm_phase","ReversePhaseEncodingTMP"]; 
% Items to search for in each folder to find correct folders/files

tic
     make_folders()
     dicom2nifti_dti() % does not uses SPM, since SPM doesn't handle the files as easily as dicm2nii
     reorientation() % this works with SPM
     rename_files()
toc

%% Create folders

function [] = make_folders()
   % Anatomy
    if ~exist(fullfile(study_dir,'1_Anatomy'),'dir') % Preparation for freesurfer
        success = mkdir(fullfile(study_dir,'1_Anatomy'));
        if ~success, error('Cannot create new subject directory.\n'); end
    end
    for i= 1:number_subjects
        if ~exist(fullfile(study_dir,subject_dir(i,:)),'dir')
            fprintf('Creating new subject directory %s.\n',subject_dir(i,:));
            success = mkdir(fullfile(study_dir,subject_dir(i,:)));
            if ~success, error('Cannot create new subject directory.\n'); end
        end
        if ~exist(fullfile(study_dir,subject_dir(i,:),'DTI'),'dir')
            mkdir(fullfile(study_dir,subject_dir(i,:),'DTI'));
        end
        
        for sfn = 1:size(subFolder_names,2)
            if ~exist(fullfile(study_dir,subject_dir(i,:),'DTI',subFolder_names(1,sfn)),'dir')
                mkdir(fullfile(study_dir,subject_dir(i,:),'DTI',subFolder_names(1,sfn)));
            end
        end
        fprintf('Creating Directories:\n Done.\n')
    end
end

%% DICOM to Nifti (with dicm2nii)
function [] = dicom2nifti_dti()
addpath(genpath('...\dcm2nii'));

    for i = 1:number_subjects

        disp('Dicom Conversion for: ');
        disp(subject_dir(i,:));

        % specification of folder structure to navigate to teh raw data
        allRawFiles = dir(dicom_dir);
        subFolders = {allRawFiles.name}; % all folders in the folder which contains raw data folder
        subRaw=char(subFolders(startsWith(subFolders,subject_dir(i,:))));
        % find the raw_data folder in these folders
        filesSubRaw=dir([dicom_dir,'\' ,subRaw]); % all folders in raw data folder
        filesSubRaw={filesSubRaw.name}; % names of all folders in raw data folder
        expFolder=char(filesSubRaw(startsWith(filesSubRaw,'BIPSY_RICHI_SMS'))); 
        dicm_dir = fullfile(dicom_dir, subRaw,expFolder,'/'); % path to all dicom files 
        nii_dir = fullfile(study_dir, subject_dir(i,:),'DTI','DWI'); % path to 

        dicm2nii(dicm_dir,nii_dir,'.nii');
        % goeas through all folders following dicm_dir and converting all
        % convertable files to nifti as either one imate or a time series

    end
    fprintf('Dicom Conversion:\n Done.\n')
end

%% Reorientation

% Only works as script if the data has already been reoriented manually and
% reoritation matrices from SPM are available

    function [] = reorientation()
        
        for i = 1:number_subjects
                
                disp('Reorientation for: ')
                disp(subject_dir(i,:))
                
                % path to DWI folder
                dti_dir = fullfile(study_dir, subject_dir(i,:), 'DTI','DWI\'); 
                
                % select data
                [func_filenames, ~] = spm_select('List',dti_dir, '.*\.nii$');
                files_for_reorient = strcat(dti_dir,func_filenames);
                
               
                % If reorientation has already been carried out, a reorientation matrix should exist in the Anatomy folder by default
                % Read this matrix (has an unusual structure)
                % Depending on where the old reorientation matrix is ​​located,
                % a new path may need to be defined (change highres_reo assignment)
                % (e.g. if you start from the beginning and the folder structure has changed)   highres_reo = '...\Reorientation_Matrices_Day2\';
                reorient_matData = dir(fullfile(highres_reo, [subject_dir(i,:),'*']));
                
                reorient_matrix =load(fullfile(reorient_matData.folder,reorient_matData.name));
                                
                jobs{1}.spm.util.reorient.srcfiles = cellstr(files_for_reorient);
                jobs{1}.spm.util.reorient.transform.transM = reorient_matrix.M;
                jobs{1}.spm.util.reorient.prefix = '';
                
                spm_jobman('run', jobs);
                clear jobs
        end
        fprintf('Reorientation:\n Done.\n')
    end

%% Generate Names and Folderstructure

 function [] = rename_files()
 
     for i = 1:number_subjects
        
        disp('Doing folder organisation for: ');
        disp(subject_dir(i,:));

        
        nii_dir_b = fullfile(study_dir, subject_dir(i,:),'DTI'); % DTI directory of the participant
        head = load(fullfile(nii_dir_b,'DWI','dcmHeaders.mat'));
        study_date = [head.h.t1_mprage_sag_p2_iso_1mm.InstanceCreationDate,'_', head.h.t1_mprage_sag_p2_iso_1mm.InstanceCreationTime(1:6)];
        % Load study date and time from header inforamtion
        % create a character vector that can be used in naming files
      
        dwi_files = dir([nii_dir_b,'/DWI']); dwi_files = {dwi_files.name}; % Listing all files in this directory
        dwi_files_clean = cellfun(@(x) erase(x,"_"),dwi_files, 'UniformOutput', false);
        dwi_files_clean = cellfun(@(x) [study_date,'_',x], dwi_files_clean, 'UniformOutput', false);
        % creating renamed names with study_date vector fro all of them
        % (cellfun)
        
        for sfn = 2:size(subFolder_names,2) %Each folder
            File_i = cell(2,1);
            File_i{1,1} = dwi_files(1,contains(dwi_files, folder_id(sfn))); 
            % identifiyng all files that have to go in the folder by their
            % specified folder_id
            File_i{2,1} = dwi_files_clean(1,contains(dwi_files, folder_id(sfn)));
            % sort correct file name to these identified files
            for fi = 1:size(File_i{1,1},2)
                source_file = fullfile(nii_dir_b,'DWI',File_i{1,1}{1,fi});
                destination_file = fullfile(nii_dir_b,subFolder_names(sfn),File_i{2,1}{1,fi});
                if ~exist(destination_file) 
                    movefile(source_file, destination_file);
                end
            end
        end
        
      % Antatomy
      source_file = fullfile(nii_dir_b,'DWI','t1_mprage_sag_p2_iso_1mm.nii');
      destination_file = fullfile(study_dir, '1_Anatomy', [subject_dir(i,:),'.nii']);
      if ~exist(destination_file)
          movefile(source_file,destination_file);
      end        
     end  
     fprintf('Folder Organisation:\n Done.\n')
 end
end