%% SCRIPT HAS TO BE RUN FROM A UNIX-PC %%
function [] = dti_unix_prepro()
% DTI Data preprocessing steps with FSL and CATO toolbox
cd('...');

addpath(genpath('.../freesurfer/matlab'));
addpath(genpath('.../CATO_3.1.2/CATO-main/src'));

study_dir = [pwd,'/DTI/Analysis'];
unix_dir = '.../';


subject_dir = ['HFG_121'];

% all Subjects
number_subjects = size(subject_dir,1);

%% Run function
tic
   FSL_merge_DTI()
   move_revIm()
   move_fs_files()
   fs_recon_all()
   CATO_struct_pipe()
toc

%% Merging DWI/DTI files

    function [] = FSL_merge_DTI()
        prepare_DWI_shell = [study_dir,'/prepareDWI_for2107_FSL6.sh '];
        % Description of the bash-File "prepareDWI_for2107_FSL6.sh"
        %     #!/bin/bash
        %     # author: Dominik Grotegerd
        %     # institution: Translational Psychiarty / University of Muenster
        %     # Version 1.1 / 2018
        %
        %     subDir=$(readlink -e $1) # uses first character string pasted after running shell from command as input for subDir
        %     # necessery to loop through code via Matalb unix()-command
        %     subBase=$(basename $subDir) # basename takes the last part of subDir (here name of participant) for subBase
        %
        %     # overwrite copy flag
        %     cpFlag="-i" # unclear what this is used for since it just defines a variable and
        %     cpFlag=""
        %
        %     # FSL configuration
        %     #/usr/local/fsl/etc/fslconf/fsl.sh # maybe needs to be changed, but is just comment not code ???
        %
        %     # Naming Conventions # Defining Names and Folders for variables that get created in shell script
        %     MagImName="$subBase"_MagnitudeImage
        %     phaImName="$subBase"_PhaseImage
        %     revImName="$subBase"_ReversedImage
        %     mergedImName="$subBase"_dwi_merged.nii.gz
        %     mergedImNameRPE="$subBase"_dwi_topuped_merged.nii.gz
        %     B0_ref="$subBase"_B0reference.nii.gz
        %     outDir=$subDir/DTI/DWI_prepared
        %
        %     bvName=$(echo $mergedImName | sed 's/.nii.gz//g') # taking the mergedImName(RPE) and cuting of the ".nii.gz"-Ending
        %     bvNameRPE=$(echo $mergedImNameRPE | sed 's/.nii.gz//g') # to create new Names fro the bval and bvec variables
        %
        %
        %     # simple check if Subject Directory/Path exists in the specified form
        %     if [ ! -d "$subDir"/DTI/ ]; then
        %         echo $subDir DTI directory does not exist
        %         exit 1
        %     fi
        %
        %     mkdir -p $outDir # create the output directory $subDir/DTI/DWI_prepared
        %     cd $outDir # set current directory to that path
        %
        %
        %     # mag image / phase image # define the magnitude and pahse image from the gradient fields
        %     magIm="$subDir"/DTI/GreField_1/$(ls "$subDir"/DTI/GreField_1/ | head -1)
        %     phaIm="$subDir"/DTI/GreField_2/*.nii
        %     revIm="$subDir"/DTI/ReversePhaseEncoding/*.nii
        %
        %     # get FieldMap Images
        %     cp $cpFlag $phaIm $phaImName.nii # copy gradient niftis to new file (name) specified in Naming Conventions
        %     cp $cpFlag $magIm $MagImName.nii
        %     cp $cpFlag $revIm $revImName.nii
        %
        %     # merge simple dwi
        %     $FSLDIR/bin/fslmerge -t "$mergedImName" "$subDir"/DTI/DTI_1/*.nii "$subDir"/DTI/DTI_2/*.nii "$subDir"/DTI/DTIB0/*nii
        %     # /usr/local/fsl/bin/fslmerge -t "$mergedImName" "$subDir"/DTI/DTI_1/*.nii "$subDir"/DTI/DTI_2/*.nii "$subDir"/DTI/DTIB0/*nii
        %     # merge the multiple DWI nifitis to one file
        %
        %     # merge bval / bvecs
        %     paste   "$subDir"/DTI/DTI_1/*bvec "$subDir"/DTI/DTI_2/*bvec  <(echo 0 0 0; echo 0 0 0; echo 0 0 0)  | column -s $'\t' -t > "$bvName".bvec
        %     paste   "$subDir"/DTI/DTI_1/*bval "$subDir"/DTI/DTI_2/*bval  <(echo 0 0 0) | column -s $'\t' -t > "$bvName".bval
        %     # paste the two bvec and beval files into one with 0 0 0 in the third column
        %     # so 3 columns 3(bvec) or 1(bval) lines
        %     # 1st col = all 3/1 lines from DTI_1 files
        %     # 2nd col = all 3/1 lines from DTI_2 files
        %     # 3rd col = 0 0 0 the three B0 Images
        %
        %     # if there's Reverse Phase Encoding
        %     if [ -d "$subDir"/DTI/ReversePhaseEncoding ]; then
        %         $FSLDIR/bin/fslmerge -t "$mergedImNameRPE"  "$mergedImName" "$subDir"/DTI/ReversePhaseEncoding/*.nii 
        %         # /usr/local/fsl/bin/fslmerge -t "$mergedImNameRPE"  "$mergedImName" "$subDir"/DTI/ReversePhaseEncoding/*.nii
        %         paste "$bvName".bvec <(echo -0 -0 -0 -0; echo -0 -0 -0 -0; echo -0 -0 -0 -0)  | column -s $'\t' -t > "$bvNameRPE".bvec
        %         paste "$bvName".bval <(echo -0 -0 -0 -0 ) | column -s $'\t' -t > "$bvNameRPE".bval
        %     fi
        %     # If a ReversePhaseEncoding-Image exists then merge the merged files from DTI_1 and DTI_2 and DTIB0 with this image # and create new files in the same way which now have the filename "dwi_topuped_merged"
        %     # for pasting add column with -0 -0 -0 -0 to represent 4 reverese(-) Images
        
        for i = 1:number_subjects
            subDir = fullfile(study_dir,subject_dir(i,:));
            unix_command = [prepare_DWI_shell, subDir];
            unix(unix_command);
            disp(['DONE: ', subject_dir(i,:)]);
        end
        
    end

%% Move the ReversedPhaseImage to the DWI_prepared folder
% only necessary before changes to prepareDWI_for2107_FSL6.sh: 22.02.2022

    function [] = move_revIm()
        bipsy_analysis = '.../DTI/Analysis/';
        for i = 1:number_subjects
            revIm_folder = [bipsy_analysis, subject_dir(i,:),'/DTI/ReversePhaseEncoding/'];
            revIm_folder_names = dir(revIm_folder);
            revIm_folder_names = {revIm_folder_names.name};
            revIm_file_pos = contains(revIm_folder_names , 'ReversePhaseEncoding');
            revIm_file = revIm_folder_names{revIm_file_pos};
            revIm_source = [revIm_folder, revIm_file];
            revIm_dest =  [bipsy_analysis, subject_dir(i,:),'/DTI/DWI_prepared/',subject_dir(i,:),'_ReversedImage.nii'];
            copyfile(revIm_source, revIm_dest);
        end
    end


%% Move files for freesurfer

    function [] = move_fs_files()
        % Freesurfer creates a symlink that windows cannot handle, all files need
        % to be copied to the Linux system/ a external drive before the recon_all
        % can start
        
        bipsy_anat = '.../DTI/Analysis/1_Anatomy';
        unix_anat = '.../Freesurfer_data';
        
        copyfile(bipsy_anat, unix_anat);
               
    end
%% freesurfer recon_all

    function [] = fs_recon_all()
        
        clc;
        cd([unix_dir,'/Freesurfer_data/']);
        
        cores = feature('numcores');
        if cores > 6 % with all 8 cores in use, probably not enough RAM
            cores = 6;
        end
        nbycore = ceil(number_subjects/cores);
        all_subjects_cell = cell(nbycore,1);
        for groups = 1:nbycore
            start = (groups-1)*cores+1;
            ending = min((groups)*cores,number_subjects);
            allsubjects = [];
            for sub_i = start:ending
                allsubjects = [allsubjects,subject_dir(sub_i,:),'.* '];
            end
            all_subjects_cell{groups,1} = allsubjects;
        end
        
        for groups = 1:8
            unix_command = ['ls ', all_subjects_cell{groups,1},...
                '| parallel --jobs ', num2str(cores),...
                ' ./1_runFreesurfer.sh'];
            unix(unix_command);
        end
        
    end
%% CATO
    function [] = CATO_struct_pipe()

%         
%         structural_pipeline(subjectDirectory,'configurationFile', configFile, 'runType', 'overwrite');
%         % seems to work despite Freesurfer 7
%         % Stuctural_pipleline command was adjusted to FREESURFER v7

          analysis_dir = [unix_dir,'CATO_Analysis/'];
          bigbipsy2 = '.../';
          

          configFile = [analysis_dir, 'Structural_configuration_merged_topup_eddy_subcortical.conf']; 
          % only necessary if structurall_pipeline runs via MATLAB and not
          % via Executables/Shell-script
                   
           for i = 1:number_subjects
              % Move and rename Folders and Files from FSL and FREESURFER 
              % has to be done no matter if CATO is run by MATLAB or by shell
              % FROM: 
              source_dwi = [bigbipsy2, 'DTI/Analysis/', subject_dir(i,:),'/DTI/DWI_prepared/'];
              source_fs = [unix_dir, 'Freesurfer_data/2_Anatomy_Reoriented/Day2/', subject_dir(i,:), '_D2/'];
              
              % TO:
              % create Subject Directory
               subjectDirectory = fullfile(analysis_dir,subject_dir(i,:));
              if ~exist(subjectDirectory,'dir') 
                 success = mkdir(subjectDirectory);
                 if ~success, error('Cannot create new subject directory.\n'); end
              end
              % create DWI folder
              dwi_folder = fullfile(analysis_dir, subject_dir(i,:),'DWI');
              if ~exist(dwi_folder,'dir')
                  mkdir(dwi_folder);
              end
              % create FreeSurfer folder
              freesurfer_folder = fullfile(analysis_dir, subject_dir(i,:),'FreeSurfer');
              if ~exist(freesurfer_folder, 'dir')
                  mkdir(freesurfer_folder)
              end
                              
              % Copy files from Bigbipsy/$HOME/Documents/Freesurfer_data
              copyfile(source_dwi, dwi_folder);
              copyfile(source_fs, freesurfer_folder);
              
              % RUN CATO pipeline via MATLAB (for each participant) 
              structural_pipeline(subjectDirectory,'configurationFile', configFile, 'runType', 'overwrite');
              % Stuctural_pipleline command was adjusted to FREESURFER v7 
           end
          
    end        
end
