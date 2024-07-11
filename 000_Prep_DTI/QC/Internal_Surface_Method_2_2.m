%% Step 2: Quality Checking Cortical Measures (FreeSurfer)
%% 2.2. Internal Surface Method

% Within Matlab
addpath('.../ENIGMA_Cortical_QC_2.0');
FS_Path = '.../Freesurfer_data/2_Anatomy_Reoriented/Day2/';
cd(FS_Path)

subjectNames= {'HFG_121'};

QC_output_dir = [FS_Path, 'QC']; 

suffix = '_D2';

mkdir(QC_output_dir);

for x = 1:size(subjectNames,1)
    b = subjectNames{x};
    try
        func_make_corticalpngs_ENIGMA_QC(QC_output_dir, b, [FS_Path, b, suffix, '/mri/orig.mgz'], [FS_Path, b, suffix, '/mri/aparc+aseg.mgz']);
    end
    display(['Done with subject: ', b, ': ', num2str(x), ' of ', num2str(size(subjectNames,1))]);
end



