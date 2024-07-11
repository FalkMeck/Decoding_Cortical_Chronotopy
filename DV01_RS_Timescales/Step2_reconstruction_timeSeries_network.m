function [] = reconstruction_timeSeries_network()
% RECONSTRUCTION_FUNCTIONAL_NETWORK Reconstruct functional connectivity.
%
%   reconstruction_functional_network(CONFIGPARAMS) reconstruct the
%   functional connectivity matrices for each template and each method
%   according to parameters specified in CONFIGPARAMS.

addpath(genpath('.../CATO_3.1.2/CATO-main/src')); % part of CATO installation

analysis_dir = '.../Analysis';

cd(analysis_dir);

subject_dir = {'HFG_121'};

% translate the necessary parts of config file to matlab
configParams.general.parameterPropertiesFile = 'functionalParameterProperties.xlsx'; % possibly add entire path
configParams.general.templates = {'aparc','lausanne120','lausanne250','lausanne500','economo','BB50human'}; % all possible templates
configParams.general.templatesDir = '.../CATO_3.1.2/templates';
configParams.general.ROIsFile = [configParams.general.templatesDir,'/TEMPLATE/ROIs_TEMPLATE.txt'];

repetitionTimeMsec = 1000; % TR = 1000ms

%% Initialization

for i = 1:numel(subject_dir)
    disp(['Processing: ',subject_dir{i}]);

    % Defining all necessary configuration Parameters

    configParams.general.subjectDir = [analysis_dir,'/',subject_dir{i},'/'];
    configParams.general.outputDir = [configParams.general.subjectDir,'fMRI_processed/'];
    configParams.functional_preprocessing.segmentationFile = [configParams.general.outputDir,subject_dir{i},'_aseg_ref.nii.gz'];

    configParams.parcellation.parcellationFile = [configParams.general.outputDir,subject_dir{i},'_TEMPLATE+aseg_ref.nii.gz'];
    configParams.parcellation.lutFile = [configParams.general.templatesDir,'/TEMPLATE/TEMPLATE.annot.ctab'];

    configParams.regressionMethod = '24pXaCompCorXVolterra_spikeReg';

    configParams.reconstruction_functional_network.methodDescription = "custom_0.01-0.1";
    configParams.reconstruction_functional_network.connectivityMatrixFile = [configParams.general.outputDir,subject_dir{i},'_connectivity_METHOD_TEMPLATE.mat'];
    configParams.reconstruction_functional_network.timeSeriesFile = [configParams.general.outputDir,subject_dir{i},'_time_series_METHOD_TEMPLATE.mat'];
    configParams.reconstruction_functional_network.minRepetitionTime = 100;
    configParams.reconstruction_functional_network.bandpass_filter.filter = true;
    configParams.reconstruction_functional_network.bandpass_filter.frequencies = [0.01,0.1];
    configParams.reconstruction_functional_network.saveTimeSeries = true;

    methods = fieldnames(configParams);
    methods = methods(contains(methods, 'reconstruction_functional_network'));

    pipelineV25Flag = false; % option for compatibility with old pipeline version.

    for iMethod = 1:length(methods)

        %% Initializate parameters

        % General parameters
        thisMethodDescription = configParams.(methods{iMethod}).methodDescription;
        % processed file after Python nuisance regression
        fmriProcessedFile = [analysis_dir,'/', subject_dir{i}, '/', subject_dir{i},'_fmri_', configParams.regressionMethod,'.nii'];

        segmentationFile = configParams.functional_preprocessing.segmentationFile;
        ROIsFile = configParams.general.ROIsFile;

        parcellationFile = configParams.parcellation.parcellationFile;
        templates = configParams.general.templates;
        lutFile = configParams.parcellation.lutFile;

        % Regression parameters
        %DELETED AS REGRESSION IN PYTHON
        
        % Bandpass filter parameters
        minRepetitionTime = configParams.(methods{iMethod}).minRepetitionTime;
        bandpass_filter = configParams.(methods{iMethod}).bandpass_filter;

        % Scrubbing parameters
        %DELETED AS SCRUBBING IS NOT POSSIBLE FOR CONTINOUS TIME SERIES
        %ALSO WE REGRESS OUT SPIKES
    
        % Time series parameters
        saveTimeSeriesFlag = configParams.(methods{iMethod}).saveTimeSeries;
        timeSeriesFile = configParams.(methods{iMethod}).timeSeriesFile;

        % Connectivity matrix parameters 
        % not really needed at the moment, but nice to have
        connectivityMatrixFile = configParams.(methods{iMethod}).connectivityMatrixFile;

        fprintf('method description: %s\n', thisMethodDescription);

        %% Prepare rs-fMRI data

        % Select voxels for processing that aree within the brain and have
        % signal in more than 90% of the timepoints
        segmentation = load_nifti(segmentationFile);
        segmentation = segmentation.vol;
        brainMask = segmentation(:) > 0;
        %frmi_full = load_nifti(fmriProcessedFile);   
        fmri = load_nifti_partially(fmriProcessedFile, brainMask);
        signalIntensities = single(fmri.partialvol);

        prevalenceMask = mean(signalIntensities ~= 0, 2) >= 0.9;

        if ~pipelineV25Flag % Maybe delete all sections 
            signalIntensities = signalIntensities(prevalenceMask, :);
            brainPrevalenceMask = false(size(brainMask));
            brainPrevalenceMask(brainMask) = prevalenceMask;
        else
            brainPrevalenceMask = brainMask;
        end

        %% Linear regression
        %DELETE, WAS DONE MORE THOROUGHLY IN PYTHON (like Ito et al., 2020)

        %% Bandpass filter

        % Get repetition time.
        assert(repetitionTimeMsec >= minRepetitionTime, ...
            ['Repetition time (%g msec) reported in fmriProcessedFile', ...
            ' is smaller than the minRepetitionTime (%g msec). ', ...
            'Transform repetition time to milliseconds ', ...
            'or adjust minRepetitionTime-parameter'], repetitionTimeMsec, minRepetitionTime);
        repetitionTimeSec = repetitionTimeMsec/1000;

        [filter_b, filter_a] = butter(2, 2*repetitionTimeSec*bandpass_filter.frequencies);

        % Use for-loop to avoid memory issues from having double. (1sec difference)
        % filteredSignal = filtfilt(filter_b, filter_a, double(signalIntensities(selectedVoxels, :)));
        filteredSignal = zeros(size(signalIntensities), 'single');
        for j = 1:size(signalIntensities, 1)
            filteredSignal(j,:) = filtfilt(filter_b, filter_a, ...
                double(signalIntensities(j, :)));
        end

        selectedTimeSeries = filteredSignal;

        clear filteredSignal

        %% Scrub
        % DELETE BECAUSE SCRUBBING IS COUTERPRODUCTIVE FOR TIME SERIES 
        
        %% Calculate average Time Series and the Correlation Matrix

        for iTemplate = 1:length(templates)
            thisTemplate = templates{iTemplate};
            fprintf('template: %s\n', thisTemplate);

            thisParcellationFile = strrep(parcellationFile, ...
                'TEMPLATE', thisTemplate);
            thisROIsFile = strrep(ROIsFile, ...
                'TEMPLATE', thisTemplate);
            thisConnectivityMatrixFile = strrep(strrep(connectivityMatrixFile, ...
                'METHOD', thisMethodDescription), ...
                'TEMPLATE', thisTemplate);
            thisTimeSeriesFile = strrep(strrep(timeSeriesFile, ...
                'METHOD', thisMethodDescription), ...
                'TEMPLATE', thisTemplate);
            thisLutFile = strrep(lutFile, ...
                'TEMPLATE', thisTemplate);

            parcellation = load_nifti(thisParcellationFile);
            parcellation = single(parcellation.vol(brainPrevalenceMask));
            ROIs = readmatrix(thisROIsFile);

            % get associated regionDescriptions
            LUT = readtable(thisLutFile,  'filetype', 'text', ...
                'ReadVariableNames', false);
            LUT.Properties.VariableNames = {'ROIs', 'regionDescriptions', ...
                'Color1', 'Color2', 'Color3', 'Other'};
            [~, indxROIs] = ismember(ROIs, LUT.ROIs);
            regionDescriptions = deblank(LUT.regionDescriptions(indxROIs));

            % Calculate time series averaged over all voxels of a region
            averageTimeSeries = zeros(length(ROIs), size(selectedTimeSeries, 2));
            for iR = 1:length(ROIs)
                averageTimeSeries(iR, :) = mean(selectedTimeSeries(...
                    parcellation == ROIs(iR), :), 1);
            end

            if saveTimeSeriesFlag
                save(thisTimeSeriesFile, ...
                    'averageTimeSeries', 'ROIs', 'regionDescriptions');
            end

            % Calculate correlation data
            [connectivity, pValues] = corr(averageTimeSeries');
            connectivity = connectivity .* ~eye(length(connectivity));
            pValues = pValues .* ~eye(length(connectivity));

            save(thisConnectivityMatrixFile, ...
                'ROIs','regionDescriptions', 'connectivity', 'pValues');

        end
    end % Maybe the iMethod-Loop is not really necessary
    disp(['Done ',char(string(i)),'/',char(string(numel(subject_dir)))])
end

%% End
fprintf('---reconstruction_functional_network finished----\n');

end
