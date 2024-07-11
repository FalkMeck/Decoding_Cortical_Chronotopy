%% Main script used to compute spectral content (intrinsic timescales)
clear all; clc; 
analysisDir = '.../Analysis/';
cd(analysisDir);
%
addpath('..../Raut_lag_code_master');
%
%% Setup
% Set parameters
subject_dir = ['HFG_121'];

numSubs = size(subject_dir,1);
method = 'custom_0.01-0.1_';
templates = {'aparc','lausanne120','lausanne250','lausanne500','economo','BB50human';
    82,114,219,448,87,78}; %
% economo is very wierd in CATO, includes left lateral ventricle as area 86
% is probably more likely, but also all networks are wrong

lags = -12:12;    % range of TR shifts; should be sufficient to allow for all autocovariance functions (ACFs) to decay below .5
tr = 1;
% MAYBE CHANGE: shorter TE have to check if all become < 0.5

motion_thresh = .25; % Same as ITO % .2;    % important: must match motion criteria used during preproc        
min_block_durn = (max(lags)+1)*tr;   % min. block duration (in seconds)

outdir = fullfile(analysisDir,'TimeScale_Raut');    % set directory for saving out images here
if ~ exist(outdir,'dir')
    mkdir(outdir);
end

maxCorr= zeros(numSubs, size(templates,2)); 

%% Loop over templates
for t = 1:size(templates,2)
    curParc = templates{1,t}; disp(['For Parcellation: ', curParc]); 
    num_nodes = templates{2,t};    % number of time series
     % initialize group matrices
    grp_acfs = single(nan(num_nodes,numel(lags),numSubs)); % keep all ACFs for all subjects for stats
    %% Loop over subjects
    for iSub = 1:numSubs
           
        tic
        subj = subject_dir(iSub,:);
        disp(['Processing ' subj]);
        
        subDir = fullfile(analysisDir,subj,'fMRI_processed'); 
        tsFile = [subDir,'/',subj, '_time_series_',method, templates{1,t},'.mat'];
        
        % initialize subject matrices
        subj_lags = single(nan(num_nodes)); % peak lags
        subj_ZL = subj_lags;   % zero-lag correlation
        subj_peak = subj_lags; % peak correlation (correlation at optimal lag)
          
        TS = importdata(tsFile); % read in time series matrix
        BOLD = single(TS.averageTimeSeries)';
        good = true(num_nodes,1); % use all nodes
        
%         % read in temporal mask/motion time series (e.g., FD or DVARS)
%         format = dlmread(fullfile(analysisDir,subj,'Movement_RelativeRMS.txt')) <= motion_thresh;
%         format = format(6:end); 
%               
%         FORMAT = create_blocks(format,min_block_durn,tr);
        
        %% Construct ACF with continous blocks
%         ACFs = single(zeros(num_nodes,numel(lags)));
%         nblocks = numel(FORMAT);
%         nframes = 0;
%         
        % De-mean time series
        run_mean = nanmean(BOLD,1);
        BOLD = bsxfun(@minus,BOLD,run_mean);
%         
%         % Loop over blocks of contiguous frames
%         for j = 1:numel(FORMAT)
%             nframes = nframes + numel(FORMAT{j});
%             FHCR = false(1,numel(format));
%             FHCR(FORMAT{j}) = true;
%             for i = 1:sum(good)
%                 ACFs(i,:) = ACFs(i,:) + squeeze(lagged_cov(BOLD(FHCR,i),BOLD(FHCR,i),max(lags)))';
%             end
%         end
%         
%         % Normalize ACFs based on entire run
%         for k = 1:numel(lags)
%             ACFs(:,k) = ACFs(:,k)/(nframes - abs(lags(k))*nblocks);
%         end
%         ACFs = bsxfun(@rdivide,ACFs,ACFs(:,lags==0));
%         
%         maxCorr(iSub,t) = max([ACFs(:,1);ACFs(:,numel(lags))]); 
            
        % Store
        %grp_acfs(good,:,iSub) = ACFs;
        
        %save([outdir,'/',subj, '_acf_',method, templates{1,t},'.mat'], 'ACFs');
        
     %   toc
        
        %% Construct ACF using the entire signal, since it is spike corrected
        ACF2s = single(zeros(num_nodes,numel(lags)));
               
        for i = 1:sum(good)
        	ACF2s(i,:) = ACF2s(i,:) + squeeze(lagged_cov(BOLD(1:end,i),BOLD(1:end,i),max(lags)))';
        end
        
        nframesAll = size(BOLD, 1); 
        
        % Normalize ACFs based on entire run
        for k = 1:numel(lags)
            ACF2s(:,k) = ACF2s(:,k)/(nframesAll - abs(lags(k)));
        end
        ACF2s = bsxfun(@rdivide,ACF2s,ACF2s(:,lags==0));
        
%         maxCorr(iSub,t) = max([ACF2s(:,1);ACF2s(:,numel(lags))]); 
            
        % Store
        grp_acfs(good,:,iSub) = ACF2s;
        
        save([outdir,'/',subj, '_acfull_',method, templates{1,t},'.mat'], 'ACF2s');
        
        toc
    end

%     acf_mean = nanmean(grp_acfs,3);
%     acf_mean = bsxfun(@rdivide,acf_mean,acf_mean(:,lags==0));
    
    acf2_mean = nanmean(grp_acfs,3);
    acf2_mean = bsxfun(@rdivide,acf2_mean,acf2_mean(:,lags==0));
    
    %% fit HWHM
    % group-wise
%     hwhm = acf_hwhm(acf_mean',tr); 
%     save('hwhm', [outdir,'/','GroupWise_hwhm_',method, templates{1,t},'.mat']);
%     
%     % subject-wise
%     hwhms = zeros(num_nodes,numSubs);
%     for iSub = 1:numSubs
%         tic
%         hwhms(:,iSub) = acf_hwhm(grp_acfs(:,:,iSub)',tr);
%         toc
%     end
%     save('hwhms', [outdir,'/','SubjectWise_hwhm_',method, templates{1,t},'.mat']);
    
        % group-wise
    hwhm2 = acf_hwhm(acf2_mean',tr); 
    save([outdir,'/','GroupWise_hwhm2_',method, templates{1,t},'.mat'],'hwhm2');
    
    % subject-wise
    hwhm2s = zeros(num_nodes,numSubs);
    for iSub = 1:numSubs
        tic
        hwhm2s(:,iSub) = acf_hwhm(grp_acfs(:,:,iSub)',tr);
        toc
    end
    save([outdir,'/','SubjectWise_hwhm2_',method, templates{1,t},'.mat'],'hwhm2s');
    dlmwrite([outdir,'/','SubjectWise_hwhm2_',method, templates{1,t},'.txt'],hwhm2s,'delimiter',',');
    
    
end
