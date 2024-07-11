%% HUB ANALYSIS RICHIE II - Single Subject

clear all;
clc;
addpath('.../BCT'); % REQUIRMENTS: BrainConnectivityToolbox

subjectNames= {'HFG_121'};

cato_folder = '.../minimally_Processed/DTI/';
cd(cato_folder);

analysis_folder = '.../';

if ~exist(fullfile(analysis_folder,'2_RCAnalysis'),'dir')
    mkdir(fullfile(analysis_folder,'2_RCAnalysis'));
end

recon_met = 'gqi_dti';
parc = {'lausanne250', 219};

numSubs = numel(subjectNames);

% Random parameters
nRandom = 2500;
rewire = 10;

group_prev_thresh = 0.6;

connectVars = {'Subj', 'connectMat', 'RC', 'randomConnectMat', 'randomRC', 'normRC'};

for p = 1:size(parc,1)
    nx = parc{p,2};
    groupMatrix = zeros(nx,nx);
    for i = 1:numSubs
        matrix_path = [cato_folder, subjectNames{i}, ...
            '/DWI_processed/', subjectNames{i}, '_connectivity_', recon_met, '_',...
            parc{p,1}, '.mat'];
        connect_mat = load(matrix_path); % load the connectivity matix
        if i == 1
            parcInfo.(parc{p,1}).weightDescriptions = connect_mat.weightDescriptions; % save names of weights onces
            parcInfo.(parc{p,1}).regionDescriptions = connect_mat.regionDescriptions; % save names of ROIs once
            parcInfo.(parc{p,1}).ROIs = connect_mat.ROIs; % save numbers of ROIs once
        end
        
        outDir = fullfile(analysis_folder, '2_RCAnalysis' ,subjectNames{i}, parc{p,1});
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        
        connect_mat = double(connect_mat.connectivity(:,:,1) >= 3);
        groupMatrix = groupMatrix + connect_mat;
    end
    % 1. Prevalence map: hard cut off
    groupMatrix = groupMatrix/numSubs;
    groupMatrix = double(groupMatrix >= group_prev_thresh);
    
    save([analysis_folder,'2_RCAnalysis/', 'groupMatrix_', parc{p,1},'.mat'], 'groupMatrix');
end
save([analysis_folder,'2_RCAnalysis/parcInfo.mat'], 'parcInfo');


%% 1. Extract binary matrices per subject, with at least 3 connections
for p = 1:size(parc,1)
    nx = parc{p,2};
    
    for i = 1:numSubs
        connect_subj_bin = cell(1,14);
        disp(subjectNames{i});
        outDir = fullfile(analysis_folder, '2_RCAnalysis' ,subjectNames{i}, parc{p,1});
        
        connect_subj_bin(1,1) = subjectNames(i);
        matrix_path = [cato_folder, subjectNames{i}, ...
            '/DWI_processed/', subjectNames{i}, '_connectivity_', recon_met, '_',...
            parc{p,1}, '.mat'];
        connect_mat = load(matrix_path); % load the connectivity matix
        
        % With NOS
        connect_subj_bin{1,2} = double(connect_mat.connectivity(:,:,1) >= 3);
        
        %% 2. Random networks per subject and calculate RC coeffient in both true and random netowrks
        Cwm = connect_subj_bin{1,2}; %get wighted connectivty matrix
        connect_subj_bin{1,3} = rich_club_bu(Cwm); % calculate weighted RC Coefficient
        
        randomConMat = zeros(nx,nx,2500);
        randomRC = cell(nRandom,1);
        for iRan = 1:nRandom % randomize network 2500 times
            randomConMat(:,:, iRan) = randmio_und(Cwm, rewire);
            randomRC{iRan,1} = rich_club_bu(randomConMat(:,:, iRan)); % calculate weighted RC Coefficient
        end
        connect_subj_bin{1,4} = randomConMat;
        connect_subj_bin{1,5} = randomRC;
        
        %% 3. RC normalization
        RC = connect_subj_bin{1,3}; % true weighted RC Coefficients
        RCRand = cell2mat(randomRC(:,1));
        RCRandMean = mean(RCRand,1); % mean random weighted RC Coefficent
        connect_subj_bin{1,6} = RC./RCRandMean; % normalized RC Coefficient
        
        % find k differences
        kComp = zeros(3,numel(RC));
        for ki = 1:numel(RC)
            kComp(1,ki) = ki;
            
            kComp_i = RC(ki);
            kComp(2,ki) = kComp_i;
            kComp_Rand = RCRand(:,ki);
            pvalue = sum(kComp_Rand >= kComp_i)/nRandom;
            
            kComp(3,ki) = pvalue;
        end
        
        connect_subj_bin{1,7} = kComp;
        
        kCompCorrected = kComp(:, ~isnan(RC));
        [kCompCorrected(5,:),~,~, kCompCorrected(4,:)] = fdr_bh(kCompCorrected(3,:),0.05,'dep','yes');
        
        connect_subj_bin{1,8} = kCompCorrected;
        
        %% 3.5 Working with the identified range of k
        % identify largest range of subsequet significant k
        % 1st case: highest and lowest to define range
        % 2nd case: intermitted by 1/ 2 non-neughboring p-values
        % accept as outlier
        % 3rd case: intermitted by more neighboring p-values
        % 1. 1.5 range of other (longest range)
        % 2. highest k in range
        % 3. exclude
        
        
        k_range = connect_subj_bin{1,8};
        kRangeSig = k_range(1, k_range(5,:)==1);
        grouped = mat2cell(kRangeSig, 1, diff([0, find(diff(kRangeSig) ~= 1), length(kRangeSig)]));
        
        if numel(grouped) > 1
            distance = [];
            for b = 2:numel(grouped)
                distance = [distance, min(grouped{b})-(max(grouped{b-1})+1)]; % how far are the clusters apart
            end
            if sum(distance) > numel(distance) % case 3: if the distnace in more than 1
                size_grouped = zeros(numel(grouped),1);
                for g = 1:numel(grouped)
                    size_grouped(g) = numel(grouped{g}); % save sized
                end
                max_g = max(size_grouped); % get max size
                rest_grouped = size_grouped(size_grouped ~= max_g); % get other sized
                % check id max is at least 1.5 times as large as the other
                if sum((max_g./rest_grouped) >= 1.5) ~= numel(rest_grouped)
                    disp(subjectNames{i}); % get max k-values
                    % case does not happen, dont waste time right now
                else
                    connect_subj_bin{1,9} =  grouped{size_grouped == max_g}; % take longest
                end
            else
                connect_subj_bin{1,9} = min(kRangeSig):max(kRangeSig);  % ignore breaks
            end
        else
            connect_subj_bin{1,9} = kRangeSig; % if just 1 group, take whole group
        end
        
        connect_subj_bin{1,10}  = length(connect_subj_bin{1,9}); % length of RC regime
        connect_subj_bin{1,11} = min(connect_subj_bin{1,9}); % min border % start of regime
        connect_subj_bin{1,12} = max(connect_subj_bin{1,9}); % max border
        connect_subj_bin{1,13} = max(connect_subj_bin{1,6}(connect_subj_bin{1,9})); % peak of norm in range
        if ~isempty(connect_subj_bin{1,13})
            connect_subj_bin{1,14} = min(find(connect_subj_bin{1,6} == connect_subj_bin{1,13},1,'first')); % peak regime position
        else
            connect_subj_bin{1,14} = [];
        end
        save([outDir, '/connect_subj_bin_', parc{p,1},'_',subjectNames{i},'.mat'], 'connect_subj_bin',  '-v7.3');
    end
end
%% 4.1 Hub Analysis Step 1
load([analysis_folder,'2_RCAnalysis/parcInfo.mat']);

varLabels = {'Region', 'ROI_CATO', 'DegreeCent',  'ClosenessCent', 'BetweenessCent',...
    'WithinModuleDegreeZ_Group', 'WithinModuleDegreeZ_Sub','ParticipationCoef_Group', 'ParticipationCoef_Sub',};
partitions = 1000;

for p = 1:size(parc,1)
    nx = parc{p,2};
    hubNum = floor(nx*.15);
    
    name = parc{p,1};
    
    ROIs = parcInfo.(name).ROIs;
    regionDescriptions = parcInfo.(name).regionDescriptions;
    
    load([analysis_folder,'2_RCAnalysis/groupMatrix_',parc{p,1},'.mat']);
    
    % module Structure
    if exist([analysis_folder,'2_RCAnalysis/repPart_group_',parc{p,1},'_prev.mat'], 'file')
        load([analysis_folder,'2_RCAnalysis/repPart_group_',parc{p,1},'_prev.mat']);
    else
        gamma = 0.5:0.1:3.5;
        [repPart_group, ~] = consensual_partition(groupMatrix, gamma, partitions);
        save([analysis_folder,'2_RCAnalysis/', 'repPart_group_', parc{p,1},'_prev.mat'], 'repPart_group');
    end
    
    for i = 1:numSubs
        outDir = fullfile(analysis_folder, '2_RCAnalysis' ,subjectNames{i}, parc{p,1});
        
        load([outDir, '/connect_subj_bin_', parc{p,1},'_',subjectNames{i},'.mat']);
        
        subMat = connect_subj_bin{1,2};
        
        %% a) top 15% degree distribution
        HubDef.(name).A.Subj(i) = subjectNames(i);
        HubDef.(name).A.degreeDist{i} = degrees_und(subMat);
        [HubDef.(name).A.degreeSort{i},...
            HubDef.(name).A.degreeSortIdx{i}] = sort(HubDef.(name).A.degreeDist{i}, 'descend');
        Cutoff = HubDef.(name).A.degreeSort{i}(hubNum);
        %         HubDef.(name).A.hubROIs{i} = ROIs(HubDef.(name).A.degreeSortIdx{i}(1:hubNum));
        %         HubDef.(name).A.hubRegions{i} = regionDescriptions(HubDef.(name).A.degreeSortIdx{i}(1:hubNum));
        HubDef.(name).A.hubROIs{i} = ROIs(HubDef.(name).A.degreeDist{i} >= Cutoff);
        HubDef.(name).A.hubRegions{i} = regionDescriptions(HubDef.(name).A.degreeDist{i} >= Cutoff);
        
        
        %% c) nodal degree >= starting point of RC regime
        
        if ~isempty(connect_subj_bin{1,9})
            startRCregime = connect_subj_bin{1,11};
            HubDef.(name).C.Subj(i) = subjectNames(i);
            HubDef.(name).C.degreeDist{i} = degrees_und(subMat);
            HubDef.(name).C.hubROIs{i} = ROIs(HubDef.(name).C.degreeDist{i} >= startRCregime);
            HubDef.(name).C.hubRegions{i}= regionDescriptions(HubDef.(name).C.degreeDist{i} >= startRCregime);
        else
            disp(subjectNames{i});
        end
        
        %% d) nodal degree >= peak of normalized RC coefficient
        
        if ~isempty(connect_subj_bin{1,9})
            peakRCregime = min(connect_subj_bin{1,14});
            HubDef.(name).D.Subj(i) = subjectNames(i);
            HubDef.(name).D.degreeDist{i} = degrees_und(subMat);
            HubDef.(name).D.hubROIs{i} = ROIs(HubDef.(name).D.degreeDist{i} >= peakRCregime);
            HubDef.(name).D.hubRegions{i} = regionDescriptions(HubDef.(name).D.degreeDist{i} >= peakRCregime);
        else
            disp(subjectNames{i});
        end
        
        %%  e) top 33% in 4 of 5 hubscores
        
        hubSc = floor(nx*(1/3));
        %meanMatNOS_thresh = groupMatrices.(chosen_met).(chosen_parc).connectGroup.NOS.(corrections{c});
        
        % module Structure
        if exist([outDir,'/repPart_',parc{p,1},'_',subjectNames{i},'.mat'], 'file')
            load([outDir,'/repPart_',parc{p,1},'_',subjectNames{i},'.mat']);
        else
            gamma = 0.5:0.1:3.5;
            [repPart_gamma, ~] = consensual_partition(subMat, gamma, partitions);
            save([outDir,'/repPart_',parc{p,1},'_',subjectNames{i},'.mat'], 'repPart_gamma');
        end
        
        HubDef.(name).E.Subj(i) = subjectNames(i);
        HubDef.(name).E.Results{i} = zeros(nx, 5);
        
        % write centrality measures in a file
        centi = zeros(numel(ROIs),8);
        centi(:,1) = ROIs;
        
        %nodal degree
        HubDef.(name).E.nodeDegree{i} = degrees_und(subMat);
        centi(:,2) = HubDef.(name).E.nodeDegree{i};
        [HubDef.(name).E.nodeDegreeSort{i},HubDef.(name).E.nodeDegreeIdx{i}] = sort(HubDef.(name).E.nodeDegree{i}, 'descend');
        Cutoff = HubDef.(name).E.nodeDegreeSort{i}(hubSc);
        %         HubDef.(name).E.Results{i}(HubDef.(name).E.nodeDegreeIdx{i}(1:hubSc),1) = ...
        %             HubDef.(name).E.Results{i}(HubDef.(name).E.nodeDegreeIdx{i}(1:hubSc),1) +1;
        HubDef.(name).E.Results{i}(HubDef.(name).E.nodeDegree{i} >= Cutoff,1) = ...
            HubDef.(name).E.Results{i}(HubDef.(name).E.nodeDegree{i} >= Cutoff,1) +1;
        
        %betweeness centrality
        HubDef.(name).E.btwCentral{i} = betweenness_bin(subMat);
        centi(:,4) =  HubDef.(name).E.btwCentral{i};
        [HubDef.(name).E.btwCentralSort{i}, HubDef.(name).E.btwCentralIdx{i}] = sort(HubDef.(name).E.btwCentral{i}, 'descend');
        Cutoff = HubDef.(name).E.btwCentralSort{i}(hubSc);
        HubDef.(name).E.Results{i}(HubDef.(name).E.btwCentral{i} >= Cutoff,2) = ...
            HubDef.(name).E.Results{i}(HubDef.(name).E.btwCentral{i} >= Cutoff,2) +1;
        
        %nodal path length
        HubDef.(name).E.nodalPL{i} = (nx-1)./sum(distance_bin(subMat),2);
        centi(:,3) = HubDef.(name).E.nodalPL{i};
        [HubDef.(name).E.nodalPLSort{i},HubDef.(name).E.nodalPLIdx{i}] = sort(HubDef.(name).E.nodalPL{i}, 'descend');
        Cutoff = HubDef.(name).E.nodalPLSort{i}(hubSc);
        HubDef.(name).E.Results{i}(HubDef.(name).E.nodalPL{i} >= Cutoff,3) = ...
            HubDef.(name).E.Results{i}(HubDef.(name).E.nodalPL{i} >= Cutoff,3) +1;
        
        %between-module participation coefficent
        HubDef.(name).E.PartCoef{i} = participation_coef(subMat,repPart_group);
        centi(:,7) =  HubDef.(name).E.PartCoef{i};
        centi(:,8) =  participation_coef(subMat,repPart_gamma);
        [HubDef.(name).E.PartCoefSort{i},HubDef.(name).E.PartCoefIdx{i}] = sort(HubDef.(name).E.PartCoef{i}, 'descend');
        Cutoff = HubDef.(name).E.PartCoefSort{i}(hubSc);
        HubDef.(name).E.Results{i}(HubDef.(name).E.PartCoef{i} >= Cutoff,4) = ...
            HubDef.(name).E.Results{i}(HubDef.(name).E.PartCoef{i} >= Cutoff,4) +1;
        
        %within-module degeree z-score
        HubDef.(name).E.ModzScore{i} = module_degree_zscore(subMat,repPart_group);
        centi(:,5) = HubDef.(name).E.ModzScore{i};
        centi(:,6) = module_degree_zscore(subMat,repPart_gamma);
        [HubDef.(name).E.ModzScoreSort{i},HubDef.(name).E.ModzScoreIdx{i}] = sort(HubDef.(name).E.ModzScore{i}, 'descend');
        Cutoff = HubDef.(name).E.ModzScoreSort{i}(hubSc);
        HubDef.(name).E.Results{i}(HubDef.(name).E.ModzScore{i} >= Cutoff,5) = ...
            HubDef.(name).E.Results{i}(HubDef.(name).E.ModzScore{i} >= Cutoff,5) +1;
        
        % Combine results
        HubDef.(name).E.Results{i}(:,6) = sum(HubDef.(name).E.Results{i}(:,1:5),2);
        
        HubDef.(name).E.ResultsAll{i} = HubDef.(name).E.Results{i}(:,6) >=4;
        
        HubDef.(name).E.hubROIs{i} = ROIs(HubDef.(name).E.ResultsAll{i});
        HubDef.(name).E.hubRegions{i} = regionDescriptions(HubDef.(name).E.ResultsAll{i});
        
        centTabI = table(regionDescriptions, centi(:,1),centi(:,2),centi(:,3),centi(:,4),centi(:,5),centi(:,6),centi(:,7),centi(:,8), 'VariableNames', varLabels);
        
        writetable(centTabI, [outDir, '/Ci_',parc{p,1},'_',subjectNames{i},'.csv']);
        
    end
    %% save new HubDef
    save([analysis_folder,'2_RCAnalysis/HubDef_Richie2_SingleSubj.mat'],'HubDef');
end
