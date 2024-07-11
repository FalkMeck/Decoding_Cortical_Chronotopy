%% Diverse HUB ANALYSIS RICHIE

clear all;
clc;
addpath('.../BCT');
addpath('.../04_DC');
addpath('.../01_RC');

subjectNames= {'HFG_121'};  

cato_folder = '.../minimally_Processed/DTI/';
cd(cato_folder);

analysis_folder = '.../';
%
parc = {'lausanne250', 219};

numSubs = numel(subjectNames);

recon_met = 'gqi_dti';

% Random parameters
nRandom = 2500;
rewire = 10;

group_prev_thresh = 0.6;

partitions = 1000;
gamma = 0.5:0.1:3.5;

single_or_group = 15; % 15 = single, 17 = group; % use singel subjetc or group matrices
% I used single subject

load([analysis_folder,'2_RCAnalysis/parcInfo.mat']);

connectVars = {'Subj', 'connectMat', 'DC', 'randomConnectMat', 'randomDC', 'normDC'};

%% 1. get matrices
for p = 1:size(parc,1)
    % Option 1: Group (did not use this)
    nx = parc{p,2};
    load([analysis_folder,'2_RCAnalysis/', 'groupMatrix_', parc{p,1},'.mat']);

    [repPart_group, ~] = consensual_partition(groupMatrix, gamma, partitions); 
    save([analysis_folder,'2_RCAnalysis/', 'repPart_group_', parc{p,1},'_prev.mat'], 'repPart_group'); 
    
    partCoeff_group = participation_coef(groupMatrix, repPart_group);
    
    for i = 1:numSubs
        %Option 2: Subject-wise
        outDir = fullfile(analysis_folder, '2_RCAnalysis' ,subjectNames{i}, parc{p,1});
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        
        load([outDir,'/connect_subj_bin_', parc{p,1},'_',subjectNames{i},'.mat']);
        diverse_subj_bin = connect_subj_bin;
        disp([parc{p,1},': ',subjectNames{i}]);
        
        %% 2. Random networks per subject and calculate DC coeffient in both true and random netowrks
        Cwm = diverse_subj_bin{1,2}; %get wighted connectiivty matrix
        
        [repPart_gamma, ~] = consensual_partition(Cwm, gamma, partitions); 
        save([outDir,'/repPart_', parc{p,1},'_',subjectNames{i},'.mat'], 'repPart_gamma'); 
               
        partCoeff = participation_coef(Cwm, repPart_gamma);
        diverse_subj_bin{1,15} = partCoeff;
        diverse_subj_bin{1,16} = unique(partCoeff);
        
        diverse_subj_bin{1,17} = partCoeff_group;
        diverse_subj_bin{1,18} = unique(partCoeff_group);
        
        % With single subject part Coeff
        diverse_subj_bin{1,3} = all_club_bu(Cwm, diverse_subj_bin{1,single_or_group}); % calculate weighted DC Coefficient
        
        randomConMat = diverse_subj_bin{1,4};
        randomDC = cell(nRandom,1);
        for iRan = 1:nRandom % randomize network 2500 times
            % With single subject part Coef
            randomDC{iRan,1} = all_club_bu(randomConMat(:,:, iRan),diverse_subj_bin{1,single_or_group}); % calculate weighted DC Coefficient
        end
        diverse_subj_bin{1,5} = randomDC;
        
        %% 3. DC normalization
        DC = diverse_subj_bin{1,3}; % true DC Coefficients
        DCRand = cell2mat(randomDC(:,1));
        DCRandMean = mean(DCRand,1); % mean random DC Coefficent
        diverse_subj_bin{1,6} = DC./DCRandMean; % normalized DC Coefficient
        
        % find k differences
        kComp = zeros(3,numel(DC));
        for ki = 1:numel(DC)
            kComp(1,ki) = ki;
            
            kComp_i = DC(ki);
            kComp(2,ki) = kComp_i;
            kComp_Rand = DCRand(:,ki);
            pvalue = sum(kComp_Rand >= kComp_i)/nRandom;
            
            kComp(3,ki) = pvalue;
        end
        
        diverse_subj_bin{1,7} = kComp;
        
        kCompCorrected = kComp(:, ~isnan(DC));
        [kCompCorrected(5,:),~,~, kCompCorrected(4,:)] = fdr_bh(kCompCorrected(3,:),0.05,'dep','yes');
        
        diverse_subj_bin{1,8} = kCompCorrected;
        
        %% 3.5 Working with the identified range of k
        % identify largest range of subsequet significant k
        % 1st case: highest and lowest to define range
        % 2nd case: intermitted by 1/ 2 non-neughboring p-values
        % accept as outlier
        % 3rd case: intermitted by more neighboring p-values
        % 1. 1.5 range of other (longest range)
        % 2. highest k in range
        % 3. exclude
        
        
        k_range = diverse_subj_bin{1,8};
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
                    max_k_group = zeros(1,numel(grouped));% case does not happen, dont waste time right now
                    for gr = 1:numel(grouped)
                        max_k_group(gr) = max(diverse_subj_bin{1,6}(grouped{gr}));
                    end
                    kRangeTrue = grouped{max_k_group == max(max_k_group)};
                    diverse_subj_bin{1,9} =  k_range(2,kRangeTrue); % take range that includes maximal normalizes clubness
                else
                    kRangeTrue = grouped{size_grouped == max_g};
                    diverse_subj_bin{1,9} =  k_range(2,kRangeTrue); % take longest
                end
            else
                kRangeTrue = min(kRangeSig):max(kRangeSig);
                diverse_subj_bin{1,9} = k_range(2,kRangeTrue);  % ignore breaks
            end
        else
            kRangeTrue = kRangeSig;
            diverse_subj_bin{1,9} = k_range(2,kRangeTrue); % if just 1 group, take whole group
        end
        
        diverse_subj_bin{1,10}  = length(diverse_subj_bin{1,9}); % length of DC regime
        diverse_subj_bin{1,11} = diverse_subj_bin{1,single_or_group+1}(min(kRangeTrue)); % min border % start of regime
        diverse_subj_bin{1,12} = diverse_subj_bin{1,single_or_group+1}(max(kRangeTrue)); % max border
        diverse_subj_bin{1,13} = max(diverse_subj_bin{1,6}(kRangeTrue)); % peak of norm in range
        if ~isempty(diverse_subj_bin{1,13})
            diverse_subj_bin{1,14} = diverse_subj_bin{1,single_or_group+1}(find(diverse_subj_bin{1,6} == diverse_subj_bin{1,13},1)); % peak regime position
        else
            diverse_subj_bin{1,14} = [];
        end
        save([outDir,'/diverse_subj_bin_', parc{p,1},'_',subjectNames{i},'.mat'], 'diverse_subj_bin',  '-v7.3');
    end
end
%% 4.1 Hub Analysis Step 1

for p = 1:size(parc,1)
    nx = parc{p,2};
    hubNum = floor(nx*.20);
    
    name = parc{p,1};
    
    ROIs = parcInfo.(name).ROIs;
    regionDescriptions = parcInfo.(name).regionDescriptions;
    
    for i = 1:numSubs
        outDir = fullfile(analysis_folder, '2_RCAnalysis' ,subjectNames{i}, parc{p,1});
        
        load([outDir,'/diverse_subj_bin_', parc{p,1},'_',subjectNames{i},'.mat']);
        
        subMat = diverse_subj_bin{1,2};
        
        %% a) top 20% participation coefficient distribution
        HubDef.(name).A.Subj(i) = subjectNames(i);
        HubDef.(name).A.partCoeff{i} = diverse_subj_bin{1,single_or_group};
        [HubDef.(name).A.partCoeffSort{i},...
            HubDef.(name).A.partCoeffSortIdx{i}] = sort(HubDef.(name).A.partCoeff{i}, 'descend');
        Cutoff = HubDef.(name).A.partCoeffSort{i}(hubNum);
        %         HubDef.(name).A.hubROIs{i} = ROIs(HubDef.(name).A.partCoeffSortIdx{i}(1:hubNum));
        %         HubDef.(name).A.hubRegions{i} = regionDescriptions(HubDef.(name).A.partCoeffSortIdx{i}(1:hubNum));
        HubDef.(name).A.hubROIs{i} = ROIs(HubDef.(name).A.partCoeff{i} >= Cutoff);
        HubDef.(name).A.hubRegions{i} = regionDescriptions(HubDef.(name).A.partCoeff{i} >= Cutoff);
        
        
        %% c) nodal degree >= starting point of DC regime
        
        if ~isempty(diverse_subj_bin{1,9})
            startDCregime = diverse_subj_bin{1,11};
            HubDef.(name).C.Subj(i) = subjectNames(i);
            HubDef.(name).C.partCoeff{i} = diverse_subj_bin{1,single_or_group};
            HubDef.(name).C.hubROIs{i} = ROIs(HubDef.(name).C.partCoeff{i} >= startDCregime);
            HubDef.(name).C.hubRegions{i}= regionDescriptions(HubDef.(name).C.partCoeff{i} >= startDCregime);
        else
            disp(subjectNames{i});
        end
        
        %% d) nodal degree >= peak of normalized DC coefficient
        
        if ~isempty(diverse_subj_bin{1,9})
            peakDCregime = min(diverse_subj_bin{1,14});
            HubDef.(name).D.Subj(i) = subjectNames(i);
            HubDef.(name).D.partCoeff{i} = diverse_subj_bin{1,single_or_group};
            HubDef.(name).D.hubROIs{i} = ROIs(HubDef.(name).D.partCoeff{i} >= peakDCregime);
            HubDef.(name).D.hubRegions{i} = regionDescriptions(HubDef.(name).D.partCoeff{i} >= peakDCregime);
        else
            disp(subjectNames{i});
        end
        
    end
    %% save new HubDef
    save([analysis_folder,'2_RCAnalysis/HubDef_Richie2_DC_Subj.mat'],'HubDef');
end


%% 4.2 Hub Analysis Step2

% identify the regions based on the hub criterion
% Bertolero et al. (2017) use a 80% cut off as for the diverse (and rich)
% clubs, as the normalized Club coefficients increases drastically after 80%


% START
for p = 1:size(parc,1)
    nx = parc{p,2};
    
    name = parc{p,1};
    
    ROIs = parcInfo.(name).ROIs;
    regionDescriptions = parcInfo.(name).regionDescriptions;
    
    for i = 1:numSubs
        outDir = fullfile(analysis_folder, '2_RCAnalysis' ,subjectNames{i}, parc{p,1});
        if ~isempty(HubDef.(name).C.Subj{1,i})
            regions.(parc{p,1}).(subjectNames{i}).conRegions = HubDef.(name).A.hubRegions{i};
            matrixPositions = 1:size(regionDescriptions,1);
            regions.(parc{p,1}).(subjectNames{i}).matrixPos = matrixPositions(ismember(regionDescriptions, regions.(parc{p,1}).(subjectNames{i}).conRegions));
            
            load([outDir,'/diverse_subj_bin_', parc{p,1},'_',subjectNames{i},'.mat']);
            
            subMat = diverse_subj_bin{1,2};
            
            for r = 1:numel(regions.(parc{p,1}).(subjectNames{i}).matrixPos)
                line = regions.(parc{p,1}).(subjectNames{i}).matrixPos(r);
                matrixPositions = 1:size(regionDescriptions,1);
                connectionPositions = matrixPositions(subMat(line,:) == 1);
                degree = sum(subMat(line,:));
                DCconnections = ismember(connectionPositions, regions.(parc{p,1}).(subjectNames{i}).matrixPos);
                DCsum= sum(DCconnections);
                regions.(parc{p,1}).(subjectNames{i}).conRegions(r,1) = regionDescriptions(line);
                regions.(parc{p,1}).(subjectNames{i}).conRegions{r,2} = degree;
                regions.(parc{p,1}).(subjectNames{i}).conRegions{r,3} = DCsum;
                regions.(parc{p,1}).(subjectNames{i}).conRegions{r,4} = DCsum/degree; %proportion of DC connections
                regions.(parc{p,1}).(subjectNames{i}).conRegions{r,5} = DCsum/degree > numel(regions.(parc{p,1}).(subjectNames{i}).matrixPos)/size(regionDescriptions,1);
            end
            if numel(regions.(parc{p,1}).(subjectNames{i}).matrixPos) > 0 % if there are rich clubs
                propRand = zeros(numel(regions.(parc{p,1}).(subjectNames{i}).matrixPos),nRandom); % create a matrix for all random networks
                RandomMats = diverse_subj_bin{1,4};
                for rand = 1:nRandom % for all random networks
                    RandomMats_i = RandomMats(:,:,rand); % get the random mat
                    for r = 1:numel(regions.(parc{p,1}).(subjectNames{i}).matrixPos) % for all regions
                        line = regions.(parc{p,1}).(subjectNames{i}).matrixPos(r); % get the line of the region
                        matrixPositions = 1:size(regionDescriptions,1);
                        connectionPositions = matrixPositions(RandomMats_i(line,:) == 1); % for all areas that are DC hubs
                        degree = sum(RandomMats_i(line,:)); % the degree of the region in this random matrix
                        DCconnections = ismember(connectionPositions, regions.(parc{p,1}).(subjectNames{i}).matrixPos);
                        DCsum= sum(DCconnections); % how many are DC connections
                        propRand(r,rand) = DCsum/degree; % proportion of DC connections
                        regions.(parc{p,1}).(subjectNames{i}).conRegions{r,6} = mean(propRand(r,:),2); % mean of proportion
                        regions.(parc{p,1}).(subjectNames{i}).conRegions{r,7} = regions.(parc{p,1}).(subjectNames{i}).conRegions{r,6} <= regions.(parc{p,1}).(subjectNames{i}).conRegions{r,4};
                    end
                end
                
                propDiff = repmat([regions.(parc{p,1}).(subjectNames{i}).conRegions{:,4}],nRandom,1)' -propRand; % difference of random and true value
                pValuesWilcox = zeros(numel(regions.(parc{p,1}).(subjectNames{i}).matrixPos),1);
                for r = 1:numel(regions.(parc{p,1}).(subjectNames{i}).matrixPos) % for all regions
                    regions.(parc{p,1}).(subjectNames{i}).conRegions{r,8} = median(propDiff(r,:));
                    [pValuesWilcox(r,1), ~] = signrank(propDiff(r,:), [],'tail', 'right'); %wilcoxon signed rank test
                end
                [h,~,~,~] = fdr_bh(pValuesWilcox);
                regions.(parc{p,1}).(subjectNames{i}).conRegions(:,9) = num2cell(h);
                
                ROI_mat = ROIs(ismember(regionDescriptions, ...
                    regions.(parc{p,1}).(subjectNames{i}).conRegions([regions.(parc{p,1}).(subjectNames{i}).conRegions{:,9}] == 1,1)));
                save([outDir,'/DC_ROI_mat_',parc{p,1},'_',subjectNames{i},'.mat'],'ROI_mat');
                
                rD_mat = regionDescriptions(ismember(regionDescriptions, ...
                    regions.(parc{p,1}).(subjectNames{i}).conRegions([regions.(parc{p,1}).(subjectNames{i}).conRegions{:,9}] == 1,1)));
                save([outDir,'/DC_rD_mat_',parc{p,1},'_',subjectNames{i},'.mat'],'rD_mat');
                
                T = cell2table(regions.(parc{p,1}).(subjectNames{i}).conRegions,'VariableNames',...
                    {'Region','Degree','DC_Connections','DC_Proportions','MoreThanProp','MeanRand_DC_Prop','MoreThanRand','MedianDiff', 'WilcoxonTest'});
                writetable(T,[outDir,'/DC_Regions_',parc{p,1},'_',subjectNames{i},'.csv']);
            end
        end
    end
end
