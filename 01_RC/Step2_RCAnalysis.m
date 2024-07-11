%% 4.2 Hub Analysis Step2
addpath('.../BCT');

subjectNames= {'HFG_121'};

cato_folder = '.../minimally_Processed/DTI/';
cd(cato_folder);

analysis_folder = '.../';

recon_met = 'gqi_dti';
parc = {'lausanne250', 219};

numSubs = numel(subjectNames);

% Random parameters
nRandom = 2500;
rewire = 10;

group_prev_thresh = 0.6;


load([analysis_folder,'/parcInfo.mat']);
load([analysis_folder,'2_RCAnalysis/HubDef_Richie2_SingleSubj.mat']);

% START
for p = 1:size(parc,1)
    nx = parc{p,2};
    
    name = parc{p,1};
    
    ROIs = parcInfo.(name).ROIs;
    regionDescriptions = parcInfo.(name).regionDescriptions;
    
    for i = 1:numSubs
       outDir = fullfile(analysis_folder, '2_RCAnalysis' ,subjectNames{i}, parc{p,1});
       load([outDir, '/connect_subj_bin_', parc{p,1},'_',subjectNames{i},'.mat']);
        
        if ~isempty(connect_subj_bin{1,9})
            regionsSub = []; 
            for l = 'ACDE'
                regionsSub = [regionsSub;HubDef.(name).(l).hubRegions{i}];
            end
            [iSub,jSub,kSub]=unique(regionsSub); % ii = unique elements, jj = position a one/last example of unique lelemnt, kk = reconstruct vector from unique elements (contains information oabout how often an elemt occurs)
            fSub=histc(kSub,1:numel(jSub));
            regions.(parc{p,1}).(subjectNames{i}).Table = table(iSub,fSub);
            regions.(parc{p,1}).(subjectNames{i}).conRegions = iSub(fSub >=3);
            matrixPositions = 1:size(regionDescriptions,1);
            regions.(parc{p,1}).(subjectNames{i}).matrixPos = matrixPositions(ismember(regionDescriptions, regions.(parc{p,1}).(subjectNames{i}).conRegions));
            
            subMat = connect_subj_bin{1,2};
            
            for r = 1:numel(regions.(parc{p,1}).(subjectNames{i}).matrixPos)
                line = regions.(parc{p,1}).(subjectNames{i}).matrixPos(r);
                matrixPositions = 1:size(regionDescriptions,1);
                connectionPositions = matrixPositions(subMat(line,:) == 1);
                degree = sum(subMat(line,:));
                RCconnections = ismember(connectionPositions, regions.(parc{p,1}).(subjectNames{i}).matrixPos);
                RCsum= sum(RCconnections);
                regions.(parc{p,1}).(subjectNames{i}).conRegions(r,1) = regionDescriptions(line);
                regions.(parc{p,1}).(subjectNames{i}).conRegions{r,2} = degree;
                regions.(parc{p,1}).(subjectNames{i}).conRegions{r,3} = RCsum;
                regions.(parc{p,1}).(subjectNames{i}).conRegions{r,4} = RCsum/degree; %proportion of RC connections
                regions.(parc{p,1}).(subjectNames{i}).conRegions{r,5} = RCsum/degree > numel(regions.(parc{p,1}).(subjectNames{i}).matrixPos)/size(regionDescriptions,1);
            end
            if numel(regions.(parc{p,1}).(subjectNames{i}).matrixPos) > 0 % if there are rich clubs
                propRand = zeros(numel(regions.(parc{p,1}).(subjectNames{i}).matrixPos),nRandom); % create a matrix for all random networks
                RandomMats = connect_subj_bin{1,4};
                for rand = 1:nRandom % for all random networks
                    RandomMats_i = RandomMats(:,:,rand); % get the random mat
                    for r = 1:numel(regions.(parc{p,1}).(subjectNames{i}).matrixPos) % for all regions
                        line = regions.(parc{p,1}).(subjectNames{i}).matrixPos(r); % get the line of the region
                        matrixPositions = 1:size(regionDescriptions,1);
                        connectionPositions = matrixPositions(RandomMats_i(line,:) == 1); % for all areas that are RC hubs
                        degree = sum(RandomMats_i(line,:)); % the degree of the region in this random matrix
                        RCconnections = ismember(connectionPositions, regions.(parc{p,1}).(subjectNames{i}).matrixPos);
                        RCsum= sum(RCconnections); % how many are RC connections
                        propRand(r,rand) = RCsum/degree; % proportion of RC connections
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
                save([outDir,'/ROI_mat_',parc{p,1},'_',subjectNames{i},'.mat'],'ROI_mat');
                    
               rD_mat = regionDescriptions(ismember(regionDescriptions, ...
                    regions.(parc{p,1}).(subjectNames{i}).conRegions([regions.(parc{p,1}).(subjectNames{i}).conRegions{:,9}] == 1,1)));
                save([outDir,'/rD_mat_',parc{p,1},'_',subjectNames{i},'.mat'],'rD_mat');
              
               T = cell2table(regions.(parc{p,1}).(subjectNames{i}).conRegions,'VariableNames',...
                    {'Region','Degree','RC_Connections','RC_Proportions','MoreThanProp','MeanRand_RC_Prop','MoreThanRand','MedianDiff', 'WilcoxonTest'});
                writetable(T,[outDir,'/Regions_',parc{p,1},'_',subjectNames{i},'.csv']);
            end
        end
    end
end
