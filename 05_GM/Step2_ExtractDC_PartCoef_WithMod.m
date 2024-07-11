%% Extract partCoef and WithinMod degZ from from the repPart used in DC 

clear all;
addpath('.../BCT'); %Brain Connecticity Toolbox
clc;

subjectNames= {'HFG_121'};

analysisFolder = '...';

RCanalysisFolder = fullfile(analysisFolder, '2_RCAnalysis');


parc = {'lausanne250', 219};

parcInfoPath = fullfile(RCanalysisFolder,'parcInfo.mat');
load(parcInfoPath);

for p = 1:size(parc,1)
    disp(['Parcellation scheme: ', parc{p,1}]);
    regDesc = parcInfo.(parc{p,1}).regionDescriptions;
     disp('Start Single Subjects...');
     for i = 1:numel(subjectNames)
        disp(['Subject: ', subjectNames{i}]);         
        % Get connect matrix
       conPath = fullfile(RCanalysisFolder,subjectNames{i},parc{p,1},['connect_subj_bin_', parc{p,1},'_',subjectNames{i},'.mat']);
         load(conPath);
         
         conMat = connect_subj_bin{1,2};
        
        repPartPath = fullfile(RCanalysisFolder,subjectNames{i},parc{p,1},['repPart_', parc{p,1},'_',subjectNames{i},'.mat']);
        load(repPartPath);
         
        results = cell(numel(regDesc),3); 
        
        results(:,1) = regDesc;
        results(:,2) = num2cell(participation_coef(conMat,repPart_gamma));
        results(:,3) = num2cell(module_degree_zscore(conMat,repPart_gamma)); 
         
        resTab = cell2table(results, 'VariableNames', {'Region', 'PartCoef_DC', 'WithinMod_DC'}); 
        
        outPath = fullfile(RCanalysisFolder,subjectNames{i},parc{p,1},['DCi_',parc{p,1},'_',subjectNames{i},'.csv']);
        writetable(resTab, outPath); 
        disp(['Done DC ', num2str(i),'/',num2str(numel(subjectNames))]);
%         
     end
%     
