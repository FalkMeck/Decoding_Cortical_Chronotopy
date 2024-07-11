%% ATLAS OVERLAP Between Glasser atlas and lausanne
spm_dir = '...\spm12\';
addpath(genpath(spm_dir));

%% Coregister a MNI template to the volumetric MMP atlas and use multiscle brain parcelator to work with that
% This templated was then parcelled according to DK/lausanne120-500
% parcellations using MSBP (Toubrier, 2019)

niiTools = '...\NIfTI_20140122';
addpath(niiTools);

msbp_dir = '...\msbp_anat\';
cd(msbp_dir); 

MMP_dir = '...\MNI_2009c_asym\';
MMP_atlas = fullfile(MMP_dir, 'MMP_in_MNI_corr.nii'); 
MMP_atlas_nii = load_untouch_nii(MMP_atlas); 

Glasser_plus = readtable('...\Glasser_plus.csv'); 
  
parc = {'scale1', 'scale2', 'scale3', 'scale4';...
    'aparc', 'lausanne120', 'lausanne250', 'lausanne500'};

NetworksGlasser = unique(Glasser_plus.Network);
Networks = [{'Region'}, NetworksGlasser',{'No_overlap'}];
OverviewNames = {'Region','Network','Network_percentage','Ji_Ito','Ji_Ito_adj','Ji_Ito_multi_perc'};

Ji_Ito_uni = {'Auditory', 'Visual1', 'Visual2', 'Somatomotor'};
Ji_Ito_multi = {'Cingulo-Opercular','Default','Dorsal-Attention','Frontoparietal',...
                'Language','Orbito-Affective','Posterior-Multimodal','Ventral-Multimodal'};
Ji_Ito_multi_adj = {'Cingulo-Opercular','Default','Frontoparietal'};

Trans_Networks = ismember(NetworksGlasser, Ji_Ito_multi); 

for p = 3

    scalePath = [msbp_dir, 'sub-01_label-L2018_desc-',parc{1,p},'_atlas.nii'];
    translator = readtable([msbp_dir, 'transS',parc{1,p}(2:end), '.csv']);
    
    scale_nii = load_untouch_nii(scalePath); 
    
    overlap_lausanne_Glasser = cell(size(translator,1),14);
    overlap_lausanne_Glasser(:,2:end) = {0};
    overlap_lausanne_Glasser = cell2table(overlap_lausanne_Glasser, 'VariableNames', Networks); 
    
    overview_lausanne_Glasser = cell(size(translator,1),numel(OverviewNames));
    overview_lausanne_Glasser = cell2table(overview_lausanne_Glasser, 'VariableNames', OverviewNames);
    
    for i = 1:size(translator,1)
        overlap_lausanne_Glasser.Region(i) = translator.region(i);
        overview_lausanne_Glasser.Region(i) = translator.region(i);
        overlap = MMP_atlas_nii.img(scale_nii.img == translator.msbpROI(i)); 
        uniOverlap = unique(overlap)';
        num_voxels_lausanne = size(overlap,1);
        for ii = 1:length(uniOverlap)
            percentage = sum(overlap == uniOverlap(ii))/num_voxels_lausanne;
            if uniOverlap(ii) ~= 0
                network_ii = Glasser_plus.Network{Glasser_plus.Num_asym == uniOverlap(ii)}; 
            else
              network_ii = Networks{end}; 
            end
            
            overlap_lausanne_Glasser.(network_ii)(i) = overlap_lausanne_Glasser.(network_ii)(i) + percentage;
            
        end
        
        network_percentages = table2array(overlap_lausanne_Glasser(i,2:13));
        
        find_best_network = [false,network_percentages == max(network_percentages),false];
        best_network = Networks(find_best_network);
        
        overview_lausanne_Glasser.Network(i) = best_network;
        overview_lausanne_Glasser.Network_percentage{i} = max(table2array(overlap_lausanne_Glasser(i,2:13)));
        
        if ismember(best_network,Ji_Ito_multi_adj) && ismember(best_network,Ji_Ito_multi)
            overview_lausanne_Glasser.Ji_Ito_adj{i} = 'Core';
            overview_lausanne_Glasser.Ji_Ito{i} = 'Transmodal';
        elseif ismember(best_network,Ji_Ito_multi) && ~ismember(best_network,Ji_Ito_multi_adj)
            overview_lausanne_Glasser.Ji_Ito_adj{i} = 'NoLabel';
            overview_lausanne_Glasser.Ji_Ito{i} = 'Transmodal';
        elseif ismember(best_network, Ji_Ito_uni)
            overview_lausanne_Glasser.Ji_Ito_adj{i} = 'Periphery';
            overview_lausanne_Glasser.Ji_Ito{i} = 'Unimodal';
        end
        
        overview_lausanne_Glasser.Ji_Ito_multi_perc{i} = sum(network_percentages(Trans_Networks));
        
    end   
    
    writetable(overlap_lausanne_Glasser, fullfile(MMP_dir,['overlap_',parc{2,p},'_Glasser.xlsx']));
    writetable(overview_lausanne_Glasser, fullfile(MMP_dir,['overview_',parc{2,p},'_Glasser.xlsx']));
end
