%% CATO connectivity Matrices Outlier control

% The OUTLIER-IDENTIFICATION is based on the Paper 10Kin1day: A Bottom-Up Neuroimaging Initiative
% Part for Outlier Identification:

%         - Connectome Reconstruction -
%         A connectome map was made by combining the (sub)cortical
%         parcellation map and the set of reconstructed fibers using
%         commonly described procedures [see (18–21)]. For each of
%         the Cammoun Desikan-Killiany parcellation maps (i.e., 14+68,
%         14+114, and 14+219 regions, respectively), the total collection
%         of reconstructed fiber streamlines was used to assess the level
%         of connectivity between each pair of (sub)cortical regions,
%         represented as the connectivity matrix CIJ. (Sub)cortical regions
%         were selected as the nodes of the reconstructed network, and
%         for each combination of region i and region j where fiber
%         streamlines touched both regions a connection (i.e., network
%         edge) was included in cell CIJ(i,j) in the connectivity matrix.
%         Five different types of strength of a connection were computed
%         and included as edge strength: (1) the number of reconstructed
%         streamlines (NOS) between region i and j, (2) the average
%         FA of the voxels traversed by the reconstructed streamlines,
%         (3) the average MD of the reconstructed streamlines, (4)
%         the average length of the reconstructed streamlines and (5)
%         streamline density computed as the number of reconstructed
%         streamlines corrected for the average volume of region i and
%         region j (18, 19).
% 
%         - Outliers -
%         A total of 15,947 connectome maps were analyzed across the
%         participating groups. Of the datasets that could be shared, 197
%         were detected as outliers (and were subsequently removed from
%         the dataset). Outliers were detected automatically by testing per
%         dataset and for each connectome map their average connection
%         strength and their distance to the group average prevalence
%         map. The average connection strength of a connectome map
%         was calculated for each of the five connection weights as the
%         mean of the strengths over all existing (nonzero) connections.
%         To measure the presence of odd connections or absence of
%         common connections in a connectome map, we constructed a
%         group prevalence matrix for each dataset, counting per node
%         pair how many times an edge was observed across the group
%         of subjects in the dataset. For each connectome map the total
%         prevalence of all observed connections and the total prevalence
%         of all non-observed connections was computed. Outliers were
%         identified as connectome maps that displayed on any of the 7
%         measures (5 weight and 2 prevalence measures) a score below Q1
%         – 2×IQR or above Q3 + 2×IQR, with Q1 and Q3 referring to
%         the first and third quartile, respectively and IQR the interquartile
%         range IQR = Q3 – Q1. This resulted in the detection of 189
%         outliers in total, which were excluded from the dataset. One
%         complete dataset (set_634413, n=584) showed across all included
%         individual sets an average lower FA / higher MD as compared
%         to the other datasets and this set was excluded from the age
%         curves shown in Figure 1. Due to the high overall sample size,
%         including or excluding this dataset did not change the shape of
%         the final plot


%% STARTING VARIABLES
clearvars; clc;

subject_dir = ['HFG_121'];

number_subjects = size(subject_dir,1); 

cato_dir = '...\DTI\Day2_withoutSubcortical\';

parcels = {'lausanne250'};
methods = {'gqi_dti'};

%% CREATE PREVALENCE MAPS

prevalence_maps = cell(size(parcels,2), size(methods,2));
for parc = 1:size(parcels,2)
    for met = 1:size(methods,2)
        for i = 1:number_subjects
            connect = load([cato_dir, subject_dir(i,:), '/DWI_processed/', subject_dir(i,:), '_connectivity_', methods{1,met}, '_' ,parcels{1,parc}, '.mat']); 
            nonzero = connect.connectivity(:,:,1) ~= 0;  
            if i == 1
                prev = zeros(size(nonzero)); 
            end
            prev = prev + nonzero; 
        end
        prevalence_maps{parc, met} = prev; 
    end    
end


%% DETECT OUTLIERS PREPARATION
detect = cell(size(subject_dir,1), size(parcels,2)+ 1); 

for i = 1:number_subjects
    detect{i,1} = subject_dir(i,:); 
    for parc = 1: size(parcels,2)
        properties = load([cato_dir, subject_dir(i,:), '/DWI_processed/', subject_dir(i,:), '_region_properties_', parcels{1,parc}, '.mat']); 
        detect_method = cell(size(methods,2),1+1+5+2+2);
        for met = 1:size(methods, 2)
            connect = load([cato_dir, subject_dir(i,:), '/DWI_processed/', subject_dir(i,:), '_connectivity_', methods{1,met}, '_' ,parcels{1,parc}, '.mat']); 
            nonzero = connect.connectivity(:,:,1) ~= 0;
            found_SL = sum(nonzero, 'all'); % number of Stremlines
            detect_method{met,1} = parcels{1,parc}; % parcellation
            detect_method{met,2} = methods{1,met}; % methods
            detect_method{met,3} = sum(connect.connectivity(:,:,1),'all')/found_SL; % mean NOS %does not matter if squareform or not, since both are x2 then
            detect_method{met,4} = sum(connect.connectivity(:,:,3),'all')/found_SL; % FA % and after that those 5 values are just one number
            detect_method{met,5} = sum(connect.connectivity(:,:,6),'all')/found_SL; % MD % which is the same regrardless of squareform()
            detect_method{met,6} = sum(connect.connectivity(:,:,2),'all')/found_SL; % length
            
            % Stremline densitiy corrected for avarage volume (not 13, but
            % corrected by number voxels, highly correlated)
            new13 = zeros(size(connect.connectivity(:,:,1)));  
            for x = 1:size(new13,1)
                RegX = properties.regionPropertiesTable.number_of_voxels(x);
                for y = 1:size(new13,2)
                    RegY = properties.regionPropertiesTable.number_of_voxels(y);
                    new13(x,y) = connect.connectivity(x,y,1)/((RegX + RegY)/2);
                end
            end
            detect_method{met,7} = sum(new13,'all')/found_SL; % Streamline desitiy avg. for voxel number/volume
            
            % Prevalence Map comparision
            % With sqaurefrom, sums of prev_curr have to be diveded by 2
            % But it is just a linar transforamtion, Median and IQR etc.
            % are just mutliplied by 2, no effect on distribution
            prev_curr = prevalence_maps{parc,met};
            detect_method{met,8} = prev_curr.*nonzero;  % prevelance map of all found connections
            detect_method{met,9} = sum(prev_curr(nonzero),'all');%/2;      % total prevalence of found connections
            detect_method{met,10}= prev_curr.*~nonzero; % prevelance map of all non-found connectiosn
            detect_method{met,11}= sum(prev_curr(~nonzero),'all');%/2;     % total prevalence of non-found connections      
        end
        detect{i,parc+1} = detect_method;
    end
end

%% OUTLIER VIA IQR

measures = {'NOS', 'FA', 'MD', 'Length', 'SL_Density', 'P_found_tot', 'P_non_found_tot'};
% all measures used for Outlier check

% Restructure dataframe to have all Subjects in a structured frame
% structred for parcealltion and then SL-identification Method
for parc = 1:size(parcels, 2) 
    for met = 1:size(methods,2) 
        A = cell(number_subjects, size(measures,2) +1);
        for i = 1:number_subjects
            A{i,1} = subject_dir(i,:); 
            A{i,2} = detect{i,parc+1}{met,3};
            A{i,3} = detect{i,parc+1}{met,4};
            A{i,4} = detect{i,parc+1}{met,5};
            A{i,5} = detect{i,parc+1}{met,6};
            A{i,6} = detect{i,parc+1}{met,7};
            A{i,7} = detect{i,parc+1}{met,9};
            A{i,8} = detect{i,parc+1}{met,11};
        end
        checkout.(parcels{1,parc}).(methods{1,met}) = A; 
    end
end

% Finding Outlier by IQR 
% Produces a Structure, strucutred for parcellation, method and mesaurement
% gives name of Suject found as outlier
% does not give direction of outlier... MAYBE CAHNGE?
for parc = 1:size(parcels, 2) 
    for met = 1:size(methods,2) 
        for mes = 1:size(measures,2)
            mesA = [checkout.(parcels{1,parc}).(methods{1,met}){:,mes+1}];
            upper = quantile(mesA, .75) + 2*iqr(mesA); 
            lower = quantile(mesA, .25) - 2*iqr(mesA); 
            out_mesA = mesA > upper | mesA < lower; 
            checkIQR.(parcels{1,parc}).(methods{1,met}).(measures{1,mes}) = ...
                checkout.(parcels{1,parc}).(methods{1,met})(out_mesA,1);
        end
    end
end
%%
save('Outlier_10kin1day.mat')