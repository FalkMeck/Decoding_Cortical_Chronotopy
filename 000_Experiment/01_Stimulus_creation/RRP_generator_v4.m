% Replacement or insertion?
clear;
cd('F:\Arbeit\Programmierung');
freqMat_diff = {[-2 -2 -2 0 -2 -2 -2; -2 -2 -2 0 -2 -2 -2; -2 -2 -2 0 -2 -2 -2; 0 0 0 0 0 0 0; -2 -2 -2 0 -2 -2 -2; -2 -2 -2 0 -2 -2 -2; -2 -2 -2 0 -2 -2 -2],...
                [-2 -6 -2 0 -2 0 -2; -6 0 -6 0 0 0 0; -2 -6 -2 0 -2 0 -2; 0 0 0 0 0 0 0; -2 0 -2 0 -2 -6 -2; 0 0 0 0 -6 0 -6; -2 0 -2 0 -2 -6 -2],...
                [-1 -6 -3 0 -1 0 -1; -5 0 -5 0 0 0 0; -3 -6 -1 0 -1 0 -1; 0 0 0 0 0 0 0; -1 0 -1 0 -1 -6 -3; 0 0 0 0 -5 0 -5; -1 0 -1 0 -3 -6 -1],...
                [0 -6 -4 -2 0 0 -2; -6 0 -6 -2 0 0 0; -4 -6 -2 -2 0 0 0; 0 0 0 0 0 0 0; 0 0 0 -2 -2 -6 -4; 0 0 0 -2 -6 0 -6; -2 0 0 -2 -4 -6 0]};
            % for creating a perfect transistion matrix after violation
            % repalcement
vio = 4; % what is the violation
level = 1; % for which level
true_seq_pos = false; % set to true for "4th" level/true order
paths = ['F:\Arbeit\Programmierung\db_all_1.mat';'F:\Arbeit\Programmierung\db_all_2.mat';'F:\Arbeit\Programmierung\db_all_3.mat'];
% M for Work PC, F for Laptop
true_seq = repmat(transpose([1 2 3 1 2 3 1 2 3 3 2 1 3 2 1 3 2 1 7 6 5 7 6 5 7 6 5 5 6 7 5 6 7 5 6 7]),8,1);
mat_dat = load(paths(level,:)); % loading path to level 2 Data
names = string(fieldnames(mat_dat)); % turning varibale names into strings
%name_of_dat = names(endsWith(names, '_ges')); % find correct string ending with_ges 
pos_beta = 1:287; % just vector from 1:143, that is needed for later
change_positions = cell(6,3); % empty postions, which gets filled wich positons for violations vor each of the 10 sequences
% 1 = positions, 2 + 4 = freqMat replacement, 3 = replaced list, 
% 5 = inserted list, 6 = freqMat insert

cp_col = 0;

for i = [1,3,5] % always balance 36 zeros over two sequences
    cp_col = cp_col +1;    
    if level == 1 && ~ true_seq_pos
        grabs = ones(1, 36); % each transistion ones
        fMd = 1; test_number = 20; % test_numebr is comparision for sum, a bit arbitrary but wirks for insertion better
    elseif level == 2 && ~ true_seq_pos
        alpha_grabs =  [0 3 0 0 0 3 3 0 3 0 0 0 0 0 0 3 0 3 3 0 0 0 3 0]; % only single only grabs
        fMd = 2; test_number = 8;
    elseif level == 3 && ~ true_seq_pos
        alpha_grabs =  [0 3 2 0 0 3 3 2 3 0 0 0 0 0 0 3 2 3 3 0 0 2 3 0]; % singele only and triplet only 
        fMd = 3; test_number = 8;
    elseif true_seq_pos
        grabs = [3 2 1 3 3 2 3 1 1 3 2 3 3 1 2 3]; fMd = 4;
    end 
    
    if ~true_seq_pos
        list1 = eval(['[mat_dat.', char(names), '(:,' char(string(i)),...
            ');mat_dat.', char(names), '(:,' char(string(i+1)),')]']); % take correct 2 sequences
    else
        list1 = true_seq; % for true sequence take true sequence
    end 
   
   level_one_from = list1(1:end-1);
   level_one_to = list1(2:end);
   freqMat_beta = zeros(7,7);

   for fic = 1:length(level_one_from)
       freqMat_beta(level_one_from(fic),level_one_to(fic)) = freqMat_beta(level_one_from(fic),level_one_to(fic)) + 1;
   end
   
   % disp(freqMat_beta);
   
  % Check if one transisiton might be missing ...
   if level == 3 && ~ true_seq_pos % grabs problem on thrid level "solution"
    nonett_ones = [1 4 5 10 11 12 13 14 15 20 21 24];
    check_for_miss = [freqMat_beta(1,1) == 0, freqMat_beta(1,5) == 0, freqMat_beta(1,7) == 0, ...
                      freqMat_beta(3,3) == 0, freqMat_beta(3,5) == 0, freqMat_beta(3,7) == 0, ...
                      freqMat_beta(5,1) == 0, freqMat_beta(5,3) == 0, freqMat_beta(5,5) == 0, ...
                      freqMat_beta(7,1) == 0, freqMat_beta(7,3) == 0, freqMat_beta(7,7) == 0];
    missing_pos = nonett_ones(check_for_miss); % ... and save wich positon is missing
    %if isempty(missing_pos)
    %   missing_pos = 1000;
    %else
    if length(missing_pos) > 1
       error('More than one missing value');
    end
   end
   
    list2 = zeros(287,1); % gets filled wich all transitions
    for ii = 1:287 % gets filled
        list2(ii) = str2double(strcat(string(list1(ii,1)),string(list1(ii+1,1))));
    end
     uni_trans = unique(list2); % saves unique transitions, not every one exists in later levels
     transisions = cell(1, size(uni_trans,1)); % save postions of each diffrent transistion
     for ut = 1:size(uni_trans,1)
        pos_log_ut = list2(:,1) == uni_trans(ut,1);
        pos_ut = pos_beta(pos_log_ut)+1;
        transisions{1,ut}= pos_ut;
     end
     succ_dream = false; % test for comparision with dream freqMat
     count = 0; % counter that it stops at some point
     while  ~succ_dream && count < 10^4
         count = count +1;
         succ_r = true; % random success varibale -> true = failure (a bit confusing ...)
         while succ_r % untill success achived
             % add 12 grabs for triplet level in Level 2
            if level == 2 && ~ true_seq_pos
                grabs = alpha_grabs;
                base_1 = [1 3 4 5]; ran_1 = base_1(sort(randperm(4,3))); % three of each row
                base_2 = [8 10 11 12]; ran_2 = base_2(sort(randperm(4,3)));
                base_3 = [13 14 15 17]; ran_3 = base_3(sort(randperm(4,3)));
                base_4 = [20 21 22 24]; ran_4 = base_4(sort(randperm(4,3)));
                add_grabs = [ran_1, ran_2, ran_3, ran_4];
                
                grabs(add_grabs) = grabs(add_grabs)+1;
                grabs_diff = (-2)*grabs;
                grabs_diff_prep = [grabs_diff(1:3) 0 grabs_diff(4) 0 grabs_diff(5:6) 0 grabs_diff(7) 0 0 0 0 grabs_diff(8:10) 0 grabs_diff(11) 0 grabs_diff(12) 0 0 0 0 0 0 0 grabs_diff(13) 0 grabs_diff(14) 0 grabs_diff(15:17) 0 0 0 0 grabs_diff(18) 0 grabs_diff(19:20) 0 grabs_diff(21) 0 grabs_diff(22:24)];
                freqMat_diff{1,fMd} = transpose(reshape(grabs_diff_prep, 7,7));
            % add 4 grapbs in nonett level in Level 3
            elseif level == 3 && ~ true_seq_pos
                missing_check = true;
                grabs = alpha_grabs;
                while missing_check % only if grabbed number is available
                   base_1 = [1 3 4 5]; ran_1 = base_1(sort(randperm(3,1))); % one of each row
                   base_2 = [8 10 11 12]; ran_2 = base_2(sort(randperm(3,1)));
                   base_3 = [13 14 15 17]; ran_3 = base_3(sort(randperm(3,1)));
                   base_4 = [20 21 22 24]; ran_4 = base_4(sort(randperm(3,1)));
                   add_grabs = [ran_1, ran_2, ran_3, ran_4];
                   missing_check = ismember(missing_pos, add_grabs);
                   possible_problem = ~isempty(intersect(add_grabs, [3 8 17 22]));
                end
                grabs(add_grabs) = grabs(add_grabs)+1;
                grabs_diff = (-2)*grabs;
                grabs_diff_prep = [grabs_diff(1:3) 0 grabs_diff(4) 0 grabs_diff(5:6) 0 grabs_diff(7) 0 0 0 0 grabs_diff(8:10) 0 grabs_diff(11) 0 grabs_diff(12) 0 0 0 0 0 0 0 grabs_diff(13) 0 grabs_diff(14) 0 grabs_diff(15:17) 0 0 0 0 grabs_diff(18) 0 grabs_diff(19:20) 0 grabs_diff(21) 0 grabs_diff(22:24)];
                freqMat_diff{1,fMd} = transpose(reshape(grabs_diff_prep, 7,7));
               
                if length(missing_pos) == 1 % concat for missing_possing correcly (grabs has to be 23, wenn there is a missing otherwise 24)
                   grabs = [grabs(1,1:(missing_pos-1)), grabs(1,(missing_pos+1):end)];
                end
            end
            rand_trans_pos = []; % empty vector whre random postions will get saved in
            if level == 3 && possible_problem
                vio_check = false;
                while ~vio_check
                    rand_trans_pos = [];
                    for r = 1:size(uni_trans,1) % draw random numbers
                        num_trans = numel(transisions{1,r});
                        rand_trans_pos = [rand_trans_pos, transisions{1,r}(1,randperm(num_trans,grabs(1,r)))]; 
                    end
                    vio_check = sum((floor((rand_trans_pos-1)/9) == (rand_trans_pos-1)/9)) == 4;
                end                    
            else
                for r = 1:size(uni_trans,1) % draw random numbers
                    num_trans = numel(transisions{1,r});
                    rand_trans_pos = [rand_trans_pos, transisions{1,r}(1,randperm(num_trans,grabs(1,r)))]; 
                end
            end
            sort_trans_pos = sort(rand_trans_pos); % sort random vector
            cumdiff = zeros((length(sort_trans_pos)-1),1); % 
            for srp = 1:(length(sort_trans_pos)-1) % calculate cumulative difference, to later check that two occurences are not directly following
                cumdiff(srp) = sort_trans_pos(srp+1)-sort_trans_pos(srp); 
            end
            succ_r = min(cumdiff) <= 1 || max(cumdiff) >= 20 || ...
                min(sort_trans_pos) < 3 || ...
                sort_trans_pos(18) > 144 ||sort_trans_pos(19) <= 147; % check condition for success 
        end
         
        freqMat_dream = freqMat_beta + freqMat_diff{1,fMd};
        freqMat_dream(:,4) = [6 6 6 0 6 6 6];
        %freqMat_dream(4,1) = [6 6 6 0 6 6 6];
        if fMd == 2 || fMd == 3
            freqMat_dream(4,:) = [6 6 6 0 6 6 6]; % unsure if that is correct, might cost problems at sum comparision, but therefore the higher number exists
        else 
            freqMat_dream(4,:) = [6 6 6 0 6 6 6];
        end
        freqMat_dream = max(freqMat_dream,0);
        list1_test = list1;
        list1_test(sort_trans_pos,1) = 4; 
        level_one_from = list1_test(1:end-1);
        level_one_to = list1_test(2:end);
        freqMat = zeros(7,7);

        for fi = 1:length(level_one_from)
            freqMat(level_one_from(fi),level_one_to(fi)) = freqMat(level_one_from(fi),level_one_to(fi)) + 1;
        end
        
        sum_dream = sum(abs(freqMat-freqMat_dream),'all'); 
        
        if true_seq_pos % || level == 1
            succ_dream = isequal(freqMat, freqMat_dream); %disp(sum_dream);
        elseif level == 1
            succ_dream = max(max(abs(freqMat-freqMat_dream))) <= 1 && sum_dream <= test_number && min(min(freqMat(:,[1 2 3 5 6 7]))) ~= 0; % a bit arbitrary
            disp([max(max(abs(freqMat-freqMat_dream))) sum_dream]);
            % freqMat
        else
            succ_dream = max(max(abs(freqMat-freqMat_dream))) <= 1 && sum_dream < test_number; % a bit arbitrary
             disp(sum_dream);% freqMat
        end
        
    end % if sum is smaller than test number and and the distance is corrected and thera are no empties
    change_positions{1,cp_col} = sort_trans_pos;
    change_positions{2,cp_col} = freqMat;
% if every thing checks out so far --> contine with rpl- and ins-freqMats    
    
    % REPLACEMENT
    list1_test = list1;
    list1_test(change_positions{1,cp_col},1) = 4; 
    
    level_one_from = list1_test(1:end-1);
    level_one_to = list1_test(2:end);
    freqMat_test = zeros(7,7);

    for fi = 1:length(level_one_from)
       freqMat_test(level_one_from(fi),level_one_to(fi)) = freqMat_test(level_one_from(fi),level_one_to(fi)) + 1;
    end
    disp(i); change_positions{3,cp_col} = list1_test; change_positions{4,cp_col} = freqMat_test; 
    freqMat_test   
    
    % INSERTION
    sort_trans_ins = change_positions{1,cp_col};

    list1_test = [list1(1:(sort_trans_ins(1)-1));vio];
    for ins = 2:36
       list1_test = [list1_test; list1((sort_trans_ins(ins-1)):(sort_trans_ins(ins)-1));vio];
    end
    list1_test =[list1_test;list1(sort_trans_ins(36):end)];


    level_one_from = list1_test(1:end-1);
    level_one_to = list1_test(2:end);
    freqMat_test = zeros(7,7);

    for fi = 1:length(level_one_from)
        freqMat_test(level_one_from(fi),level_one_to(fi)) = freqMat_test(level_one_from(fi),level_one_to(fi)) + 1;
    end
    disp(i); change_positions{5,cp_col} = list1_test; change_positions{6,cp_col} = freqMat_test; 
    freqMat_test
end

%% Save and Export all Change_positions structures
if true_seq_pos
    save(['Change_positions_True_2707.mat'], 'change_positions');
else
    save(['Change_positions_LvL',char(string(fMd)),'_2807.mat'], 'change_positions');
end
% list1_test = true_seq;
%     
%  list1_test= [list1(1:6);4;list1(7:12);4;list1(13:38);4;list1(39:57);4;list1(58:68);4;list1(69:73);4;list1(74:77);4;list1(78:84);4;list1(85:92);...
% 4;list1(93:99);4;list1(100:102);4;list1(103:121);4;list1(122:132);4;list1(133:137);4;list1(138:141);4;list1(142:144);4;list1(145:154);...
% 4;list1(155:157);4;list1(158:164);4;list1(165:167);4;list1(168:170);4;list1(171:174);4;list1(175:186);4;list1(187:191);4;list1(192:203);...
% 4;list1(204:208);4;list1(209:222);4;list1(223:230);4;list1(231:241);4;list1(242:246);4;list1(247:250);4;list1(251:266);4;list1(267:276);...
% 4;list1(277:279);4;list1(280:282);4;list1(283:286);4;list1(287:288)];
% 
% level_one_from = list1_test(1:end-1);
%     level_one_to = list1_test(2:end);
%     freqMat = zeros(7,7);
% 
%     for fi = 1:length(level_one_from)
%        freqMat(level_one_from(fi),level_one_to(fi)) = freqMat(level_one_from(fi),level_one_to(fi)) + 1;
%     end
%     disp(i);
%     freqMat