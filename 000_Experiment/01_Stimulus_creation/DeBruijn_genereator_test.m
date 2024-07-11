% W. Owen Brimijoin (2021). de Bruijn sequence generator (https://www.mathworks.com/matlabcentral/fileexchange/28015-de-bruijn-sequence-generator), MATLAB Central File Exchange. Retrieved May 5, 2021.
clear;
for all = 6 %2:6
    db_all = zeros(36,4); % empty frame for all 4 
    warning('off','all'); % so the dbgenerator doesn't show the warning each time
    for db = 1:size(db_all,2) % for all 4 sequences fpr one 144 sequence
        if db > 1 % creates 1st Seqeunce withiout checking
            general = true; % Defintions for later
            counter = 0; %to count to 1.000.000
            comp_beta = 1000; % just a high number for the sum comp later
            db_best_so_far = 0;
            while general && counter <= 10^6 % general is true, if Sequence generation failed
                counter = counter +1; % counting repetiotions
                db_all(:,db) = transpose(debruijn_generator(6,2)); % generate
                db_new = db_all(:,db); % save dor comparison 
                total_circ_test = []; total_comp_test = [];
                if db == 4 %for 4th repetion compare with 2+3 and nocht 1,2+3, since it seemed impossible to find a solotion then
                    k_start = 2; 
                else 
                    k_start = 1; 
                end
                for k = k_start:(db-1) % starting from 1 or 2 respectively
                    db_old = db_all(:,k); % compare old (one of the already genreated ones) and teh newely genreated order
                    db_equal = NaN((size(db_new,2)-1),1); %to calculate all times the two are equela when the one is movesd one space at a time
                    for i = 1:(size(db_old,2)-1) % circle test
                        db_shift = circshift(db_new,i);
                        db_equal(i,1) = isequal(db_old, db_new);
                    end
                    circ_test = sum(db_equal) > 0; % circ_test = ture is bad

                    comp = 0;
                    for l = 3:9 %test for same 3er or 9er seqeunces within the two sequence --> thus no similar structures (rarely fins something over 4er, 9er enough)
                        for i = 1:(size(db_old,1)-(l-1)) % every piece old
                            for j = 1:(size(db_new, 1)-(l-1)) % vs every piece new
                                equality = isequal(db_old(i:(i+(l-1)),1),db_new(j:(j+(l-1)),1)); % search for equality
                                if equality
                                    comp = comp +1*(l/3); % if equal is found add to comp counter (with weigh, but it doesn't actually matter, because it is compared with 0)
                                end
                            end
                        end 
                    end
                    comp_test = comp > 0; % comp test = true = bad
                    total_circ_test = [total_circ_test,circ_test]; % bind all circ test for all old vs. new comparisions together
                    total_comp_test = [total_comp_test,comp_test]; % fro all comp tests as well 
                    if sum(total_comp_test) < comp_beta && sum(total_circ_test) == 0 %check if ok and save best sequence so far
                        db_best_so_far = db_new;
                        comp_beta = sum(total_comp_test); 
                    end
                end
            general = sum(total_circ_test) > 0 || sum(total_comp_test) > 0; % check if it worked, general = true is bad
            if counter == 10^6 && general == true % if it tried 1.000.000 times, use best so far, but never happened
                db_all(:,db) = db_best_so_far;
            end 
            end
        else 
            db_all(:,db) = transpose(debruijn_generator(6,2));% for first order
        end
    end
save(['db_all_1_',char(string(all)),'.mat'],'db_all'); % save this matrix
end

%% save the matrices (after all 6 versions has been produced)
db_LvL1 = zeros(144,6); 
for i = 1:6
    db_all_1_beta = load(['db_all_1_',char(string(i)),'.mat']);
    db_LvL1(:,i) = [db_all_1_beta.db_all(:,1);db_all_1_beta.db_all(:,2);db_all_1_beta.db_all(:,3);db_all_1_beta.db_all(:,4)]; % bind all 4 together into one vector (144 x 1)
    for j = 1:size(db_LvL1,1)
        if db_LvL1(j,i) >=4 % rename 
           db_LvL1(j,i) = db_LvL1(j,i)+ 1;
        end
    end
end

%% Save Matrix of all together as one file
save('db_all_1.mat','db_LvL1');
writematrix(db_LvL1, 'E:\Arbeit\Programmierung\order_1_1405.txt','Delimiter','tab');