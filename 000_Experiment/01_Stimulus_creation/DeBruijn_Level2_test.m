% W. Owen Brimijoin (2021). de Bruijn sequence generator (https://www.mathworks.com/matlabcentral/fileexchange/28015-de-bruijn-sequence-generator), MATLAB Central File Exchange. Retrieved May 5, 2021.
%db_total_LvL2 = zeros(48,6);
%fails = []; 
%for all = 2:6 % also did it one by one
all = 6;
db_all_2_beta = zeros(16,3);
warning('off','all');
for db = 1:size(db_all_2_beta,2)
if db > 1 
    general = true;
    counter = 0;
    comp_beta = 1000; 
    db_best_so_far = 0;
    while general && counter <= 10^6
        counter = counter +1; % count to 1.000.000
        db_all_2_beta(:,db) = transpose(debruijn_generator(4,2));
        db_new = db_all_2_beta(:,db); 
        total_circ_test = []; total_comp_test = []; % empty vectors for the end comparison of all k tests
        for k = 1:(db-1)
            db_old = db_all_2_beta(:,k);
            db_equal = NaN((size(db_new,2)-1),1); 
            for i = 1:(size(db_old,2)-1)
                db_shift = circshift(db_new,i);
                db_equal(i,1) = isequal(db_old, db_new);
            end
            circ_test = sum(db_equal) > 0;

            comp = 0;
            for l = 3:9
                for i = 1:(size(db_old,1)-(l-1))
                    for j = 1:(size(db_new, 1)-(l-1))
                        equality = isequal(db_old(i:(i+(l-1)),1),db_new(j:(j+(l-1)),1)); 
                        if equality
                            comp = comp +1*(l/3);
                        end
                    end
                end 
            end
            comp_test = comp > 0;
            total_circ_test = [total_circ_test,circ_test];
            total_comp_test = [total_comp_test,comp_test];
            if sum(total_comp_test) < comp_beta && sum(total_circ_test) == 0
                db_best_so_far = db_new;
                comp_beta = sum(total_comp_test); 
            end
        end
    general = sum(total_circ_test) > 0 || sum(total_comp_test) > 0;
    if counter == 10^6 && general == true
        db_all_2_beta(:,db) = db_best_so_far;
        fails = [fails, all];
    end 
    end
else 
    db_all_2_beta(:,db) = transpose(debruijn_generator(4,2));
end
end
save(['db_all_2_',char(string(all)),'.mat'],'db_all_2_beta');
    %db_144_2 = [db_all_2_beta(:,1);db_all_2_beta(:,2);db_all_2_beta(:,3)];
   % db_total_LvL2(:,all) = db_144_2;
%end
%% Rename all 6 produced orders and save them 
%db_LvL2 = zeros(144,6);
for i = 1:6
    db_all_2_beta = load(['db_all_2_',char(string(i)),'.mat']);
    db_144_2 = [db_all_2_beta.db_all_2_beta(:,1);db_all_2_beta.db_all_2_beta(:,2);db_all_2_beta.db_all_2_beta(:,3)];

    for d = 1:size(db_144_2,1)
        if db_144_2(d,1) == 1
            db_144_2(d,1) = 123;
        elseif db_144_2(d,1) == 2
             db_144_2(d,1) = 321;
        elseif db_144_2(d,1) == 3
             db_144_2(d,1) = 765;
        else 
          db_144_2(d,1) = 567;
        end
    end

    trip_order_db_beta = [];  
    for ij = 1:size(db_144_2,1)
        N = db_144_2(ij,1);
        m = floor(log10(N));  
        D = mod(floor(N ./ 10 .^ (m:-1:0)), 10);
        trip_order_db_beta = [trip_order_db_beta, D]; 
    end 
    trip_144_db = transpose(trip_order_db_beta);
    db_LvL2(:,i) = trip_144_db;
end

%% Save and Export
save('db_all_2.mat','db_LvL2');
writematrix(db_LvL2, 'E:\Arbeit\Programmierung\order_3_0705.txt','Delimiter','tab');