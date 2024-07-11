% W. Owen Brimijoin (2021). de Bruijn sequence generator (https://www.mathworks.com/matlabcentral/fileexchange/28015-de-bruijn-sequence-generator), MATLAB Central File Exchange. Retrieved May 5, 2021.
db_total_LvL3 = zeros(16,6);
fails = []; 
warning('off','all');
for db = 1:size(db_total_LvL3,2)
if db > 1 
    general = true;
    %counter = 0;
    %comp_beta = 1000; 
    %db_best_so_far = 0;
    while general %&& counter <= 10^6
        %counter = counter +1; 
        db_total_LvL3(:,db) = transpose(debruijn_generator(4,2));
        db_new = db_total_LvL3(:,db); 
        total_circ_test = []; % total_comp_test = [];
        for k = 1:(db-1)
            db_old = db_total_LvL3(:,k);
            db_equal = NaN((size(db_new,2)-1),1); 
            for i = 1:(size(db_old,2)-1)
                db_shift = circshift(db_new,i);
                db_equal(i,1) = isequal(db_old, db_new);
            end
            circ_test = sum(db_equal) > 0;

%             comp = 0; % the within check is not necessary because 1
%             DB-Srequence is already 144 
%             for l = 3:9
%                 for i = 1:(size(db_old,1)-(l-1))
%                     for j = 1:(size(db_new, 1)-(l-1))
%                         equality = isequal(db_old(i:(i+(l-1)),1),db_new(j:(j+(l-1)),1)); 
%                         if equality
%                             comp = comp +1*(l/3);
%                         end
%                     end
%                 end 
%             end
%             comp_test = comp > 0;
            total_circ_test = [total_circ_test,circ_test]; % only circ test
%             total_comp_test = [total_comp_test,comp_test];
%             if sum(total_comp_test) < comp_beta && sum(total_circ_test) == 0
%                 db_best_so_far = db_new;
%                 comp_beta = sum(total_comp_test); 
%             end
        end
    general = sum(total_circ_test) > 0 ;% || sum(total_comp_test) > 0;
%     if counter == 10^6 && general == true
%         db_all_2(:,db) = db_best_so_far;
%         fails = [fails, all];
%     end 
    end
else 
    db_total_LvL3(:,db) = transpose(debruijn_generator(4,2));
end
end
save('db_all_3_beta.mat','db_total_LvL3');
%% Rename and save
for w = 1:size(db_total_LvL3,2)
    for d = 1:size(db_total_LvL3,1)
        if db_total_LvL3(d,w) == 1
            db_total_LvL3(d,w) = 123123123;
        elseif db_total_LvL3(d,w) == 2
             db_total_LvL3(d,w) = 321321321;
        elseif db_total_LvL3(d,w) == 3
             db_total_LvL3(d,w) = 765765765;
        else 
          db_total_LvL3(d,w) = 567567567;
        end
    end
end  

none_order_db = zeros(144,6);
for w = 1: size(db_total_LvL3,2)
    none_order_db_beta = [];
    for ij = 1:size(db_total_LvL3,1)
        N = db_total_LvL3(ij,w);
        m = floor(log10(N));  
        D = mod(floor(N ./ 10 .^ (m:-1:0)), 10);
        none_order_db_beta = [none_order_db_beta, D]; 
    end 
    none_order_db(:,w) = transpose(none_order_db_beta); 
end

%% Save and Export
save('db_all_3.mat','none_order_db');
writematrix(none_order_db, 'E:\Arbeit\Programmierung\order_9_0605.txt','Delimiter','tab');