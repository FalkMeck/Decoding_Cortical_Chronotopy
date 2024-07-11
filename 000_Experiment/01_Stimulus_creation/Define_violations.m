%% define violations
clear; cd('F:\Arbeit\Programmierung\NewSeq');
levels = ['LvL1';'LvL2';'LvL3';'True'];
rpl_num = 0:143;
vio_total = cell(2,4);
for i = 1:size(levels,1)
    fii = fopen([cd(), '\db_input_','ins','_',levels(i,:),'.txt']);
    d_ins = textscan(fii, '%f%f%f%f%f%f', 'Delimiter', 'Tab');
    fii = fclose(fii);
    fir = fopen([cd(), '\db_input_','rpl','_',levels(i,:),'.txt']);
    d_rpl = textscan(fir, '%f%f%f%f%f%f', 'Delimiter', 'Tab');
    fir = fclose(fir);
    vio_level_ins = cell(1,6);
    vio_level_rpl = cell(1,6);
    for j = 1:size(d_ins,2)
    % INSERTION
    % only counts non violation number, vilation have double number
    % writes out violation  number and devides them bie repspective
    % seqences lengths (1,3,9,36) -> thus knowing if violation comes after
    % a sequence breake or breaking within a sequence
        ins_beta = [d_ins{1,j}(1:end) zeros(size(d_ins{1,j},1),(2+i))];
        count_non_vio = 0;
        ins_vio = ins_beta(:,1) == 4; 
        for k = 1:size(ins_beta,1) % only count when, it is no vio, so vios are the same number twice
            if ins_beta(k,1) ~= 4
                count_non_vio = count_non_vio +1;
            end
            if ins_vio(k)
                ins_beta(k,2) = count_non_vio; % only write aour these double numbers --> because insertet after 
            end
        end
        while_count = 0; 
        while while_count <= i
            while_count = while_count +1;
            if while_count == 1
                ins_beta(ins_vio,while_count+2) = 1;
            elseif while_count == 2
                % always: take rows in which the divintion works without
                % rest -> rows between natural sequence breaks
                ins_beta(ins_beta(:,2)/3 == floor(ins_beta(:,2)/3)& ins_beta(:,2) ~= 0, while_count+2) = 1;
            elseif while_count == 3
                ins_beta(ins_beta(:,2)/9 == floor(ins_beta(:,2)/9)& ins_beta(:,2) ~= 0, while_count+2) = 1;
            elseif while_count == 4
                ins_beta(ins_beta(:,2)/36 == floor(ins_beta(:,2)/36)& ins_beta(:,2) ~= 0, while_count+2) = 1;
            end
        end

        ins_beta(:,i+3) = sum(ins_beta(:,3:(i+2)),2);
        vio_level_ins{1,j}=ins_beta;

   % REPLACMENT
   % position - 1/ sequences lengths --> tells if it is at border  
        rpl_beta = [d_rpl{1,j}(1:end) zeros(size(d_rpl{1,j},1),(2+i))];
        rpl_beta(rpl_beta(:,1)==4,2) = rpl_num(rpl_beta(:,1)==4);
        
        rpl_vio = rpl_beta(:,1) == 4; 
        rpl_beta(rpl_vio,2) = rpl_num(rpl_vio); % rpl num = 0:143 --> here come the -1

        while_count = 0; 
        while while_count <= i
            while_count = while_count +1;
            if while_count == 1
                rpl_beta(rpl_vio,while_count+2) = 1;
            elseif while_count == 2
                rpl_beta(rpl_beta(:,2)/3 == floor(rpl_beta(:,2)/3)& rpl_beta(:,2) ~= 0, while_count+2) = 1;
            elseif while_count == 3
                rpl_beta(rpl_beta(:,2)/9 == floor(rpl_beta(:,2)/9)& rpl_beta(:,2) ~= 0, while_count+2) = 1;
            elseif while_count == 4
                rpl_beta(rpl_beta(:,2)/36 == floor(rpl_beta(:,2)/36)& rpl_beta(:,2) ~= 0, while_count+2) = 1;
            end
        end

        rpl_beta(:,i+3) = sum(rpl_beta(:,3:(i+2)),2);
        vio_level_rpl{1,j}=rpl_beta;
    end
    vio_total{1,i} = vio_level_ins;
    vio_total{2,i} = vio_level_rpl;
end

 for i = 1:2
     for j = 1:4
         for k = 1:6
             disp(sum(vio_total{i,j}{1,k}(:,end) < 3));
         end
    end
end


for ie = 1:4
    txt_matrix_rpl = zeros(144,6);
    txt_matrix_ins = zeros(162,6); 
    for je = 1:6
        txt_matrix_rpl(:,je) = vio_total{2,ie}{1,je}(:,ie+3);
        txt_matrix_ins(:,je) = vio_total{1,ie}{1,je}(:,ie+3);
    end
    writematrix(txt_matrix_rpl,['violations_rpl_', levels(ie,:), '.txt'],'Delimiter','tab');
    writematrix(txt_matrix_ins,['violations_ins_', levels(ie,:), '.txt'],'Delimiter','tab');
end