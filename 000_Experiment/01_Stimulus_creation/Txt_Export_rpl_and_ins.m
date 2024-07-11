clear; cd('F:\Arbeit\Programmierung');
LvL_name = ['LvL1';'LvL2';'LvL3';'True'];

for i = 1:4
    %load(['Change_positions_', LvL_name(i,:), '_2304.mat']); 
    %load(['Change_positions_', LvL_name(i,:), '_1405.mat']); 
    %load(['Change_positions_', LvL_name(i,:), '_1307.mat']);
    load(['Change_positions_', LvL_name(i,:), '_2707.mat']);
    txt_matrix_rpl = zeros(144,6);
    txt_matrix_ins = zeros(162,6); 
    for j = [1,3,5]
        txt_matrix_rpl(:,j) = change_positions{3,(j/2+0.5)}(1:144,1);
        txt_matrix_ins(:,j) = change_positions{5,(j/2+0.5)}(1:162,1);
    end
    for j = [2,4,6]
        txt_matrix_rpl(:,j) = change_positions{3,(j/2)}(145:288,1);
        txt_matrix_ins(:,j) = change_positions{5,(j/2)}(163:324,1);
    end
    writematrix(txt_matrix_rpl,['F:\Arbeit\Programmierung\NewSeq\db_input_rpl_', LvL_name(i,:), '.txt'],'Delimiter','tab');
    writematrix(txt_matrix_ins,['F:\Arbeit\Programmierung\NewSeq\db_input_ins_', LvL_name(i,:), '.txt'],'Delimiter','tab');
end
