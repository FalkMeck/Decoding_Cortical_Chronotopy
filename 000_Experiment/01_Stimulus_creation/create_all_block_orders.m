%% create n possible diffrent block orders with variying starting numbers
clear; cd("F:\Arbeit\Programmierung");
n = 50; % 50 diffrent orders
block_orders = cell(1,4); % empty cell array
for s = 3
    block_order_s = zeros(24,n); % 24 blocks times n 
    for k = 1:n
        if k == 1
            % block_order_s_k = create_block_order(s);% create block order with own create_block_order function
            block_order_s_k = create_block_order3(s);% create block order with own create_block_order function
            block_order_s(1:end,k) = block_order_s_k;
        else 
           success_k = false;
           success_test_k = zeros(1,k-1);
           while ~success_k
            %block_order_s_k = create_block_order(s); % create block order with own create_block_order function
            block_order_s_k = create_block_order3(s);
            for e = 1:(k-1) % compare if two orders are the same
                if isequal(block_order_s_k,block_order_s(1:end,e))
                    success_test_k(1,e) = 1;
                end
            end
            success_k = sum(success_test_k) == 0; % only success, if its not equal to anything
           end
           block_order_s(1:end,k) = block_order_s_k; % save Block order
        end
    end
    block_orders{1,s} = block_order_s; % save all block otder for this star point in cell array
end

%% export all as Txt files
stimuli_folder = 'F:\Arbeit\Programmierung\NewBlock';
n = 50;
for s = 3 % svae each block order under \block_order_s_nn.txt
    for k = 1:n
        filename = [stimuli_folder, '\block_order_', char(string(s)),'_',sprintf('%02d', k),'.txt'];
        fileID = fopen(filename,'w');
        fprintf(fileID,'%d\n',block_orders{1,s}(1:end,k));
        fclose(fileID);
    end
end

%% block_order 3
function db_final = create_block_order3(starting_value)
    warning('off');
    success_total = false;
    while ~success_total
        db_total = zeros(1,32);

        success1 = false;
        while ~ success1
          db_beta = debruijn_generator(4,2);
          success1 = db_beta(1) == starting_value && ...
              db_beta(end) ~= starting_value && ...
              db_beta(5) == starting_value && ...
              length(unique(db_beta(2:4))) == 3 && ~ismember(starting_value, db_beta(2:4));
        end
        db_total(1,1:16) = db_beta;

        success2 = false;
        while ~success2
          db_beta = debruijn_generator(4,2);
          success2 = db_beta(1) == starting_value && ....
              db_beta(end) ~= starting_value && ...
              ~isequal(db_beta, db_total(1,1:16));
        end 
        db_total(1, 17:32) = db_beta;
        db_total_short = db_total([true diff(db_total)~=0]);

        success_total = length(unique(db_total_short(5:9))) == 4 && ...
                        length(unique(db_total_short(10:14))) == 4 && ...
                        length(unique(db_total_short(15:19))) == 4 && ...
                        length(unique(db_total_short(20:24))) == 4;
    end
    db_final = db_total_short;
    return
end

%% block_order 2
% function db_final = create_block_order2(starting_value)
%     warning('off');
%     db_total = zeros(1,32);
%     
%     success1 = false;
%     while ~ success1
%       db_beta = debruijn_generator(4,2);
%       success1 = db_beta(1) == starting_value && ...
%           db_beta(end) ~= starting_value && ...
%           db_beta(5) == starting_value && ...
%           length(unique(db_beta(2:4))) == 3 && ~ismember(starting_value, db_beta(2:4));
%     end
%     db_total(1,1:16) = db_beta;
%     
%     success2 = false;
%     while ~success2
%       db_beta = debruijn_generator(4,2);
%       success2 = db_beta(1) == starting_value && ....
%           db_beta(end) ~= starting_value && ...
%           ~isequal(db_beta, db_total(1,1:16));
%     end 
%     db_total(1, 17:32) = db_beta;
% 
%     db_final = db_total([true diff(db_total)~=0]);
%     return
% end

%% block_order 1
% function db_final = create_block_order(starting_value)
%     warning('off');
%     db_total = zeros(1,32);
%     for i = 1:2
%         success = false;
%         while ~ success
%             db_beta = debruijn_generator(4,2);
%             success = db_beta(1) == starting_value && db_beta(end) ~= starting_value;
%             if i == 2 && isequal(db_total(1:16),db_beta)
%                 success = false;
%             end
%         end
%         db_total(1, ((i-1)*16+1):((i-1)*16+16)) = db_beta;
%     end
%     db_final = db_total([true diff(db_total)~=0]);
%     return
% end