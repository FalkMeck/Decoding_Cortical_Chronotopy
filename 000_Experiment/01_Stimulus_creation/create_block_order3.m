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