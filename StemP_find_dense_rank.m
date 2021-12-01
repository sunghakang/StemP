% Created by Mengyi 2/18/2021
function [rank] = StemP_find_dense_rank(energy,i)
rank = i;
    while (i > 1) 
        if (cell2mat(energy(i,2)) == cell2mat(energy(i-1,2)))
            rank = rank -1;
        end
        i = i-1;
    end
end

