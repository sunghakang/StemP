function [i] = find_rank(energy,i)
    while (i > 1) && (cell2mat(energy(i,2)) == cell2mat(energy(i-1,2)))
     i = i - 1;
    end

end

