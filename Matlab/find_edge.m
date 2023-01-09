function [edge] = find_edge_only_first_layer(vertex, pseudoknots)
% 
% last update: 13/01/2020


edge = {};

[w, ~] = size(vertex);

% Not Consider pseudoknots
if pseudoknots == 0 
    for i = 1:w
        a = cell2mat(vertex(i,1));
        for j = i+1 : w
            b = cell2mat(vertex(j,1));
            if a(2) < b(1) || b(2) < a(1) ||(a(end-1)<b(1) && b(2)<a(end))
                s = [i,j];
                edge{end+1, 1} = s;
            end
        end
    end

% Consider pseudoknots
else

    for i = 1:w
        a = cell2mat(vertex(i));
        for j = i+1 : w
            b = cell2mat(vertex(j));
            if a(2) < b(1) || ...
                    (a(end-1) < b(1) && b(2) < a(end)) ||...
                    (a(end-1) < b(1) &&...
                    b(end-1) < a(end) && a(2) < b(end))
                s = [i,j];
                edge{end+1, 1} = s;
            end
        end
    end
    
end


end
     


