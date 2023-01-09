function [edge] = StemP_find_edge_5s(vertex, index1, index2)
% Find edges between three domains
edge = {};

[w, ~] = size(vertex);

% Not Consider pseudoknots

for i = 1:index1
    a = cell2mat(vertex(i));
    for j = index1+1 : w
        b = cell2mat(vertex(j));
        
        if a(2) < b(1) ||...
                (a(end-1) < b(1) && b(2) < a(end))
            s = [i,j];
            edge{end+1, 1} = s;
            
            
        else if b(2) < a(1) ||...
                    (b(end-1) < a(1) && a(2) < b(end))
                s = [i,j];
                edge{end+1,1} = s;
            end
        end
        
    end
end

for i = index1+1:index2
    a = cell2mat(vertex(i));
    for j = index2+1 : w
        b = cell2mat(vertex(j));
        
        if a(2) < b(1) ||...
                (a(end-1) < b(1) && b(2) < a(end))
            s = [i,j];
            edge{end+1, 1} = s;
            
            
        else if b(2) < a(1) ||...
                    (b(end-1) < a(1) && a(2) < b(end))
                s = [i,j];
                edge{end+1,1} = s;
            end
        end
        
    end
end

end



