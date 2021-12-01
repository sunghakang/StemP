function [alignment] = vertexToalignment(rna,vertex)

v = vertex;
    alignment = '.';

    for i = 1 : length(rna)-1
        alignment = strcat(alignment,'.');
    end

    for i = 1 : length(v)
                    
        % imod3 =1
        if rem(i,3) == 1
        stack= cell2mat(v(i));
        for j = 1 : length(stack)
            if rem(j,2) == 1
                alignment(stack(j)) = '(';
            elseif rem(j,2) == 0
                alignment(stack(j)) = ')';
            end
        end
        
        % i is even
        else if rem(i,3) == 2
        stack= cell2mat(v(i));
        for j = 1 : length(stack)
            if rem(j,2) == 1
                alignment(stack(j)) = '[';
            elseif rem(j,2) == 0
                alignment(stack(j)) = ']';
            end
        end 
        
            else
            stack= cell2mat(v(i));
        for j = 1 : length(stack)
            if rem(j,2) == 1
                alignment(stack(j)) = '{';
            elseif rem(j,2) == 0
                alignment(stack(j)) = '}';
            end
        end     
            
            end
        end
        
    end    
end

