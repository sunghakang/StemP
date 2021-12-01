function [alignment_all] = pred_alignment(rna,energy,vertex,pred_number)

alignment_all = {};

for k = 1:pred_number
    
    vertex_pred={};

    clique = cell2mat(energy(k,1));

 
    for i = 1:length(clique)
        vertex_pred(end+1,1) = vertex(clique(i),1);
    end

    v=vertex_pred;
     
        

    alignment = '.';

    for i = 1 : length(rna)-1
        alignment = strcat(alignment,'.');
    end

    for i = 1 : size(v,1)
                    
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
        
        
        
        % imod3 = 2
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
    alignment_all(end+1,1)={alignment};   
end

end



