function [vertex] = find_vertex_InARow0310(rna,number1,number2,s_min,s_max)

stack_length = 2;
vertex = {};


for i = 1 : length(rna)
    for j = (i+2) : length(rna)
        starting = i;
        ending = j;
        addvertex = [];
        
        while (base_pairing(rna(starting),rna(ending)) ...
                || wobble_pairing(rna(starting),rna(ending)))...
                && (starting < ending) && length(addvertex) < number2*2
            x = [starting,ending];
            addvertex = [addvertex, x];
            starting = starting + 1;
            ending = ending - 1;
            if length(addvertex)>= number1 * 2
                distance = addvertex(2) - addvertex(1);
                stem_length = length(addvertex)/2;
                ratio = distance / stem_length;
                if  ratio >= s_min && ratio <=s_max 
                    vertex(end+1, 1:4) = {addvertex, stem_length, distance, ratio};
                end
                
            end
        end
        
    end
    
end


end

