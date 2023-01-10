function [vertex] = find_vertex_nPlusn0309(rna,s_min,s_max,l1, l2, gap1, gap2)

vertex = {};
total_length = l1+l2;

for i = 1 : length(rna)
    for j = (i+2) : length(rna)
        starting = i;
        ending = j;
        addvertex = [];
        k=0;
        while (base_pairing(rna(starting),rna(ending)) ...
                || wobble_pairing(rna(starting),rna(ending)))...
                && (starting < ending) && length(addvertex) < total_length *2
            
            k = k+1;
            x = [starting,ending];
            addvertex = [addvertex, x];
            
            if k ~= l1
            starting = starting + 1;
            ending = ending - 1;
            else 
                for k1 = 1:gap1
                    for k2 = 1:gap2 

                        if  ending-k2>=1 &&...
                            (base_pairing(rna(starting+k1),rna(ending-k2)) ...
                                || wobble_pairing(rna(starting+k1),rna(ending-k2)))...
                                
                            continue
                        end
                    end
                end
                starting = starting +1+gap1;
                ending = ending -1-gap2;
            end
            if starting > length(rna) || ending <= starting 
                break
            end
            
        end

        
        
        if length(addvertex) == total_length * 2
            distance = addvertex(2) - addvertex(1);
          stem_length = length(addvertex)/2;
            ratio = distance / stem_length;
            if s_min == 0 && s_max == 0
                vertex(end+1, 1:4) = {addvertex, stem_length, distance, ratio};
            else if ratio >= s_min && ratio <=s_max
                    vertex(end+1, 1:4) = {addvertex, stem_length, distance, ratio};
                end
            end
        end
    end
       
end

end

