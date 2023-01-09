function [vertex] = find_vertex_OnePlusSixPlusOne(rna,s_min,s_max)

vertex = {};
total_length = 8;
gapindex1 = 1;
gapindex2 = 7;

for i = 1 : length(rna)
    for j = (i+2) : length(rna)
        starting = i;
        ending = j;
        addvertex = [];
        k=0;
        while ending>0 && (base_pairing(rna(starting),rna(ending)) ...
                || wobble_pairing(rna(starting),rna(ending)))...
                && (starting < ending) && length(addvertex) < total_length *2
            
            k = k+1;
            x = [starting,ending];
            addvertex = [addvertex, x];
            
            if k == gapindex1
                starting = starting + 2;
                ending = ending - 2;
            else if k == gapindex2
                    starting = starting +3;
                    ending = ending -1;
                else
                    starting = starting +1;
                    ending = ending -1;
                end
            end
            
        end
        
        
        
        if length(addvertex) == total_length * 2
            distance = addvertex(2) - addvertex(1);
            stem_length = length(addvertex)/2; % modifed
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

