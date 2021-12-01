function vertex = find_vertex_PDB(rna,stack_length, wobble, uu, ratio_min, ratio_max)

% 
% Usage
% vertex = find_vertex(rna, stack_length, wobble, ratio_min, ratio_max)
% input: rna - rna seuqence 
%        stack_length - length of base pairing. >= stack_length
%        wobble - 0 for Watson-Clark base pairing and 
%                 1 for both W-C and wobble pairing
%        uu     - 0 for U-U base pairing and 
%                 1 for both W-C and wobble pairing
%        ratio_min - staring ratio level
%                    0 default no limit for minimum of ratio
%        ratio_max - ending ratio level
%                    0 default no limit for maximum of ratio
% output: vertex - information of constructed verices 
%                  column 1 : position of rna squence 
%                  column 2 : stack length
%                  column 3 : ending - starting sequence position 
%                  column 4 : ratio
%

vertex = {};
vertex1 = {};

% 0. find all possible vertex whose length is greater than stack_length 
% watson clarke pairing and wobble pairing(G-U)

if (wobble == 1) && (uu == 0)
for i = 1 : length(rna)
    for j = (i+3) : length(rna)
        starting = i;
        ending = j;
        addvertex = [];
        while (base_pairing(rna(starting),rna(ending)) || wobble_pairing(rna(starting),rna(ending)))...
                && (starting < ending)
            x = [starting,ending];
            addvertex = [addvertex, x];
            starting = starting + 1;
            ending = ending - 1;
        end

        if length(addvertex) >= stack_length * 2
            distance = addvertex(2) - addvertex(1);
            stem_length = length(addvertex)/2;
            ratio = distance / stem_length;
            vertex1(end+1, 1:4) = {addvertex, stem_length, distance, ratio};

        end
    end
end

else  if (wobble == 0) && (uu == 0) % without wobble_paring
for i = 1 : length(rna)
    for j = (i+3) : length(rna)
        starting = i;
        ending = j;
        addvertex = [];
        while (base_pairing(rna(starting),rna(ending)) && (starting < ending))
            x = [starting,ending];
            addvertex = [addvertex, x];
            starting = starting + 1;
            ending = ending - 1;
        end

        if length(addvertex) >= stack_length * 2
            distance = addvertex(2) - addvertex(1);
            stem_length = length(addvertex)/2;
            ratio = distance / stem_length;
            vertex1(end+1, 1:4) = {addvertex, stem_length, distance, ratio};

        end
    end
end

    


else  if (wobble == 1) && (uu == 1) % without wobble_paring
for i = 1 : length(rna)
    for j = (i+3) : length(rna)
        starting = i;
        ending = j;
        addvertex = [];
        while (base_pairing(rna(starting),rna(ending))  || wobble_pairing(rna(starting),rna(ending))...
                 || uu_pairing(rna(starting),rna(ending)))...
                 && (starting < ending)
             
            x = [starting,ending];
            addvertex = [addvertex, x];
            starting = starting + 1;
            ending = ending - 1;
        end

        if length(addvertex) >= stack_length * 2
            distance = addvertex(2) - addvertex(1);
            stem_length = length(addvertex)/2;
            ratio = distance / stem_length;
            vertex1(end+1, 1:4) = {addvertex, stem_length, distance, ratio};

        end
    end
end



else
for i = 1 : length(rna)
    for j = (i+3) : length(rna)
        starting = i;
        ending = j;
        addvertex = [];

        while (base_pairing(rna(starting),rna(ending))==1  || uu_pairing(rna(starting),rna(ending))==1 )...
                 && (starting < ending)
            x = [starting,ending];
            addvertex = [addvertex, x];
            starting = starting + 1;
            ending = ending - 1;
        end

        if length(addvertex) >= stack_length * 2
            distance = addvertex(2) - addvertex(1);
            stem_length = length(addvertex)/2;
            ratio = distance / stem_length;
            vertex1(end+1, 1:4) = {addvertex, stem_length, distance, ratio};
        end
    end
end


end

end

end


[w, ~] = size(vertex1);

%information of vertices
if ratio_min == 0 && ratio_max == 0
    for i = 1 : w
        a = cell2mat(vertex1(i));
        distance = a(2) - a(1);
        stack_length = length(a)/2;
        ratio = distance / stack_length;
        vertex(end+1, 1:5) = {a, stack_length, distance, ratio,i};
    end
    
else
    for i = 1 : w
        a = cell2mat(vertex1(i));
        distance = a(2) - a(1);
        stack_length = length(a)/2;
        ratio = distance / stack_length;   
        if ratio <= ratio_max && ratio > ratio_min
            vertex(end+1, 1:5) = {a, stack_length, distance, ratio,i};
        end
    end
end

