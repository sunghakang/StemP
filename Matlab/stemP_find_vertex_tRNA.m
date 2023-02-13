function vertex = stemP_find_vertex_tRNA(rna,stack_length,d1,d2,s2,s3)
% Author: Kumbit Hwang
% update: 01/15/2016, 11/30/2017, 2021
%rna = 'GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA'
%stack_length = 3;
vertex = {};
vertex1 = {};
% s1 is the ratio bound for the acceptor stem
% usually its 3.41
% s2 = 3
% s3 = 4.6 --> 4.7
% d1 = 12
% d2 = 17 --> 4.8


% 0. find all possible vertex whose length is greater than stack_length 
% watson clarke pairing and wooble pairing(G-U)

for i = 1 : length(rna)
    for j = (i+3) : length(rna)
        starting = i;
        ending = j;
        addvertex = [];
        
            
        while ((base_pairing(rna(starting),rna(ending)) || wobble_pairing(rna(starting),rna(ending)))...
                && (starting < ending)) == 1
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

    
% 
% 
%information of vertices
for i = 1 : size(vertex1,1)
    
    a = cell2mat(vertex1(i));
    distance = a(2) - a(1);
    stem_length = length(a)/2;
    ratio = distance / stem_length;
    
    if distance > length(rna)/2
        % open-ended modified
%         ratio = distance / ( stem_length + 22) +1;
%         
%         while ratio < s1
%         stem_length = stem_length - 1;
%         ratio = distance / ( stem_length + 22) +1;
%         end   
%         a = a(1:2*stem_length); 
        
        distance = -distance +length(rna)+2*stem_length-2;
        ratio = distance/stem_length;
        a = a(1:2*stem_length);
        if ratio <3 
            vertex(end+1, 1:4) = {a, stem_length, distance, ratio};
        end
        
    else
        
        while ratio < 3
            stem_length = stem_length - 1;
            ratio = distance ./ stem_length;
        end
        a = a(1: 2*stem_length);
        if (distance >= d1 && distance <= d2 ) ...
            && ratio >= s2 && ratio <= s3 
            vertex(end+1, 1:4) = {a, stem_length, distance, ratio};
        end
    
    end
    
%     if (distance >= d1 && distance <= d2 ) ...
%             && ratio >= s2 && ratio <= s3 && stem_length >= stack_length
%             vertex(end+1, 1:4) = {a, stem_length, distance, ratio};
%     end
    
end

end











