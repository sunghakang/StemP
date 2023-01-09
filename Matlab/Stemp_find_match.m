function  [energy,answer] = Stemp_find_match(vertex)
% edited by Mengyi on 2017/9/30
energy = {};
 s = size(vertex);



for i = 1 : s(1)
    stem_length = cell2mat(vertex(i,2));
    energy(end+1, 1:2) = {i,stem_length};
end

    energy = sortrows(energy,-2);

answer = energy(1);