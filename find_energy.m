function [energy,answer] = find_energy(cliques,vertex)
% 
% Usage
%   [energy,answer] = find_energy(cliques,vertex)
% input: maximal_cliques -
%        vertex - 
% output: energy - 
%         answer - 
%

energy = {};

for i = 1 : length(cliques)
    clique_set = cell2mat(cliques(i));
    sum_length = 0;
    for j = 1 : length(clique_set)
        length_stack = length(cell2mat(vertex(clique_set(j)))) / 2;
        sum_length = sum_length + length_stack;
    end
    energy(end+1, 1:2) = {cliques{i,1},sum_length};
end

energy = sortrows(energy,-2);
answer = energy(1);
