function [v] = CliqueToVertex(vertex,clique)
% Author: Mengyi
% First update: 10/16/2019
% Last update: 10/16/2019
v={};
clique=cell2mat(clique);
[~,l]= size(clique);

for i = 1 : l
    v(end+1,1)=vertex(clique(i),1);
end

end

