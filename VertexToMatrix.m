% Convert data type of vertex
% First update 10/16/2019 Mengyi
function [vertexmatrix] = VertexToMatrix(vertex)

addvertexcol1 = [];
addvertexcol2 = [];

for i = 1:length(vertex)
    nun = vertex(i);
    nun = cell2mat(nun);
    for j = 1:2:length(nun)-1
        addvertexcol1(end+1) = nun(j);
    end
    for j = 2:2:length(nun)
        addvertexcol2(end+1) = nun(j);
    end
end

vertexmatrix = [ addvertexcol1; addvertexcol2 ];
vertexmatrix = vertexmatrix';   
   

end

