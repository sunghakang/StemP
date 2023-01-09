function [rna,vertexmatrix,TrueVertex,vertexoutput,p,w,uu] = readpdbseq(inputArg1)
% Author: Mengyi
% First update: 10/16/2019
% Last update: 10/16/2019
%clear; close all; clc;
%%
data = [];

%%
 [~,txt,raw] = xlsread('PDB_data2.xlsx',inputArg1);

 
 initialindex =  raw(1,3);
 initialindex = cell2mat(initialindex);
 
 k = raw(1,4);
 k = cell2mat(k);
 
 p = raw(2,3);
 p = cell2mat(p);
 
 w = raw(2,4);
 w = cell2mat(w);
 
 uu = raw(2,5);
 uu = cell2mat(uu);
 
 rna = txt(1,5);
 rna = str2mat(rna);

 firstColumn = txt(:,1);
 secondColumn = txt(:,2);
 
 
 FirstColumn = cell2mat(firstColumn);
 SecondColumn = cell2mat(secondColumn);
 
 %%  
 newFirst = [];
 newFirst_row = [];
 newSecond = [];
 newSecond_row = [];
 l = length(cell2mat(firstColumn(1,1)));
%%  

for i = 1:length(firstColumn)
    newFirst_row = [];
    newSecond_row = [];
    for j = 1:k

        newFirst_row = [FirstColumn(i,l-j+1) newFirst_row];
        newSecond_row = [SecondColumn(i,l-j+1) newSecond_row];
    end
     newFirst = [newFirst;newFirst_row];
     newSecond = [newSecond; newSecond_row];
     
end


%% 
pairingdata = [ str2num(newFirst) str2num(newSecond)];
pairingdata = pairingdata - initialindex;
v = pairingdata;
TrueVertex = v;

%% 
vertexmatrix = {};
vertex1 = {};

stack_length = 1;

i = 1;
while i <= length(rna)-3
    i_0 = i;
    for j = (i+3) : length(rna)
        starting = i;
        ending = j;
        addvertex = [];
        while (find_basepairingture(v,starting,ending) && (starting < ending))
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
            
            vv = vertex1(end,1);
            vv = cell2mat(vv);
            vv = vv(end-1);
            i = vv + 1;
        end
        
    end
    
    if i == i_0
        i = i + 1;
    end
end
    
vertexoutput = sortrows(vertex1,4);
vertexmatrix = vertexoutput(:,1);


%writetable(rna,vertexoutput, 'inputArg1.txt')

end

