function [p] = datacheckbyvertex2(v1,v2,l,type)
% Author: Mengyi
% First update: 12/21/2020
% Last update: 12/21/2020
% compute accuracy between true strucutre and predicted structure
% v1 : True matching base pairs
% v2 : prediced matching base pairs
% TP : correctly matched
% FP : matched --> incorrect matched / alone
% TN : correctly unmatched
% FN : single  --> incorrect matched

same_alignment = 0;
TP = 0; TN = 0;
FP = 0; FN = 0;
% l is the length of rna seq

%% the following part is to delete the nonpairing line from v1

index=[];
for k = 1 : length(v1)
    if v1(k,1)==0 || v1(k,2)==0
        index=[index,k];
    end
end
v1(index,:)=[];

index=[];
for k = 1 : length(v2)
    if v2(k,1)==0 || v2(k,2)==0
        index=[index,k];
    end
end
v2(index,:)=[];


%%
for i = 1: l
    [r1,c1] = find(v1==i);
    if ~isempty(r1)
        r1 = r1 (1,:); %assume there are duplicated vertex from the collecing data
    end
    if ~isempty(c1)
        c1 = c1 (1,:);
    end
    
    
    [r2,c2] = find(v2==i);
    
    if ~isempty(r2)
        r2 = r2 (1,:);
    end
    if ~isempty(c2)
        c2 = c2 (1,:);
    end

    
    if isempty(r1) == isempty(r2)
        if isempty(r1) == 1
            same_alignment = same_alignment + 1;
            TN = TN +1;
        elseif v1(r1,3-c1)==v2(r2,3-c2)
            same_alignment = same_alignment + 1;
            TP = TP + 1;
        else
            FP = FP + 1;
        end
    elseif isempty(r1) == 1
        FN = FN + 1;
    else
        FP = FP +1;
    end
end

acc = same_alignment / l;
sens = TP / (TP + FN);
ppv = TP / (TP + FP);
mcc = sqrt(sens * ppv);
if type == "sens"
    p = sens;
elseif type == "ppv"
    p = ppv;
elseif type == "mcc"
    p = mcc;
elseif type == "acc"
    p = acc;
end

end


