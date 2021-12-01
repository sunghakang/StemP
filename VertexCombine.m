function [VertexNew] = VertexCombine(VertexLayerOne,VertexLayerTwo,s_min,s_max)
[r1,c1] = size(VertexLayerOne);
[r2,c2] = size(VertexLayerTwo);
VertexNew = {};
for i = 1:r1
    v1 = cell2mat(VertexLayerOne(i,1));
    for j = 1:r2
        v2= cell2mat(VertexLayerTwo(j,1));
        if v1(end-1) < v2(1) && v1(end) >v2(2)
            v_new = [v1,v2];
            distance = v_new(2) - v_new(1);
            stem_length = length(v_new)/2;
            ratio = distance / stem_length;
            if (s_min ==0 && s_max ==0 ) || ( ratio <= s_max && ratio >= s_min)
            VertexNew(end+1, 1:4) = {v_new, stem_length, distance, ratio}; 
        end
    end
end



end

