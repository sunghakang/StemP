function [pred_number,pred_number2,predpercentage,predpercentage2,alignment_all] = StemP_pred_accuracy(rna,energy, vertex,TrueVertex)

leng_rna=length(rna);
highest_pred = 0;
pred_number = 0;

[length_energy,~] = size(energy);

for i = 1 : length_energy
    PredVertexcell = CliqueToVertex(vertex,energy(i));
    PredVertex = VertexToMatrix(PredVertexcell);
    pred = datacheckbyvertex2(TrueVertex,PredVertex,leng_rna,"mcc");
    energy(i, 3) = {pred};
    
    if highest_pred < pred
        pred_number = i;
        pred_number2 = find_rank(energy,i);
        highest_pred = pred;
        phe_vertex = PredVertex;
    end
    
end

predpercentage = datacheckbyvertex2(TrueVertex,phe_vertex,leng_rna,"mcc");
predpercentage2 = energy{1,3};
alignment_all = pred_alignment(rna,energy,vertex,pred_number);

end

