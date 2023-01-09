function pred_alignments = stemP_pred_5srRNA(rna)
pred_alignments = [];
[vertex,index1,index2] = stemP_find_vertex_rRNA5S_Archaeal(rna);
disp(['Number of vertex found  : ', num2str(size(vertex,1))])
if isempty(vertex)
    disp('No vertex found');
    return
end
edge = StemP_find_edge_5s(vertex,index1,index2);
disp(['Number of edges found   : ', num2str(size(edge,1))])
if isempty(edge)
    disp('No vertex found');
end
if ~isempty(edge)% case 1 : when there is no edge
    [~,adjacency_matrix] = find_graph(edge);
    % 2. find cliques in the graph
    MC = Bron_Kerbosch(adjacency_matrix,2);
    cliques = find_cliques(MC);
    disp(['Number of cliques found : ', num2str(size(cliques,1))])
    
    % 3. find maximum matching from cliques
    [energy,~]= find_energy(cliques,vertex);
else % case 2
    [energy,~] = find_maxmatching(vertex);
end

pred_alignments = get_top_alignments(rna,energy,vertex);

end
