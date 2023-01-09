function pred_alignments = stemP_pred_tRNA(rna,L,d1,d2,s2,s3)
pred_alignments = [];
vertex = stemP_find_vertex_tRNA(rna,L,d1,d2,s2,s3);
disp(['Number of vertex found  : ', num2str(size(vertex,1))])
if isempty(vertex)
    disp('No vertex found');
    return
end
edge = find_edge(vertex,0);
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

