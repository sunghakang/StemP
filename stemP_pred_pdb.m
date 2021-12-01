function Pred_alignments = stemP_pred_pdb(rna,L,w,uu,LS1,LS2, pseudoknot)
% Predict PDB sequence -- stemP


% find vertex and corresponding edges for given sequence
vertex = StemP_find_vertex_PDB(rna,L,w,uu,LS1,LS2);
disp(['Number of vertex found  : ', num2str(size(vertex,1))])
edge = find_edge(vertex,pseudoknot);
disp(['Number of edges found   : ', num2str(size(edge,1))])
if ~isempty(edge) % case 1( there are some edges exist)
    
    [~,adjacency_matrix] = find_graph(edge);
    
    % 2. find cliques in the graph
    MC = Bron_Kerbosch(adjacency_matrix,2);
    cliques = find_cliques(MC);
    disp(['Number of cliques found : ', num2str(size(cliques,1))])
    % 3. find maximum matching from cliques
    [energy,~]= find_energy(cliques,vertex);
    
    % 4. compare the predictions with real data.
    Pred_alignments =  get_top_alignments(rna,energy,vertex);
    
    
else % case 2 ( there is no edge, then compare the energy ofeach vertex)
    Pred_alignments =  get_top_alignments(rna,energy,vertex);
end

end

