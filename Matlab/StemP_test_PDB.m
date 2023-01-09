% Predict a secondary structure of PDB
% vertexmatrix : the cell type of vertex of true folding
% vertexoutput : the double type of vertex of true folding
% w : existence of wooble pairing
% u : existence of U-U pairing
% pseudoknot: Y/N existance of pseudoknot
% w: Y/N existance of wobble pairing G-U
% uu: Y/N existance of U-U pairing
clc
clear
L = 3;
LS1 = 0;
LS2 = 0;
w = 0;
uu = 0;
pseudoknot = 0;
% Case 1: The sequence is in the library
%         read true folding of PDB from excel.

seq_name = '1MSY'; % seq_name = '1MSY';
[rna,vertexmatrix,TrueVertex,vertexoutput,pseudoknot,w,uu] ...
    = readpdbseq(seq_name);

% Case 2: The sequence is not in the library, directly input the sequence
% rna = 'UGCUCCUAGUACGUAAGGACCGGAGUG';
% rna = 'GGCACAGAAGAUAUGGCUUCGUGCC';



% record CPU
tic;

if exist('seq_name','var')
    group1 = {'2AB4', '1MSY', '1L2X', '2AP5', '1F1T', '1XJR', '1DK1',...
        '1MMS', '1KXK', '2DU3', '2OIU'};
    group2 = {'1RLG',  '2OZB', '1F1T', '1SO3', '1XJR', '3E5C'};
    if ismember(seq_name, group1)
        LS1 = 2;LS2 = 20;
    end
    if ismember(seq_name, group2)
        L = 2;
    end
    
    True_alignment = vertexToalignment(rna,vertexmatrix);
end

vertex = StemP_find_vertex_PDB(rna,L,w,uu,LS1,LS2);
edge = find_edge(vertex,pseudoknot);

% Find cliques in the graph (2 cases acoording to the existence of edges)

if ~isempty(edge) % case 1( there are some edges exist)
    [graph,adjacency_matrix] = find_graph(edge);
    
    % 2. find cliques in the graph
    MC = Bron_Kerbosch(adjacency_matrix,2);
    cliques = find_cliques(MC);
    
    % 3. find maximum matching from cliques
    [energy,answer]= find_energy(cliques,vertex);
    
    % 4. compare the predictions with real data.
    if exist('TrueVertex', 'var') == 1
        [row_number_k,standard_rank,highest_MCC,first_MCC,alignments] ...
            = StemP_pred_accuracy(rna,energy, vertex,TrueVertex);
        Dense_rank = StemP_find_dense_rank(energy,row_number_k);
        top_k_energy = transpose(energy(1:row_number_k,2));
        Pred_alignment = alignments{end};
    else
        Pred_alignment =  get_top_alignments(rna,energy,vertex);
        Pred_alignment = Pred_alignment{1, 1} ;
        row_number_k =1;
    end
    
else % case 2 ( there is no edge, then compare the energy ofeach vertex)
    if exist('TrueVertex', 'var') == 1
        [energy,answer] = Stemp_find_match(vertex);
        [row_number_k,standard_rank,highest_MCC,first_MCC,alignments] ...
            = StemP_pred_accuracy(rna,energy, vertex,TrueVertex);
        Dense_rank = StemP_find_dense_rank(energy,row_number_k);
        top_k_energy = transpose(energy(1:row_number_k,2));
        Pred_alignment = alignments{end};
    else
        Pred_alignment = get_top_alignments(rna,energy,vertex);
        Pred_alignment = Pred_alignment{1, 1} ;
        row_number_k =1;
    end
end
cputime = toc;


%

fprintf('JOB COMPELETE \n');
fprintf(['CPU :           ',num2str(cputime), ' seconds \n'])

if exist('True_alignment','var') == 1
    fprintf(['Sequence name : ', seq_name, '\n'])
    fprintf(['Sequence      : ', rna, '\n'])
    fprintf(['True alignment: ', True_alignment, '\n'])
    fprintf(['Prediction    : ', Pred_alignment, '\n'])
    fprintf(['Row number    : ', num2str(row_number_k), '\n'])
    fprintf(['SCR (1224)    : ', num2str(standard_rank), '\n'])
    fprintf(['Dense rank    : ', num2str(Dense_rank), '\n'])
    fprintf(['Highest MCC   : ', num2str(highest_MCC), '\n'])
    fprintf(['First MCC     : ', num2str(first_MCC), '\n'])
else
    fprintf(['Sequence name : ', 'UNKNOWN', '\n'])
    fprintf(['Sequence      : ', rna, '\n'])
    fprintf(['True alignment: ', 'UNKNOWN', '\n'])
    fprintf(['Prediction    : ', Pred_alignment, '\n'])
    fprintf(['Row number    : ', num2str(row_number_k), '\n'])
    fprintf(['SCR (1224)    : ', num2str(row_number_k), '\n'])
    fprintf(['Dense_rank    : ', num2str(row_number_k), '\n'])
end
% [length(rna),dense_rank,first_MCC,row_number_k,standard_rank,pred_number3,cputime,w,uu,cell2mat(top_k_energy)];