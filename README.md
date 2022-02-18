StemP: an algorithm to predicti RNA secondary structure

Copyright 2022, All Rights Reserved Code by Mengyi Tang, Kumbit Hwang, and Sung Ha Kang For Paper, "StemP: A fast and deterministic Stem-graph approach for RNA and protein folding prediction


main.m

     -- This code allows you to input the squence such as 'AGUGUGCG' and it's type, output the alignment(s) that has the highest envergy.


StemP_test_PDB.m              

     -- This code allows you to input the name of a pdb sequence that's in the libray, or the rna sequence that is not in the library. Output the corresponding alignments that has the highest MCC values with the true folding(if the sequence is in the library).

     -- All the other functions are helper functions.


I will keep updating this file.

base_pairing.m

     -- A helper function that determintes whether two bases forms canonical base pairs.


Bron_Kerbosch.m

     -- A helper function that outputs Maximal Cliques by using Bron-Kerbosch Algorithm


CliqueToVertex.m

     -- A helper function that outputs a set of vertex in a clique.


datacheckbyvertex2.m

     -- A helper function that output the sensitiviy/ppv/mcc/acc given true matched bases and predicted matched bases.


find_basepairingtrue.m

     -- A helper function that returns whether oor not m-th base matches n-th base.

find_cliques.m

     -- A helper function that takes the graph as input and returns a column of cliques.


find_edges.m

     -- A helper function that takes the vertex and pseudoknots(0/1) as input and output a set of edges.

find_energy.m

     -- A helper function that takes cliques and correponding bases in each vertex as input and ouput a list of energy in the descending order together with the corresponding maximum energy

find_graph.m

     -- A helper function that takes the edges as input and create a graph with adjacency matrix as output

find_rank.m

     -- A helper function that takes the row number of a descending energy list and ouput the Standard Rank (SDR).

get_top_alignments.m

     -- A helper function that outputs the unique or multiple prediced alignments that has the highest energy

PDB_data2.xlsx

     -- The excel file that stores the PDB true alignments from the Protein Data Bank.

pred_alignment.m

     -- A helper function that outputs the predicted alignments of energy >= the energy of an alignment that has highest MCC value compared to the true folding structure

readpdbseq.m

     -- A helper function that helps read the true folding info from PDB_data2.xlsx.

StemP_find_dense_rank.m

     -- A helper function that helps find the dense rank number given by row number.


StemP_find_edge_5s.m

     -- A helper function that finds edges between three domains in 5s rRNA

