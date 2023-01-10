# StemP: A fast and deterministic Stem-graph approach for RNA secondary structure prediction

## Table of contents
- [Summary](#Summary)
- [Run StemP on a short RNA or unknown type sequence](#Run-StemP-on-a-short-RNA-or-unknown-type-sequence)
- [Run StemP on a tRNA sequence](#Run-StemP-on-a-tRNA-sequence)
- [Run StemP on a 5s rRNA sequence](#Run-StemP-on-a-5s-rRNA-sequence)
- [A list of functions with descriptions](#A-list-of-functions-with-descriptions)

## Summary <a name="Summary"></a>

We provide source code of StemP in Matlab and some examples of implementing StemP on various types of sequences using Matlab. 

## Run StemP on a short RNA (or unknown type) sequence <a name="shortrna"></a>
Run `Matlab/main.m` with specified RNA sequence in Line 7. Your output in the command window should be: 
```
JOB STARTS
Number of vertex found  : 5
Number of edges found   : 4
Number of cliques found : 4
JOB COMPELETES 
CPU :           0.71935 seconds 
Sequence      : UGCUCCUAGUACGUAAGGACCGGAGUG
Prediction    : 
    {'..((((...[[[]]]......))))..'}
    {'...((((..[[[]]]))))........'}
```

## Run StemP on a tRNA sequence <a name="trna"></a>
Uncomment line 10-11 in `Matlab/main.m` and run it. Your output in the command window should be:
```
JOB STARTS
Number of vertex found  : 30
Number of edges found   : 225
Number of cliques found : 300
JOB COMPELETES 
CPU :           0.090718 seconds 
Sequence      : GACCUGUUAGUUUAAUGGUAAAACGGUAGCCUCCAGAGUGAUUGAUACUGGUUCGAUUCCGGUAUAGGUCC
Prediction    : 
    {'(((((((..[[[[.......]]]].........{{{...........}}}(((.......)))))))))).'}
    {'(((((((..[[[[.......]]]]............{{{.......}}}.(((.......)))))))))).'}
```

## Run StemP on a 5s rRNA sequence <a name="5srna"></a>
Uncomment line 13-14 in `Matlab/main.m` and run it. Your output in the command window should be:
```
JOB STARTS
Number of vertex found  : 9
Number of edges found   : 20
Number of cliques found : 13
JOB COMPELETES 
CPU :           0.057272 seconds 
Sequence      : AUGCCGACGGUCAUAGGACGGGGGAAACACCCGGACUCAUUCCGAACCCGGAAGUUAAGCCCCGUUCCGUCCCGCACAGUACUGUGUUCCGAGAGGGCACGGGAACUGCGGGAACCGUCGGCUU
Prediction    : 
    {'.....((((((....[[[[[[[[.....[[[[[[.............]]]]..]]....]]]]]]]]..{{{{{{{.....{{{{{{.{{....}}}}}}}}....}}}}}}})))))).....'}
    {'.....((((((....[[[[[[[[.....[[[[[[.............]]]]..]]....]]]]]].]].{{{{{{{.....{{{{{{.{{....}}}}}}}}....}}}}}}})))))).....'}
```

## A list of functions with descriptions. <a name="functions"></a>

`Matlab/main.py` allows you to input a sequence in line 7, such as 'AGUGUGCG' and its type (`'PDB', 'tRNA', '5S-rRNA-Archaeal'`), and output the alignment(s) with the highest energy. If a sequence has an unknown type, the input sequence type can be considered `PDB` as if it's a general sequence with an unknown structure.

`StemP_test_PDB.m`         

This code allows you to input the name of a short RNA sequence in the library or the RNA sequence that is not in the library. Output the corresponding alignments with the highest MCC values with the actual (true) folding (if the sequence is in the library).

We also provide all the helper functions for StemP.

`base_pairing.m`

This function determines whether two bases form canonical base pairs.


`Bron_Kerbosch.m`

This function outputs Maximal Cliques by using Bron-Kerbosch Algorithm.


`CliqueToVertex.m`

This function outputs a set of vertice in a clique.


`datacheckbyvertex2.m`

This function outputs the sensitivity / ppv / mcc / acc given true matched bases and predicted matched bases.


`find_basepairingtrue.m`

This function returns whether or not the m-th base matches the n-th base.

`find_cliques.m`

This function takes the graph as input and returns a column of cliques.


`find_edges.m`

This function takes the vertex and pseudoknots(0/1) as input and outputs a set of edges.

`find_energy.m`

This function takes cliques and corresponding bases in each vertex as input and outputs a list of energy in descending order together with the related maximum energy

`find_graph.m`

This function takes the edges as input and creates a graph with an adjacency matrix as output.

`find_rank.m`

This function takes the row number of a descending energy list and outputs the Standard Rank (SDR).

`get_top_alignments.m`

This function outputs the unique or multiple predicted alignments with the highest energy.

`PDB_data2.xlsx`

This file stores the PDB true alignments from the Protein Data Bank.

`pred_alignment.m`

This function outputs the predicted alignments of energy >= the energy of an alignment with the highest MCC value compared to the true folding structure.

`readpdbseq.m`

This function helps read the true folding info from PDB_data2.xlsx.

`StemP_find_dense_rank.m`

This function helps find the dense rank number given by the row number.


`StemP_find_edge_5s.m`

This function finds edges between three domains in 5s rRNA.