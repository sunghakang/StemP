
# StemP: A fast and deterministic Stem-graph approach for RNA secondary structure prediction

Last update: January, 2023.

This repository provides source code of StemP, a fast algorithm to predict secondary strucutre of an RNA sequences. If you find this useful in your research, please consider cite our [paper](https://arxiv.org/abs/2201.05724):

``` 
@article{tang2022stemp,
  title={StemP: A fast and deterministic Stem-graph approach for RNA and protein folding prediction},
  author={Tang, Mengyi and Hwang, Kumbit and Kang, Sung Ha},
  journal={arXiv preprint arXiv:2201.05724},
  year={2022}
}
```

We provide both Python and Matlab versions of the source code. They are listed in folders `Python` and `Matlab` respectively. We also provide examples of implementing StemP on various sequences, including short RNA sequences with predetermined parameters, tRNA, 5s rRNA, and unknown sequences. 
```
STEMP
├───Matlab
└───Python
    ├───configs
    └───utils
```

## Table of contents
- [Abstract](#Abstract)
- [Run StemP Using Python](#Run-StemP-Using-Python)
   - [Dependencies](#dependencies)
   - [Input arguments](#Input-arguments)
   - [Run StemP](#Run-StemP)
   - [Sample Output](#Sample-Output)
- [Run StemP Using Python](#Run-StemP-Using-Python)
   - [Run StemP on a short RNA or unknown type sequence](#Run-StemP-on-a-short-RNA-or-unknown-type-sequence)
   - [Run StemP on a tRNA sequence](#Run-StemP-on-a-tRNA-sequence)
   - [Run StemP on a 5s rRNA sequence](#Run-StemP-on-a-5s-rRNA-sequence)
   - [A list of functions with descriptions](#A-list-of-functions-with-descriptions)


## Abstract

We propose a new deterministic methodology to predict the secondary structure of RNA sequences.  *What information of stem is important for structure prediction, and is it enough ?* The proposed simple deterministic algorithm uses  minimum stem length, Stem-Loop score, and co-existence of stems, to give good structure predictions for short RNA and tRNA sequences.}
The main idea is to consider all possible stem with certain stem loop energy and strength to predict RNA secondary structure.
We use graph notation, where stems are represented as vertexes, and co-existence between stems as edges.  This full Stem-graph presents all possible folding structure, and we pick sub-graph(s) which give the best matching energy for structure prediction.  Stem-Loop score adds structure information and speeds up the computation.  
The proposed method can predict secondary structure even with pseudo knots.  
One of the strengths of this approach is  the simplicity and flexibility of the algorithm, and it gives a  deterministic answer. 
Numerical experiments are done {on various sequences from Protein Data Bank and the Gutell Lab} using a laptop and results take only a few seconds.  


## Run StemP Using Python 

### Dependencies 
In order to run Stem (Python), it is required to have `Python3` and libraries `numpy, time, argparse, yaml` installed.

### Input arguments 
We provide separate `yaml` files in the folder `Python/configs` to store input arguments. 

### Run StemP  
Type  `cd Python`  in your terminal. Then type command `python main.py --config configs/$config file name$` to run StemP. For example, to run an example of a short rna sequence specified in `configs/example_short_seq.yaml`, use command `python main.py --config configs/example_short_seq.yaml`. You can also change the given sequence and parameters in `configs/example_short_seq.yaml`.

### Sample Output  
* Command: `python main.py --config configs/example_short_seq.yaml`
```
JOB STARTS ...
Number of vertices found :  6
Number of edges found  :  0
Number of cliques found:  6
Running time:  0.0  seconds.
JOB ENDS.
RNA                   :
UCCGAAGUGCAACGGGAAAAUGCACU
Predicted alignment(s):
.....((((((.........))))))
```

* Command: `python main.py --config configs/example_unknown.yaml`
```
JOB STARTS ...
Number of vertices found :  5
Number of edges found  :  4
Number of cliques found:  4
Running time:  0.0019669532775878906  seconds.
JOB ENDS.
RNA                   :
UGCUCCUAGUACGUAAGGACCGGAGUG
Predicted alignment(s):
..((((...[[[]]]......))))..
...((((..[[[]]]))))........
```

* Command: `python main.py --config configs/example_trna.yaml`
```
JOB STARTS ...
Number of vertices found :  30
Number of edges found  :  225
Number of cliques found:  300
Running time:  0.03787088394165039  seconds.
JOB ENDS.
RNA                   :
GACCUGUUAGUUUAAUGGUAAAACGGUAGCCUCCAGAGUGAUUGAUACUGGUUCGAUUCCGGUAUAGGUCC
Predicted alignment(s):
(((((((..[[[[.......]]]].........{{{...........}}}(((.......)))))))))).
(((((((..[[[[.......]]]]............{{{.......}}}.(((.......)))))))))).
```

* Command: `python main.py --config configs/example_5s_rrna.yaml`
```
JOB STARTS ...
Number of vertices found :  9
Number of edges found  :  20
Number of cliques found:  13
Running time:  0.1258082389831543  seconds.
JOB ENDS.
RNA                   :
AUGCCGACGGUCAUAGGACGGGGGAAACACCCGGACUCAUUCCGAACCCGGAAGUUAAGCCCCGUUCCGUCCCGCACAGUACUGUGUUCCGAGAGGGCACGGGAACUGCGGGAACCGUCGGCUU    
Predicted alignment(s):
.....((((((....[[[[[[[[.....[[[[[[.............]]]]..]]....]]]]]]]]..{{{{{{{.....{{{{{{.{{....}}}}}}}}....}}}}}}})))))).....    
.....((((((....[[[[[[[[.....[[[[[[.............]]]]..]]....]]]]]].]].{{{{{{{.....{{{{{{.{{....}}}}}}}}....}}}}}}})))))).....    
```


## Run StemP Using Matlab <a name="stempmatlab"></a>


### Run StemP on a short RNA (or unknown type) sequence <a name="shortrna"></a>
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

### Run StemP on a tRNA sequence <a name="trna"></a>
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

### Run StemP on a 5s rRNA sequence <a name="5srna"></a>
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

### A list of functions with descriptions. <a name="functions"></a>

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

