# StemP (Python)

Code author: Mengyi Tang Rajchel

## Table of contents
1. [Summary](#Summary)
2. [Dependencies](#dependencies)
3. [Input arguments](#Input-arguments)
4. [Run StemP](#Run-StemP)
5. [Sample Output](#Sample-Output)

## Summary <a name="Summary"></a>

We provide source code of StemP using `Python` and some examples of implementing StemP on various types of sequences, including short RNA sequences with predetermined parameters, tRNA, 5s rRNA and unknown sequences.

## Dependencies <a name="dependencies"></a>
In order to run Stem (Python), it is required to have `Python3` and libraries `numpy, time, argparse, yaml` installed.

## Input arguments <a name="input"></a>
We provide separate `yaml` files in the folder `Python/configs` to store input arguments. 

## Run StemP  <a name="runstemp"></a>
Type  `cd Python`  in your terminal. Then type command `python main.py --config configs/$config file name$` to run StemP. For example, to run an example of a short rna sequence specified in `configs/example_short_seq.yaml`, use command `python main.py --config configs/example_short_seq.yaml`. You can also change the given sequence and parameters in `configs/example_short_seq.yaml`.

## Sample Output <a name="outputs"></a>
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