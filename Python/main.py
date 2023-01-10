from model import stemP_rna_structure_pred
import argparse
import yaml

parser = argparse.ArgumentParser(description='Implementation of StemP')
parser.add_argument('--config', default='configs/example_short_seq.yaml')

"""
Main function of StemP to predict RNA secondary structures.

Copyright 2023, All Rights Reserved.
Code by Mengyi Tang Rajchel
For Paper, "StemP: A fast and deterministic Stem-graph approach for RNA secondary structure prediction"
by Mengy Tang, Kumbit Hwang and Sung Ha Kang

"""


def main():
    global args
    args = parser.parse_args()
    with open(args.config) as f:
        config = yaml.safe_load(f)

    for name, value in config['sequence_information'].items():
        setattr(args, name, value)

    rna = args.rna
    seq_type = args.seq_type

    para = config['parameters']

    print('JOB STARTS ...')

    pred_alignments = stemP_rna_structure_pred(rna, seq_type, para)

    print('JOB ENDS.')
    print('RNA                   : ')
    print(rna)
    print('Predicted alignment(s): ')
    for i in range(len(pred_alignments)):
        print(pred_alignments[i])

    return


if __name__ == '__main__':
    main()