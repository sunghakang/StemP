% Run this file to predict short RNA sequences, tRNA, or 5s rRNA sequences on some examples.

clc
clear
close all

rna          = 'UGCUCCUAGUACGUAAGGACCGGAGUG';
seq_type     = 'PDB';
 
% rna = 'GACCUGUUAGUUUAAUGGUAAAACGGUAGCCUCCAGAGUGAUUGAUACUGGUUCGAUUCCGGUAUAGGUCC';
% seq_type = 'tRNA';
% 
% rna = 'AUGCCGACGGUCAUAGGACGGGGGAAACACCCGGACUCAUUCCGAACCCGGAAGUUAAGCCCCGUUCCGUCCCGCACAGUACUGUGUUCCGAGAGGGCACGGGAACUGCGGGAACCGUCGGCUU';
% seq_type = '5S-rRNA-Archaeal';




disp('JOB STARTS');
tic;
% Set up default parameters
switch seq_type
    case 'PDB'
        L               = 3;     % minimum length of as Stem
        LS1             = 0;     % Lower bound of Stem-Loop Score
        LS2             = 0;     % Upper bound of Stem-Loop Score
        w               = 0;     % 0/1 existence of wooble pairing
        uu              = 0;     % 0/1 existence of U-U pairing
        p               = 0;     % 0/1 existence o0f Pseudo-knot
        Pred_alignments = stemP_pred_pdb(rna,L,w,uu,LS1,LS2,p);
    case 'tRNA'
        s2              =  3;   % the lowerbound of Stem-Loop score
        s3              =  5.4; % the upperbound of Stem-Loop score 
                                % s3 is 4.7 or 5.4 as default value depending 
                                % on the organism; 
        d1              = 12.0; % minimum distance of an open-ended stem
        d2              = 18.0; % maximum distance of an open-ended stem
        L               =  3.0; % minimum length of a possible stem
        Pred_alignments = stemP_pred_tRNA(rna,L,d1,d2,s2,s3);
    case '5S-rRNA-Archaeal'
        Pred_alignments = stemP_pred_5srRNA(rna);
end
cputime = toc;


% print the  prediction results
fprintf('JOB COMPELETES \n');
fprintf(['CPU :           ',num2str(cputime), ' seconds \n'])
fprintf(['Sequence      : ', rna, '\n'])
fprintf(['Prediction    : \n'])
disp(Pred_alignments)
% Remark: There might be multiple outputs in the case when various cliques have the largest energy.
