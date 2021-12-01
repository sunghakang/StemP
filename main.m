% This file is created for predicting PDB, tRNA, or 5s rRNA sequences
clc
clear
close all

rna          = 'UGCUCCUAGUACGUAAGGACCGGAGUG';
seq_type     = 'PDB';
% 
% rna = 'GACCUGUUAGUUUAAUGGUAAAACGGUAGCCUCCAGAGUGAUUGAUACUGGUUCGAUUCCGGUAUAGGUCC';
% seq_type = 'tRNA';

rna = 'AUGCCGACGGUCAUAGGACGGGGGAAACACCCGGACUCAUUCCGAACCCGGAAGUUAAGCCCCGUUCCGUCCCGCACAGUACUGUGUUCCGAGAGGGCACGGGAACUGCGGGAACCGUCGGCUU';
seq_type = '5S-rRNA-Archaeal';




disp('JOB STARTS');
tic;
switch seq_type
    case 'PDB'
        % Default parameters
        L               = 3;     % minimum length of as Stem
        LS1             = 0;     % Lower bound of Stem-Loop Score
        LS2             = 0;     % Upper bound of Stem-Loop Score
        w               = 0;     % 0/1 existence of wooble pairing
        uu              = 0;     % 0/1 existence of U-U pairing
        p               = 0;     % 0/1 existence o0f Pseudo-knot
        Pred_alignments = stemP_pred_pdb(rna,L,w,uu,LS1,LS2,p);
    case 'tRNA'
        % Default parameters
        s2              =  3;   % the lowerbound of Stem-Loop score
        s3              =  5.4; % the upperbound of Stem-Loop score 
                                % s3 is 4.7 or 5.4 as default value depends 
                                % on the organism; 
        d1              = 12.0; % d of accemptor stem satisfies d >= d1
        d2              = 18.0; % d of accemptor stem satisfies d <= d2
        L               =  3.0; % minimum length of as Stem
        Pred_alignments = stemP_pred_tRNA(rna,L,d1,d2,s2,s3);
    case '5S-rRNA-Archaeal'
        % Default refined parameters is built-in
        Pred_alignments = stemP_pred_5srRNA(rna);
end
cputime = toc;


% print the  prediction results
fprintf('JOB COMPELETES \n');
fprintf(['CPU :           ',num2str(cputime), ' seconds \n'])
fprintf(['Sequence      : ', rna, '\n'])
fprintf(['Prediction    : \n'])
disp(Pred_alignments)
% there might be multiple output depends on if there are multiple top1
% energy
