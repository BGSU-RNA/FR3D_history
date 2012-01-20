
% define the triple you want to construct

%Pair1 = 'tHS';
%Pair2 = 'cWW';

%Base1 = 'A';
%Base2 = 'G';
%Base3 = 'U';

% How to use it:
% F = zMakeTriple('tHS','cWW','A','G','U');
% zWritePDB(F,'tHS_cWW_AGU_triple.pdb');

function [File] = zMakeTriple(Pair1,Pair2,Base1,Base2,Base3)

clear File

% retrieve examplar basepairs

[N1,N2] = zGetExemplar(Pair1,Base1,Base2);
[M2,M3] = zGetExemplar(Pair2,Base2,Base3);

% rotate the second pair so that the middle bases (M2, N2) superimpose

[s,t] = size(M2.Fit);

M2.Fit = (M2.Fit - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;

[s,t] = size(M3.Fit);

M3.Fit = (M3.Fit - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;

[s,t] = size(M2.Sugar);
M2.Sugar = (M2.Sugar - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;
M3.Sugar = (M3.Sugar - ones(s,1)*M2.Center)*M2.Rot*N2.Rot' + ones(s,1)*N2.Center;

% add the three nucleotides to a File data structure

File.NT(1) = N1;
File.NT(2) = N2;

%M3.Center = 
File.NT(3) = M3;
%File.NT(4) = M2;
File.NumNT = 3;

% display in a figure window

VP.Sugar = 1;
zDisplayNT(File,1:3,VP)
