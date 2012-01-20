% zCompare3DAlignment displays the nucleotides in the indicated rows of the
% alignment A, from File1 and File2

% One way to get A is to run these commands:
% [a,b,c] = xlsread('23S_Ec_Tt_Struct_alignment_12_4_06_JS.xls')


function [void] = zCompare3DAlignment(File1,File2,A,rows)

NTList1 = cat(2,A(rows,1));
NTList2 = cat(2,A(rows,2));

VP.Sugar = 1;

zSuperimposeNucleotides(File1,NTList1,File2,NTList2,VP);
