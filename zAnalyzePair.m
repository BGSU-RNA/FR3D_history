% zAnalyzePair(N1,N2) computes distances, angles, and classification codes
% for nucleotides N1 and N2, then returns a pair data type

function [Pair] = zAnalyzePair(N1,N2)

  CL = zClassLimits;                              % read ClassLimits matrix
  Pair = zAnalyzePairFast(N1,N2,CL);