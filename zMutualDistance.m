% zMutalDistance(A,L) finds the mutual distances between the rows of A and
% returns a sparse matrix D in which all entries are less than L

function [D] = zMutualDistance(A,L)

D = squareform(pdist(A));

D = sparse(D .* (D < L));
