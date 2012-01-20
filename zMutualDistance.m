% zMutalDistance(A,L) finds the mutual distances between the rows of A and
% returns a sparse matrix D in which all entries are less than L

function [D] = zMutualDistance(A,L)

a = sum((A.*A)');                  % sum of squares of each row

N = length(a);

B = A * A';                        % inner products of rows of A

C = repmat(a,N,1);

D = C + C' - 2*B;

D = sparse(D .* (D < L^2));

D = sqrt(D);
