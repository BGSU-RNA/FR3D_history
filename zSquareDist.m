% zDist find the mutual distances between the rows of A

function [D] = zDist(A)

a = sum((A.*A)');                  % sum of squares of each row

N = length(a);

B = A * A';                        % inner products of rows of A

C = repmat(a,N,1);

%D = ones(N,1) * a  +  a' * ones(1,N)  -  2*B;
D = C + C' - 2*B;