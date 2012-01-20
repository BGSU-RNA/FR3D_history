% zDistance(A,B) finds the distances between the rows of A and of B

function [D] = zDistance(A,B)

a = sum((A.*A)');                  % sum of squares of each row of A
b = sum((B.*B)');                  % sum of squares of each row of B

M = length(a);                     % number of rows in A
N = length(b);                     % number of rows in B

X = repmat(b,M,1);                 % repeat a in each row
Y = repmat(a',1,N);                % repeat b in each column

G = A * B';                        % inner products of rows of A

D = sqrt(X + Y - 2*G);             % |u-v| = |u|^2 + |v|^2 - 2 u . v
