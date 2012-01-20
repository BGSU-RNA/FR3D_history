
function [void] = pcolor(M)

[s,t] = size(M);

a = zeros(s,1);
b = zeros(1,t+1);

N = [[M a]; b];

pcolor(N);
