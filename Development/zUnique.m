% zUniqueRows(M) returns the unique rows of M in b and counts of them in t
function [b,t] = zUniqueRows(M)

N = sortrows(M);

t(1)   = 1;
n = N(1,:);
c = 1;
b(c,:) = N(c,:);

for r = 2:length(N(:,1)),
  if all(n == N(r,:)),                  % same as previous row
    t(c) = t(c) + 1;                    % add a count for it
  else
    n = N(r,:);                         % new current row
    c = c + 1;
    b(c,:) = n;
    t(c) = 1;

if 0 > 1,
  b
  t
  pause
end

  end
end

[y,i] = sort(-t);

b = b(i,:);                             % most populous first
t = t(i);

% [b,i,j] = unique(M,'rows');


