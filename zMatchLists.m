
function [i,j] = zMatchLists(A,B)

% assume A and B are sorted vectors

a = 1;
b = 1;
i = [];
j = [];

while (a <= length(A)) && (b <= length(B)),
  if A(a) < B(b)
    a = a + 1;
  elseif A(a) > B(b)
    b = b + 1;
  else
    i = [i; a];
    j = [j; b];
    a = a + 1;
    b = b + 1;
  end
end

