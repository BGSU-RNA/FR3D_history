% zOrderGroups takes the result of a cluster analysis and puts the
% instances in order, consistent with the tree, so that similar branches are
% near each other in the list

% Y is the square matrix giving mutual distances
% Z is the linkage, telling what coalesces with what
% q is a permutation

function [q] = zOrderGroups(Y,Z,DD)

[s,t] = size(DD);

q = 1:s;                                % simplest ordering, doesn't do well

for i = 1:s,
  Group{i} = i;                         % initial groups have one element
end

gc = s;                                 % group counter

for z = 1:length(Z(:,1)),
  g  = Group{Z(z,1)};
  h  = Group{Z(z,2)};
  gg = fliplr(g);

  S(1) = Score(DD([g h],[g h]));
  S(2) = Score(DD([h g],[h g]));
  S(3) = Score(DD([gg h],[gg h]));
  S(4) = Score(DD([h gg],[h gg]));

  [y,i] = sort(S);

  switch i(1)
    case 1, Group{gc+1} = [g h];
    case 2, Group{gc+1} = [h g];
    case 3, Group{gc+1} = [gg h];
    case 4, Group{gc+1} = [h gg];
  end

  gc = gc + 1;
end

q = Group{end}


return




DDD = DD(q,q);                          % re-order according to p
DDDD = zeros(s+1,t+1);
DDDD(1:s,1:t) = DDD;

figure
pcolor(-DDDD)
shading flat
axis ij
view(2)

% -----------------------------------------------------------------------

function [S] = Score(M)

S = 2*sum(diag(M,1)) + sum(diag(M,2));