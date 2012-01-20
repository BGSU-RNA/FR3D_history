% zDistanceToExemplars(Exemplar,Pair) computes the distance to each
% exemplar for the given Pair

function [c,d,i] = zDistanceToExemplars(Exemplar,Pair)

pc = Pair.Paircode;

k = 1;

for j = 1:length(Exemplar(:,pc)),
  if ~isempty(Exemplar(j,pc).Class),
    c(k) = Exemplar(j,pc).Class;
    d(k) = zPairDiscrepancy(Exemplar(j,pc),Pair);
  else
    c(k) = 99;
    d(k) = 99999999;
  end
  k = k + 1;
end

[a,i] = sort(d);

c = c(i);
d = d(i);
