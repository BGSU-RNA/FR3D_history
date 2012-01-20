
load PDBInfo

[y,i] = sort(n(:,2));            % sort by number of nucleotides

n = n(i,:);
t = t(i,:);

s = 1105;

match = [];

for b = 1:5,
  match(b,b) = 0.5;
  for c = (b+1):5,
    match(b,c) = nw(t{b+s,9},t{c+s,9},0.99999,2)/min(n(b+s,2),n(c+s,2));
  end
end

match = match + match';

n((s+1):(s+5),2)'

match