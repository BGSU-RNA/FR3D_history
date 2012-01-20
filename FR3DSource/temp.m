
warning off

[y,i] = sort(n(:,2));             % sort by number of nucleotides
t = t(i,:);
n = n(i,:);

for i = 1:length(n(:,1)),
  fprintf('BP/Nucl ratio %7.2f  PDBID: %s  Source: %s\n', n(i,3)/n(i,2), t{i,1}, t{i,8});
end
