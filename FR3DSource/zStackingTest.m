
[i,j] = find(triu((abs(File.Edge) > 20) .* (abs(File.Edge) < 24)));

% idea:  when you have two nucleotides, i and j, and you want to measure
% the degree of stacking overlap, calculate zStackingOverlap twice, once in
% each order, and average them.

for k = 1:length(i),
  so(k) = (zStackingOverlap(File.NT(i(k)),File.NT(j(k))) + zStackingOverlap(File.NT(j(k)),File.NT(i(k))))/2;
end

[y,p] = sort(so);

VP.AtOrigin = 1;

figure(1)
clf

hist(so,30)

figure(2)

for n = 1:length(i),
  fprintf('Stacking overlap is %5.2f, %s%s %s\n', so(p(n)), File.NT(i(p(n))).Base, File.NT(j(p(n))).Base, zEdgeText(File.Edge(i(p(n)),j(p(n)))));
  clf
  zDisplayNT(File, [i(p(n)) j(p(n))],VP);
  view(2);
  drawnow
  pause
end
