% pAdjustScore adjusts the 4x4 substitution matrix according to additional interactions the bases make

function [P] = pAdjustScore(File,i1,i2,P,method,ExclList)

if nargin < 6,
  ExclList = [];
end

g = File.Edge(i1,:);                         % interactions made by i1
g = g .* (abs(g) > 0) .* (abs(g) < 13);      % basepairs only
g(ExclList) = zeros(1,length(ExclList);      % ignore these interactions
g(i2) = 0;                                   % this one too
j = find(g);

for a = 1:length(j),
  Q = pIsoScore(File.Edge(i1,j(a)),File.NT(i1).Code,File.NT(j).Code,method)
  q = sum(Q,1) * ones(1,4)
  P = P .* q;
end

if length(j) > 0,
  pause
end

g = File.Edge(i2,:);                         % interactions made by i1
g = g .* (abs(g) > 0) .* (abs(g) < 13);      % basepairs only
g(ExclList) = zeros(1,length(ExclList);      % ignore these interactions
g(i1) = 0;                                   % this one too
j = find(g);

for a = 1:length(j),
  Q = pIsoScore(File.Edge(i2,j(a)),File.NT(i2).Code,File.NT(j).Code,method)
  q = sum(Q,1) * ones(1,4)
  P = P .* q';
end

if length(j) > 0,
  pause
end


P = P / sum(P);

