% pAdjustScore adjusts the 4x4 substitution matrix according to additional interactions the bases make

function [P] = pAdjustScore(File,i1,i2,P,method,ExclList)

Verbose = 0;

if nargin < 6,
  ExclList = [];
end

g = File.Edge(i1,:);                         % interactions made by i1
g = g .* (abs(g) > 0) .* (abs(g) < 13);      % basepairs only
g(ExclList) = zeros(1,length(ExclList));      % ignore these interactions
g(i2) = 0;                                   % this one too
j = find(g);

for a = 1:length(j),
  Q = pIsoScore(File.Edge(i1,j(a)),File.NT(i1).Code,File.NT(j(a)).Code,method);
  q = sum(Q,2) * ones(1,4);
  P = P .* q;

if Verbose > 0,
  fprintf('LR Basepair %s%4s %s%4s %s\n', File.NT(i1).Base, File.NT(i1).Number, File.NT(j(a)).Base, File.NT(j(a)).Number, zEdgeText(File.Edge(i1,j(a))));

  P
  Q
  q
  P = P / sum(sum(P))
  pause
end

end

g = File.Edge(i2,:);                         % interactions made by i1
g = g .* (abs(g) > 0) .* (abs(g) < 13);      % basepairs only
g(ExclList) = zeros(1,length(ExclList));      % ignore these interactions
g(i1) = 0;                                   % this one too
j = find(g);

for a = 1:length(j),
  Q = pIsoScore(File.Edge(i2,j(a)),File.NT(i2).Code,File.NT(j(a)).Code,method);
  q = sum(Q,2) * ones(1,4);
  P = P .* q';

if Verbose > 0,
  fprintf('LR Basepair %s%4s %s%4s %s\n', File.NT(i2).Base, File.NT(i2).Number, File.NT(j(a)).Base, File.NT(j(a)).Number, zEdgeText(File.Edge(i2,j(a))));

  P
  Q
  q
  P = P / sum(sum(P))
  pause
end
end

P = P / sum(sum(P));

% -------------------------------------- Check for GU packing interaction

if isfield(File,'Motifs'),
if strfind(File.NT(i1).Nucl