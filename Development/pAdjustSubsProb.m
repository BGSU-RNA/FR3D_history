% pAdjustScore adjusts the 4x4 substitution matrix according to additional interactions the bases make

function [R] = pAdjustSubsProb(File,i1,i2,P,method,ExclList)

Verbose = 1;

if nargin < 6,
  ExclList = [];
end

NT1 = File.NT(i1);
NT2 = File.NT(i2);

R = P;

g = File.Edge(i1,:);                         % basepair interactions made by i1
g = g .* (abs(g) > 0) .* (abs(g) < 13);      % basepairs only
g(ExclList) = zeros(1,length(ExclList));     % ignore these interactions
g(i2) = 0;                                   % this one too
j = find(g);

for a = 1:length(j),
  if Verbose > 1,
    fprintf('Long-range basepair %s%4s %s%4s %s\n', NT1.Base, NT1.Number, File.NT(j(a)).Base, File.NT(j(a)).Number, zEdgeText(File.Edge(i1,j(a))));

%    R
  end

  Q = pIsoScore(File.Edge(i1,j(a)),NT1.Code,File.NT(j(a)).Code,method);
  q = sum(Q,2) * ones(1,4);
  R = R .* q;                                % adjust base 1 probabilities
  
  if Verbose > 1,
    Q
    q
    R = R / sum(sum(R))
    pause
  end

end

g = File.Edge(i2,:);                         % basepair interactions made by i1
g = g .* (abs(g) > 0) .* (abs(g) < 13);      % basepairs only
g(ExclList) = zeros(1,length(ExclList));     % ignore these interactions
g(i1) = 0;                                   % this one too
j = find(g);

for a = 1:length(j),
  if Verbose > 1,
    fprintf('Long-range basepair %s%4s %s%4s %s\n', NT2.Base, NT2.Number, File.NT(j(a)).Base, File.NT(j(a)).Number, zEdgeText(File.Edge(i2,j(a))));

    R
  end

  Q = pIsoScore(File.Edge(i2,j(a)),NT2.Code,File.NT(j(a)).Code,method);
  q = sum(Q,2) * ones(1,4);
  R = R .* q';                               % adjust base 2 probabilities

  if Verbose > 1,
    Q
    q
    R = R / sum(sum(R))
    pause
  end
end

% -------------------------------------------- Adjust for BPh interactions

g = File.BasePhosphate(i1,:);                % BPh made by i1, i1 is the base
g = g .* (abs(g) > 0) .* (abs(g) < 100);     % true BPh only
g(i1) = 0;                                   % ignore self interactions
j = find(g);                                 % interaction partners

for a = 1:length(j),                         % loop through BPh inter
  bph = File.BasePhosphate(i1,j(a));

  if Verbose > 0,
    fprintf('BPh interaction %s%4s %s%4s %s\n', NT1.Base, NT1.Number, File.NT(j(a)).Base, File.NT(j(a)).Number, zBasePhosphateText(bph,1));
    R
  end

  Q = pBPhSpecificity(bph,1);                % 1 says it is in a basepair
  q = Q * ones(1,4);
  R = R .* q;                                % adjust base 1 probabilities

  if Verbose > 1,
    Q
    q
    R = R / sum(sum(R))
%    pause
  end
end

g = File.BasePhosphate(i2,:);                % BPh made by i2, i2 is the base
g = g .* (abs(g) > 0) .* (abs(g) < 100);     % true BPh only
g(i2) = 0;                                   % ignore self interactions
j = find(g);                                 % interaction partners

for a = 1:length(j),                         % loop through BPh inter
  bph = File.BasePhosphate(i2,j(a));
  if Verbose > 0,
    fprintf('BPh interaction %s%4s %s%4s %s\n', NT1.Base, NT1.Number, File.NT(j(a)).Base, File.NT(j(a)).Number, zBasePhosphateText(bph,1));
    R
  end

  Q = pBPhSpecificity(bph,1);                % 1 says it is in a basepair
  q = Q * ones(1,4);
  R = R .* q';                               % adjust base 2 probabilities


  if Verbose > 1,
    Q
    q
    R = R / sum(sum(R))
%    pause
  end
end



R = R / sum(sum(R));                         % normalize substitution probs

