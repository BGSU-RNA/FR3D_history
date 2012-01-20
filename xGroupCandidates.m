% xGroupCandidates clusters sequences and displays groups together
% The first time, run it as 
% [Group,Group2] = xGroupCandidates(File,Search);

function [Group,Group2] = xGroupCandidates(File,Search,NumGroups,Group,Group2)

Query      = Search.Query;
Candidates = Search.Candidates;

if nargin < 3,
  NumGroups = min(ceil(length(Candidates)/2),8);
end

[s,t] = size(Candidates);
N = Query.NumNT;

if ~isfield(Query,'LocWeight'),
  Query.LocWeight = ones(1,N);
end  

if ~isfield(Query,'AngleWeight'),
  Query.AngleWeight = ones(1,N);
end  

s = min(s,200);                              % don't do too many of these

if nargin < 4,
  Dist = zeros(s);                           % distances between candidates

  for k=1:s,
    f1    = Candidates(k,N+1);
    c1.NT = File(f1).NT(Candidates(k,1:N));

    c1.Centers         = cat(1,c1.NT.Center);
    c1.WeightedCenter  = Query.LocWeight * c1.Centers / Query.NumNT;
    c1.CenteredCenters = c1.Centers-ones(N,1)*c1.WeightedCenter;
    c1.WeightedCenteredCenters = diag(Query.LocWeight)* c1.CenteredCenters;
    c1.LocWeight       = Query.LocWeight;
    c1.NumNT           = Query.NumNT;
    c1.LDiscCutoff     = Inf;
    c1.AngleWeight     = Query.AngleWeight;
    for j=1:k-1,
      f2 = Candidates(j,N+1);
      c2 = File(f2).NT(Candidates(j,1:N));     
      Dist(j,k) = abs(xDiscrepancyFast(c1,c2));
      Dist(k,j) = Dist(j,k);
    end
  end

  [Group,Group2] = pCluster(Dist);
end

[y,i] = sort(Group2(NumGroups,:));

%diary(['xFindMotif_' num2str(NumGroups) '_group.txt']);

% ----------------------------------- Sort by centrality

  fprintf('Candidates sorted by centrality within these candidates:\n');

  [z,j] = sort(sum(Dist));

  S.Query       = Query;
  S.Candidates  = Candidates(j,:);
  S.Discrepancy = z;

  for k=1:s,
    g = j(k);
    f = Candidates(g,N+1);                 % file number
    fprintf('%15s', File(f).Filename);
    for jj=1:N
      fprintf('%3s',File(f).NT(Candidates(g,jj)).Base);    
      fprintf('%4s',File(f).NT(Candidates(g,jj)).Number);    
    end
    if N == 2,
      fprintf('   C1*-C1*: %8.4f', norm(File(f).NT(Candidates(g,1)).Sugar(1,:) - ...
                            File(f).NT(Candidates(g,2)).Sugar(1,:)));
    end
    fprintf('  %14.4f', z(k));
    fprintf('\n');
  end
  fprintf('\n');

% ----------------------------------- Form into groups

fprintf('Candidates grouped according to mutual discrepancy, %4d groups\n', NumGroups);

for g = 1:NumGroups,
  fprintf('\n');
  j = find(Group(NumGroups,i) == g);

  for k = 1:length(j),
    g = i(j(k));
    f = Candidates(g,N+1);
    fprintf('%15s', File(f).Filename);
    if Query.Geometric > 0,
      fprintf('%11.4f',Search.Discrepancy(g));
    end
    for jj=1:N
      fprintf('%3s',File(f).NT(Candidates(g,jj)).Base);    
      fprintf('%4s',File(f).NT(Candidates(g,jj)).Number);    
    end
    if N == 2,
      fprintf('   C1*-C1*: %8.4f', norm(File(f).NT(Candidates(g,1)).Sugar(1,:) - ...
                            File(f).NT(Candidates(g,2)).Sugar(1,:)));
    end
    fprintf('\n');
  end
end

  xDisplayCandidates(File,S,1);

%diary off
