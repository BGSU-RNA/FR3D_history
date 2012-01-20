% This program needs revisions!  This version has garbage from a specific problem!

% pGroupCandidates clusters sequences and displays groups together

function [Group,Group2] = pGroupCandidates(File,Model,Candidates,NumGroups,Group,Group2)

if nargin < 4,
  NumGroups = min(ceil(length(Candidates)/2),8);
end

[s,t] = size(Candidates);
N = t-2;

if nargin < 6,
  Dist = zeros(s);

  for k=1:s,
    f1    = Candidates(k,N+1);
    c1.NT = File(f1).NT(Candidates(k,1:N));

    c1.Centers         = cat(1,c1.NT.Center);
    c1.WeightedCenter  = Model.LocWeight * c1.Centers / Model.NumNT;
    c1.CenteredCenters = c1.Centers-ones(N,1)*Model.WeightedCenter;
    c1.WeightedCenteredCenters = diag(Model.LocWeight)* c1.CenteredCenters;
    c1.LocWeight       = Model.LocWeight;
    c1.NumNT           = Model.NumNT;
    c1.LDiscCutoff     = Model.LDiscCutoff;
    c1.AngleWeight     = Model.AngleWeight;
    for j=1:k-1,
      f2 = Candidates(j,N+1);
      c2 = File(f2).NT(Candidates(j,1:N));     
      Dist(j,k) = xDiscrepancy(c1,c2);
      Dist(k,j) = Dist(j,k);
    end
  end

  [Group,Group2] = pCluster(Dist);
end

[y,i] = sort(Group2(NumGroups,:));

%diary(['xFindMotif_' num2str(NumGroups) '_group.txt']);

fprintf('Candidates grouped according to mutual discrepancy, %4d groups\n\n', NumGroups);

for g = 1:NumGroups,
  fprintf('\n');
  j = find(Group(NumGroups,i) == g);

  for k = 1:length(j),
    g = i(j(k));
    f = Candidates(g,N+1);
    fprintf('%15s', File(f).Filename);
    if Model.Geometric > 0,
      fprintf('%11.4f',Candidates(g,N+2));
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

%diary off
