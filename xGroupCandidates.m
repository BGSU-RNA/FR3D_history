% xGroupCandidates clusters sequences and displays groups together

function [Search] = xGroupCandidates(File,Search,NumGroups)

if nargin < 3,
  NumGroups = min(ceil(length(Candidates)/2),8);
end

Search = xMutualDiscrepancy(File,Search);
Done   = find(Search.DiscComputed);            % ones already computed

[Group,Group2] = pCluster(Dist);

[y,i] = sort(Group2(NumGroups,:));

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

