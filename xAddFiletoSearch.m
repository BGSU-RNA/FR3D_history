% xAddFiletoSearch(File,Search) adds data on nucleotides found in Search
% It assumes that File is in the proper order for the index in Candidates

function [Search] = xAddFiletoSearch(File,Search)

Query       = Search.Query;
Candidates  = Search.Candidates;

[s,t]       = size(Candidates);
N           = Query.NumNT;

if (s==0),
  Search.CandidateFilenames{1} = '';
  Search.File(1).Filename = '';
else

for i = 1:s,
  f = double(Candidates(i,N+1));        % file number for this candidate
  Search.CandidateFilenames{f} = Search.Filenames{f};
end

[y,p] = sort(double(Candidates(1,1:N)));    % put in increasing order

if isfield(Query,'MaxDiffMat'),
  MaxDiff = diag(Query.MaxDiffMat(p,p),1);
else
  MaxDiff = Inf*ones(1,N-1);
end

% ---------------------------- Calculate maximum gaps between cand. nucleotides

maxinsert = zeros(1,N-1);
for c = 1:s,
  maxinsert = max(maxinsert,abs(diff(double(Candidates(c,1:N))))-1);
end

% ---------------------------- Add nucleotide information

if ~isempty(File),
  for f = 1:max(Candidates(:,N+1)),
    Search.File(f).Edge = sparse(zeros(1,1));
    Search.File(f).BasePhosphate = sparse(zeros(1,1));
  end

  for i = 1:s,
    f = double(Candidates(i,N+1));             % file number for this candidate

    Search.CandidateFilenames{f} = File(f).Filename;
    Search.File(f).Filename = File(f).Filename;
    Search.File(f).NumNT    = File(f).NumNT;
    Search.File(f).Info     = File(f).Info;

    Indices = Candidates(i,1:N);               % indices of nucleotides

    for j = Indices,
      Search.File(f).NT(j) = File(f).NT(j);
      for k = Indices,
        Search.File(f).Edge(j,k) = File(f).Edge(j,k);
        Search.File(f).Edge(k,j) = File(f).Edge(k,j);
        Search.File(f).BasePhosphate(j,k) = File(f).BasePhosphate(j,k);
        Search.File(f).BasePhosphate(k,j) = File(f).BasePhosphate(k,j);
      end
    end

    % include intervening nucleotides, if only a few, for alignments

    for n = 1:(N-1),                      
      if (MaxDiff(n) < Inf) | (maxinsert(n) < 5),   % if only few insertions
      if Indices(n+1) - Indices(n) > 1,             % increasing order
        for i = (Indices(n)+1):(Indices(n+1)-1),
          Search.File(f).NT(i) = File(f).NT(i);
        end
      elseif Indices(n+1) - Indices(n) < -1,        % decreasing order
        for i = (Indices(n)-1):-1:(Indices(n+1)+1),
          Search.File(f).NT(i) = File(f).NT(i);
        end
      end
      end
    end


  end
end

end