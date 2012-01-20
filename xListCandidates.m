% xListCandidates prints a candidate list to the screen

% It may be run directly from Matlab using the command:
%   [File] = xListCandidates([],Search);
% and after this,
%   [File] = xDisplayCandidates(File,Search);


function [File] = xListCandidates(File,Search,NumToOutput)

if isempty(File),
  [File,SIndex] = zAddNTData(Search.Filenames,2);   % load PDB data
else
  [File,SIndex] = zAddNTData(Search.Filenames,2,File); % add PDB data if needed
end

File = File(SIndex);                   % re-order file numbers

Query       = Search.Query;
Candidates  = Search.Candidates;

[s,t]       = size(Candidates);
N           = Query.NumNT;

if N == 2,
  CP = zeros(1,s);
end

if nargin < 3,
  NumToOutput = Inf;                    % limit on number printed to screen
end

% -------------------------------------- print header line
if isfield(Search,'AvgDisc'),
  fprintf('  Filename Avg Discrep  ');
  for i=1:N,
    fprintf('%7d ', i);
  end
elseif Query.Geometric > 0,
  fprintf('  Filename Discrepancy ');
  for i=1:N,
    fprintf('%7d ', i);
  end
else
  fprintf('  Filename Number Nucleotides');
end

c = 'Chains                                             ';
fprintf('%s', c(1:N));

for i=1:N,
  for j=(i+1):N,
    fprintf('   %d-%d', i,j);
  end
end

c = 'Configuration                                      ';
fprintf(' %s', c(1:N));

for i=1:N,
  for j=(i+1):N,
    fprintf('   %d-%d', i,j);
  end
end

if N == 2,
  fprintf(' Pair data');
end

fprintf('\n');

% -------------------------------------- list candidates

Config = {'A' , 'S'};

for i=1:min(s,NumToOutput),
  f = double(Candidates(i,N+1));               % file number for this candidate
  fprintf('%10s', File(f).Filename);

  Indices = Candidates(i,1:N);                 % indices of nucleotides

  if isfield(Search,'AvgDisc'),
    fprintf('%12.4f',Search.AvgDisc(i));
  elseif Query.Geometric > 0,
    fprintf('%12.4f',Search.Discrepancy(i));
  else
    fprintf('%6d',Search.Discrepancy(i));      % original candidate number
  end

  for j=1:N,
    fprintf('%3s',File(f).NT(Indices(j)).Base);    
    fprintf('%5s',File(f).NT(Indices(j)).Number);    
  end

  fprintf(' ');

  for j=1:N,
    fprintf('%s',File(f).NT(Indices(j)).Chain);
  end

  for k=1:length(Indices),
    for j=(k+1):length(Indices),
      fprintf('%6s', zEdgeText(File(f).Edge(Indices(k),Indices(j))));
    end
  end

  fprintf(' ');

  for k=1:length(Indices),
    fprintf('%c', Config{File(f).NT(Indices(k)).Syn+1});
  end

  for k=1:length(Indices),
    for j=(k+1):length(Indices),
      fprintf('%6d', abs(double(Indices(k))-double(Indices(j))));
    end
  end

  if N == 2,                        % special treatment for basepairs
    CP(i) = norm(File(f).NT(Candidates(i,1)).Sugar(1,:) - ...
                          File(f).NT(Candidates(i,2)).Sugar(1,:));
    fprintf('   C1*-C1*: %8.4f', CP(i));
    NT1 = File(f).NT(Candidates(i,1));
    NT2 = File(f).NT(Candidates(i,2));
    Edge  = full(File(f).Edge(Candidates(i,1),Candidates(i,2)));
    fprintf(' %s ', zEdgeText(Edge));
    fprintf('%7.1f ', Edge);
    SA = {'A', 'S'};
    fprintf('%c', SA{1+File(f).NT(Candidates(i,1)).Syn});
    fprintf('%c', SA{1+File(f).NT(Candidates(i,2)).Syn});
  end
  fprintf('\n');

end

% -------------------------------------- Additional notifications and info
if (Query.Geometric > 0),
  if (Query.RelCutoff > Query.DiscCutoff) && ~isfield(Search,'AvgDisc'),
    fprintf(fidOUT,'Some motifs with discrepancy between %7.4f and %7.4f might not appear above\n\n', Query.DiscCutoff, Query.RelCutoff);
  end
end

if s > NumToOutput,
  fprintf('Only the first %d candidates were listed.\n', NumToOutput);
end

if N == 2,
  figure(1)
  hist(CP,30)
  fprintf('Average C1''-C1'' distance is: %8.4f\n', mean(CP));
end
