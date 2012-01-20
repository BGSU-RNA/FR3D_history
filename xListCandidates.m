% xListCandidates prints a candidate list to the screen

function [] = xListCandidates(File,Search,NumToOutput)

Model       = Search.Query;
Candidates  = Search.Candidates;
Discrepancy = Search.Discrepancy;
[s,t]       = size(Candidates);
N           = Model.NumNT;

if N == 2,
  CP = zeros(1,s);
end

if nargin < 3,
  NumToOutput = 30;                    % limit on number printed to screen
end

Notify = 1;

% ------------------------------- Print to screen

if Model.Geometric > 0,
  fprintf('       Filename Discrepancy Nucleotides ... Chains\n');
else
  fprintf('       Filename Number Nucleotides ... Chains\n');
end

for i=1:s,
  if (Discrepancy(i) >= 0) & (i <= NumToOutput),
    f = double(Candidates(i,N+1));
    fprintf('%15s', File(f).Filename);
    if Model.Geometric > 0,
      fprintf('%11.4f',Discrepancy(i));
    else
      fprintf('%4d',Discrepancy(i));            % original candidate number
    end
    for j=1:N,
      fprintf('%3s',File(f).NT(Candidates(i,j)).Base);    
      fprintf('%4s',File(f).NT(Candidates(i,j)).Number);    
    end

    fprintf(' ');
    for j=1:N,
      fprintf('%s',File(f).NT(Candidates(i,j)).Chain);
    end

    if N == 2,                        % special treatment for basepairs
      CP(i) = norm(File(f).NT(Candidates(i,1)).Sugar(1,:) - ...
                            File(f).NT(Candidates(i,2)).Sugar(1,:));
      fprintf('   C1*-C1*: %8.4f', CP(i));
      NT1 = File(f).NT(Candidates(i,1));
      NT2 = File(f).NT(Candidates(i,2));
      Edge  = File(f).Edge(Candidates(i,1),Candidates(i,2));
      fprintf('%s ', zEdgeText(Edge));
    end
    fprintf('\n');
  end
end

if (Model.Geometric > 0) & (Model.RelCutoff > Model.DiscCutoff),
  fprintf(fidOUT,'Some motifs with discrepancy between %7.4f and %7.4f might not appear above\n\n', Model.DiscCutoff, Model.RelCutoff);
end

if s > NumToOutput,
  fprintf('Only the first %d candidates were listed.\n', NumToOutput);
end

if N == 2,
  figure(1)
  hist(CP,30)
  fprintf('Average C1''-C1'' distance is: %8.4f\n', mean(CP));
end
