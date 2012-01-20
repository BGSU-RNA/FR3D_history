% xListCandidates(Search) prints a candidate list to the screen
% The optional argument NumToOutput limits the list's length
% The optional argument WheretoOutput has this effect:
%   Value 1 : prints a wide listing to the Matlab command window
%   Value 2 : prints a wide listing to an Editbox
%   Value 3 : prints a narrow listing to an Editbox

% It may be run directly from Matlab using the command:
%   xListCandidates(Search);

function [void] = xListCandidates(Search,NumToOutput,WheretoOutput)

File        = Search.File;

Query       = Search.Query;
Candidates  = Search.Candidates;

[s,t]       = size(Candidates);
N           = Query.NumNT;

if s == 0,
  fprintf('There are no candidates to list\n');
  return
end

if N == 2,
  CP = zeros(1,s);
end

if nargin < 2,
  NumToOutput = Inf;                    % limit on number printed to screen
end

if nargin < 3,
  WheretoOutput = 1;
end

%  For the PC compiled version, use this code:

%if nargin < 3,
%  WheretoOutput = 2;
%  xListCandidates(Search,NumToOutput,3);
%end

% -------------------------------------- print header line

Text{1} = '';

if isfield(Search,'AvgDisc'),
  Text{1} = [Text{1} sprintf('  Filename Avg Discrep ')];
  for i=1:N,
    Text{1} = [Text{1} sprintf('%7d ', i)];
  end
elseif Query.Geometric > 0,
  Text{1} = [Text{1} sprintf('  Filename Discrepancy ')];
  for i=1:N,
    Text{1} = [Text{1} sprintf('%7d ', i)];
  end
else
  Text{1} = [Text{1} sprintf('  Filename Number Nucleotides')];
end

c = 'Chains                                             ';
Text{1} = [Text{1} sprintf('%s', c(1:N))];

if WheretoOutput < 3,
  for i=1:N,
    for j=(i+1):N,
      Text{1} = [Text{1} sprintf('%6s', [num2str(i) '-' num2str(j)])];
    end
  end
  
  c = 'Configuration                                      ';
  Text{1} = [Text{1} sprintf(' %s', c(1:N))];
  
  for i=1:N,
    for j=(i+1):N,
      Text{1} = [Text{1} sprintf('%6s', [num2str(i) '-' num2str(j)])];
    end
  end
end   

if N == 2,
  Text{1} = [Text{1} sprintf(' Pair data')];
end
  
% -------------------------------------- list candidates

Config = {'A' , 'S'};

for i=1:min(s,NumToOutput),

  Text{i+1} = '';

  f = double(Candidates(i,N+1));               % file number for this candidate
  Text{i+1} = [Text{i+1} sprintf('%10s', File(f).Filename)];

  Indices = Candidates(i,1:N);                 % indices of nucleotides

  if isfield(Search,'AvgDisc'),
    Text{i+1} = [Text{i+1} sprintf('%12.4f',Search.AvgDisc(i))];
  elseif Query.Geometric > 0,
    Text{i+1} = [Text{i+1} sprintf('%12.4f',Search.Discrepancy(i))];
  else
    Text{i+1} = [Text{i+1} sprintf('%6d',Search.Discrepancy(i))];      % original candidate number
  end

  for j=1:N,
    Text{i+1} = [Text{i+1} sprintf('%3s',File(f).NT(Indices(j)).Base)];    
    Text{i+1} = [Text{i+1} sprintf('%5s',File(f).NT(Indices(j)).Number)];    
  end

  Text{i+1} = [Text{i+1} sprintf(' ')];

  for j=1:N,
    Text{i+1} = [Text{i+1} sprintf('%s',File(f).NT(Indices(j)).Chain)];
  end

  if WheretoOutput < 3,
    for k=1:length(Indices),
      for j=(k+1):length(Indices),
          Text{i+1} = [Text{i+1} sprintf('%6s', zEdgeText(File(f).Edge(Indices(k),Indices(j))))];
      end
    end
    
    Text{i+1} = [Text{i+1} sprintf(' ')];
    
    for k=1:length(Indices),
      Text{i+1} = [Text{i+1} sprintf('%c', Config{File(f).NT(Indices(k)).Syn+1})];
    end
    
    for k=1:length(Indices),
      for j=(k+1):length(Indices),
        Text{i+1} = [Text{i+1} sprintf('%6d', abs(double(Indices(k))-double(Indices(j))))];
      end
    end
  end
    
  if N == 2,                        % special treatment for basepairs
    CP(i) = norm(File(f).NT(Candidates(i,1)).Sugar(1,:) - ...
                          File(f).NT(Candidates(i,2)).Sugar(1,:));
    Text{i+1} = [Text{i+1} sprintf('   C1*-C1*: %8.4f', CP(i))];
    NT1 = File(f).NT(Candidates(i,1));
    NT2 = File(f).NT(Candidates(i,2));
    Edge  = full(File(f).Edge(Candidates(i,1),Candidates(i,2)));
    Text{i+1} = [Text{i+1} sprintf(' %s ', zEdgeText(Edge))];
    Text{i+1} = [Text{i+1} sprintf('%7.1f ', Edge)];
    SA = {'A', 'S'};
    Text{i+1} = [Text{i+1} sprintf('%c', SA{1+File(f).NT(Candidates(i,1)).Syn})];
    Text{i+1} = [Text{i+1} sprintf('%c', SA{1+File(f).NT(Candidates(i,2)).Syn})];
  end

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

% -------------------------------------- Display the listing

if WheretoOutput > 1,
  mEditbox(Text,'Listing of Candidates');
else
  for i=1:length(Text),
    fprintf('%s\n',Text{i});
  end
end
