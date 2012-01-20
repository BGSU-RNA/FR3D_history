% pMakeModelFromSearchSaveFile(Search) creates an SCFG/MRF Node variable corresponding to the model in Search

% pMakeModelFromSearchSaveFile('LIB00002 IL 2008-03-20_23_29_25-Sarcin_13_flanked_by_cWW_in_1s72')
% Search = 'LIB00002 IL 2008-03-20_23_29_25-Sarcin_13_flanked_by_cWW_in_1s72';

% load LIB00014_IL_tSH-tSH-tHS-tHS.mat
% pMakeModelFromSearchSaveFile(Search,'IL',1);

function [Node,Truncate] = pMakeModelFromSearchSaveFile(Search,Type,Verbose)

if nargin < 3,
  Verbose = 0;
end

% ----------------------------------- Load Search from filename, if applicable

if strcmp(class(Search),'char'),
  load(['MotifLibrary' filesep Search],'Search','-mat');
end

% ----------------------------------- Gather basic information about the search

[L,N] = size(Search.Candidates);        % L = num instances; N = num NT
N = N - 1;                              % number of nucleotides

f = Search.Candidates(:,N+1);           % file numbers of motifs

File = Search.File(f(1));               % file of query motif
NTNumber = double(Search.Candidates(1,1));     % index of first NT
LastNTNumber = double(Search.Candidates(1,N)); % index of last NT



% ----------------------------------- Find locations of truncations

Direction = 1;                          % which order to put nucleotides
[y,p] = sort(Direction*double(Search.Candidates(1,1:N))); % 

if isfield(Search.Query,'MaxDiffMat'),
  MaxDiff = diag(Search.Query.MaxDiffMat(p,p),1);
else
  MaxDiff = Inf*ones(1,N-1);
end

maxinsert = zeros(1,N-1);
for c = 1:L,
  maxinsert = max(maxinsert,abs(diff(double(Search.Candidates(c,1:N))))-1);
end

Truncate = [];

for n = 1:(N-1),
  if (MaxDiff(n) == Inf) || (maxinsert(n) > 5),   % if only few insertions
    Truncate = [Truncate n+1];                    % truncate the strand here
  end
end

% ---------------------------- Find consensus interaction list

for a = 1:N,                                    % first NT of possible pair
  for b = 1:N,                                  % second NT of possible pair
    e = [];                                     % record observed edges
    for c = 1:L,                                % run through candidates
      i = Search.Candidates(c,a);               % index of first nucleotide
      j = Search.Candidates(c,b);               % index of second nucleotide
      e = [e Search.File(f(c)).Edge(i,j)];      % append observed interaction
    end

    e = fix(e);                                 % round subcategories

    for d = 1:length(e),
      if any(e(d) == [-1 -2 -7 -8]),            % don't distinguish sign here
        e(d) = -e(d);
      end
    end

    if max(abs(e)) > 0 && min(abs(e)) < 30,     % there was some bp interaction
%  e
      F.Edge(a,b) = mode(e);                        % use the most common one
  % a more sophisticated method for determining the consensus is needed
  % count how many times each one occurs, giving maybe 0.5 for near
    end
  end
end

% ---------------------------- Make the model for the consensus structure

i = Search.Candidates(1,1:N);                   % indices of query motif

if Verbose > 0,
  fprintf('Query motif interactions:\n');
  full(fix(File.Edge(i,i)))

  fprintf('Consensus interactions:\n')
  full(F.Edge)                                    % consensus interactions
end

% File.Edge(i,i) = F.Edge;                  % substitute consensus
F.NT = File.NT(Search.Candidates(1,1:N));   % use the first candidate as model

%Node = pMakeNodes(File,NTNumber,LastNTNumber,Truncate);
Node = pMakeNodes(F,1,N,Truncate);          % make the SCFG/MRF model

for m = 1:length(Node),
%  Node(m)
end

% ---------------------------- Shift indices down to 1

% This should not be necessary since the model is based on F, not File

if 0 > 1,
a = Node(1).LeftIndex;
b = Node(1).RightIndex-length(Node);

for n = 1:length(Node),
%  Node(n).LeftIndex  = Node(n).LeftIndex - a + 1;
%  Node(n).RightIndex = Node(n).RightIndex - b + 1;
%  Node(n).RightIndex = max(Node(n).RightIndex,max(Node(n).LeftIndex)+1);
%  Node(n).MiddleIndex = Node(n).MiddleIndex - a + 1;
end
end

% ---------------------------- Set parameters for the nodes from instances

for n = 1:length(Node),
  switch Node(n).type
  case 'Initial'

    % we should probably allow for stray bases at the beginning

  case 'Basepair'
    a = Node(n).LeftIndex;                   % which NT of the query motif
    b = Node(n).RightIndex;                  % right NT of the query motif
    Score = pConsensusPairSubstitution(a,b,f,Search.File,F,Node(n).Delete,L,Search,Verbose);

%fprintf('Original substitution probabilities\n');
%Node(n).SubsProb

%fprintf('Consensus substitution probabilities\n');
%Score

    Node(n).SubsProb = Score;

    % ----------------------------- tally insertions on the left
    inscount = [];
    letter = [0 0 0 0];                      % record which bases occur
    if n < length(Node),
      for c = 1:L,
        inscount(c) = Node(n+1).LeftIndex(1)-Node(n).LeftIndex(1)-1;
        
      end
    end

    lld = ones(1,max(inscount)+2);           % Dirichlet distribution
    for c = 1:L,
      lld(inscount(c)+1) = lld(inscount(c)+1) + 1;
    end

    Node(n).leftLengthDist = lld / sum(lld);    % normalize

    % ----------------------------- tally insertions on the right
    inscount = [];
    if n < length(Node),
      for c = 1:L,
        inscount(c) = Node(n).RightIndex(1)-Node(n+1).RightIndex(1)-1;
      end
    end

    rld = ones(1,max(inscount)+2);           % Dirichlet distribution
    for c = 1:L,
      rld(inscount(c)+1) = rld(inscount(c)+1) + 1;
    end

    Node(n).rightLengthDist = rld / sum(rld);    % normalize

  case 'Cluster'
    Indices = [Node(n).LeftIndex(Node(n).Left) ...
               Node(n).RightIndex(Node(n).Right)];
    for ii = 1:length(Node(n).IBases(:,1)),
      a = Indices(Node(n).IBases(ii,1));
      b = Indices(Node(n).IBases(ii,2));
      Score = pConsensusPairSubstitution(a,b,f,Search.File,F,0,L,Search,Verbose);
      Node(n).SubsProb(:,:,ii) = Score;
      if Verbose > 0,
        fprintf('\n');
      end
    end  

  case 'Junction'

  end

  fprintf('\n')
end
