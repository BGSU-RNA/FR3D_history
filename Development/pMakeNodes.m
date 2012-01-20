% pMakeNodes(File,NTNumber,LastNTNumber,Truncate,Interact,Node,n) makes a secondary structure node model based on the Edge interaction matrix in File, starting at NTNumber and ending at LastNTNumber.  It assigns various nodes consistent with this secondary structure.  Truncate indicates where to put * hairpins.  Interact, Node, and n are optional parameters specified when pMakeNodes is called by itself.

function [Node] = pMakeNodes(File,Verbose,NTNumber,LastNTNumber,Truncate,Interact,Node,n)

if nargin < 2,
  Verbose = 1;
end

if nargin < 3,
  NTNumber = 1;
end

load PairExemplars

method = 2;                       % method for assigning pair subst probs

cdepth  = 10;                      % how far to look ahead for a cluster
jcdepth = 4;                      % how far to look for a junction cluster

if nargin < 5,
  Truncate = [];
end

if nargin < 8,
  n=0;                            % current node number
end

% -------------------------- if File is a text string (filename), load the file

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

% ----------------- if NTNumber is a cell array of numbers, look up the indices

if strcmp(class(NTNumber),'char'),
  NTNumber = {NTNumber};
end

if strcmp(class(NTNumber),'cell'),
  NTNumber = zIndexLookup(File,NTNumber);
end

% ----------------- if Truncate is a cell array of numbers, look up the indices

if strcmp(class(Truncate),'char'),
  Truncate = {Truncate};
end

if strcmp(class(Truncate),'cell'),
  if isempty(Truncate),
    Truncate = [];
  else
    Truncate = zIndexLookup(File,Truncate);
  end
end

% ------------------------------------------ prepare to identify motifs

HasMotif = zeros(1,length(File.NT));
if isfield(File,'Nucl'),
  for i = 1:length(File.NT),
    if ~isempty(File.Nucl(i).Motif),
      HasMotif(i) = 1;
    end
  end
end

% ------------------------------------------ Set key variables

N = length(File.NT);                       % number of nucleotides in File
E = abs(fix(File.Edge));                   % don't distinguish subcategories
G = E .* (E < 13) .* (E ~= 0);             % consider basepairing only
                                           % don't consider bifurcated now
H = G .* (File.Crossing == 0);             % nested pairs only

if nargin < 4,
  LastNTNumber = N;
elseif strcmp(class(LastNTNumber),'cell'),
  LastNTNumber = zIndexLookup(File,LastNTNumber);
end

% ------------------------------------------ Store indices of interacting bases

if nargin < 6,
 for a = 1:N,                              % loop through nucleotides
  k = find(G(a,:));                        % find indices of interacting bases
  [y,L] = sort(E(a,k));                    % sort by edge interaction category
  Interact{a}.Categ = abs(File.Edge(a,k(L)));   % store categories
  Interact{a}.Index = k(L);                % store indices of interacting bases
 end
end

% ------------------------------------------ Set up initial values of counters
a  = NTNumber;                             % first index; current index
A  = a;                                    % previous interacting base on left
AA = a;                                    % previous cWW base on left

B  = LastNTNumber;                         % next base on right
BB = a;                                    % previous cWW base on right

if Verbose > 0,
  fprintf('Loop %4s %4s\n', File.NT(a).Number, File.NT(B).Number);
end

% Initial node creation -------------------------------------------------

n = n+1;                                   % move to next node
Node(n).type      = 'Initial';             % node type
Node(n).nextnode  = n+1;                   % index of next node in tree
Node(n).lpar      = 0.01;                  % left insertion parameter
Node(n).rpar      = 0.01;                  % right insertion parameter
Node(n).LeftIndex = a;
Node(n).RightIndex= B;
Node(n).Insertion = [];                    % make sure this field exists
Node(n).id        = '';                    % make sure this field exists
Node(n).leftLengthDist = [];
Node(n).leftLetterDist = [];
Node(n).rightLengthDist = [];
Node(n).rightLetterDist = [];
Node(n).LeftLetter = '';
Node(n).RightLetter = '';
Node(n).Comment = '';
Node(n).NumLoops = [];
Node(n).Edge = [];
Node(n).Delete = [];
Node(n).SubsProb = [];
Node(n).Z = [];
Node(n).MiddleIndex = [];
Node(n).P = [];
Node(n).PIns = [];
Node(n).subtype = [];
Node(n).Left = [];
Node(n).Right = [];
Node(n).InsertionComment = '';
Node(n).IBases = [];
Node(n).InteractionComment = '';


                                           % probe for insertions
while sum(G(a,a+1:B)) == 0,                % no interaction w/in this loop
  Node(n).lpar = Node(n).lpar + 1;         % increase mean number of insertions
  a = a + 1;                               % move to next base
end

while sum(G(a:B-1,B)) == 0,                % no interaction w/in this loop
  Node(n).rpar = Node(n).rpar + 1;         % increase mean number of insertions
  B = B - 1;                               % move to next base
end

Node(n).leftLengthDist  = subPoisson(Node(n).lpar);
Node(n).leftLetterDist  = [0.25 0.25 0.25 0.25];
Node(n).rightLengthDist = subPoisson(Node(n).rpar); 
Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];

Node(n).LeftLetter  = cat(2,File.NT(Node(n).LeftIndex:(a-1)).Base);
Node(n).RightLetter = cat(2,File.NT((B+1):Node(n).RightIndex).Base);

Node(n).Comment = [' // Initial node ' File.NT(a).Base File.NT(a).Number ' - ' File.NT(B).Base File.NT(B).Number ' from junction ' Node(n).id];

if Verbose > 0,
  fprintf('%3d Initial %s:%s and %s:%s\n', n, File.NT(Node(n).LeftIndex).Number, File.NT(a).Number, File.NT(Node(n).RightIndex).Number, File.NT(B).Number);
end

% ---------------------------------------------------------------------------

EndLoop = 0;                               % flag for the end of the loop

while (EndLoop == 0) & (a <= LastNTNumber), % while not the end of the loop,

  b = Interact{a}.Index(1);                % index of what a interacts with

  ii = 2;
 
  while (b < a) && (ii <= length(Interact{a}.Index)),
    fprintf('Nucleotide %s%s interacts with %s%s\n',File.NT(a).Base,File.NT(a).Number,File.NT(b).Base,File.NT(b).Number);
    fprintf('Trying another interaction for %s%s =================================================================================\n', File.NT(a).Base,File.NT(a).Number);
    b = Interact{a}.Index(ii);
    ii = ii + 1;
  end    

  if (a < b),                            % if b comes after a

    % ---------------------------------- Check for junction

    % check to see if a and b are now in a new nested loop

    r = a;                               % leftmost index of a cWW
    rr= a;                               % start of current loop

    while sum(H(r,(r+1):B) > 0) == 0 && r < B,
                                         % if r does not make a nested pair,
      r = r + 1;
    end
    
    s = Interact{r}.Index(1);            % what it interacts with
    t = s+1;                             % next after that
    u = B;                               % end of current known loop

    if (sum(sum(H(t:u,t:u) == 1)) > 0) && ...
       (sum(sum(H(r:s,r:s) == 1)) > 0),
                    % there are helices between r and s and between t and u
      pMakeNodesJunction                   % identify and set up junctions
      return                               % nothing left to do!

    else                                   % not a junction

      % ---------------------------------- Identify basepair or cluster

      BBB = max(1,B-cdepth);
      aaa = min(length(File.NT),a+cdepth);

      if HasMotif(a),   % ---------------------- Insert motif model when needed
 
        pMakeNodesMotif

      elseif ((G(a,B) > 0) && sum(sum(G(a,BBB:(B-1)))) == 0 ...
                       && sum(sum(G((a+1):aaa,b:B))) == 0), 
                      % a and B interact, but not also with other nearby bases

        if b ~= B
          disp('b and B are different +++++++++++++++++++++++++++++++++++');
          b = B;
        end

        pMakeNodesBasepair

      else     % a and B also interact with nearby bases - use a cluster node

        pMakeNodesCluster

      end                                          % basepair or cluster

      % ------------------- check for truncation and hairpin and if not, probe for insertions

      if (EndLoop == 0),
        if ismember(a,Truncate) || ismember(a-1,Truncate) || isempty(File.NT(a).Base),

          pMakeNodesTruncate

        elseif (a == B) || ((sum(sum(G(a:B,a:B))) == 0)),

          pMakeNodesHairpin

        else                                 % probe for insertions

          pMakeNodesProbeForInsertions

        end                                  % hairpin or insertions
      end                                    % if EndLoop == 0
    end                                      % junction and junction cluster
  else
    fprintf('Nucleotide %s%s interacts with %s%s\n',File.NT(a).Base,File.NT(a).Number,File.NT(b).Base,File.NT(b).Number);
    fprintf('Skipping this nucleotide and moving to the next ================================================================================\n');
    a = a + 1;
  end                                   % if (a < b)

  if (ismember(a,Truncate) || isempty(File.NT(a).Base)) && (EndLoop == 0) , % cap with a * hairpin

    pMakeNodesTruncate

  end
end                                       % while (EndLoop == 0) & (a <= N),

% ---------------------------------- Poisson distribution for lengths -------

function [d] = subPoisson(m)

n = max(3,2*m);

d = exp(-m) * (m .^ (0:n)) ./ factorial(0:n);

d = d / sum(d);                     % normalize
