% pMakeNodes(File,NTNumber,LastNTNumber,Truncate,Interact,Node,n) makes a secondary structure node model based on the Edge interaction matrix in File, starting at NTNumber and ending at LastNTNumber.  It assigns various nodes consistent with this secondary structure.  Truncate indicates where to put * hairpins.  Interact, Node, and n are optional parameters specified when pMakeNodes is called by itself.

function [Node] = pMakeNodes(File,Verbose,NTNumber,LastNTNumber,Truncate,Interact,Node,n)

if nargin < 2,
  Verbose = 1;
end

if nargin < 3,
  NTNumber = 1;
end

load PairExemplars

method = 4;                       % method for assigning pair subst probs

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
H = (G ~= 0) .* max(File.Crossing == 0, abs(G) == 1) ;

                                           % 1 for nested pairs, 0 otherwise

J = abs(G .* (File.Crossing >  0));        % long-range basepairs only

DelProb = 0.01;                            % nominal deletion probability
                                           % for basepairs
TertiaryFreeNode = 0;                      % first node in this stem making
                                           % no tertiary intearctions beyond it
                                           
if nargin < 4,
  LastNTNumber = N;
elseif strcmp(class(LastNTNumber),'cell'),
  LastNTNumber = zIndexLookup(File,LastNTNumber);
elseif strcmp(class(LastNTNumber),'char'),
  LastNTNumber = zIndexLookup(File,{LastNTNumber});
end

% ------------------------------------------ Store indices of interacting bases

if nargin < 6,
 for a = 1:N,                              % loop through nucleotides
  k = find(G(a,:));                        % find indices of interacting bases
%  k = find(H(a,:));                        % find indices of interacting bases
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

pMakeNodesNewNode;                         % set up blank node with all fields
Node(n).type      = 'Initial';             % node type
Node(n).nextnode  = n+1;                   % index of next node in tree
Node(n).LeftIndex = a;                     % index of first base on left
Node(n).RightIndex= B;                     % index of first base on right

pMakeNodesProbeForInsertions;              % probe for insertions, each strand

% ---------------------------------------------------------------------------

EndLoop = 0;                               % flag for the end of the loop

while (EndLoop == 0) & (a <= LastNTNumber), % while not the end of the loop,

    % ---------------------------------- Check for junction

    % check to see if a is now in a new nested loop

    r = a;                               % leftmost index of a cWW
    rr= a;                               % start of current loop

    while sum(H(r,(r+1):B)) == 0 && r < B, % if r does not make a nested pair,
      r = r + 1;
    end
    
    s = Interact{r}.Index(1);            % what it interacts with
    t = s+1;                             % next after that
    u = B;                               % end of current known loop

    if (sum(sum(G(t:u,t:u)==1)) > 0) && (sum(sum(G(r:s,r:s)==1)) > 0),
            % there are nested cWW pairs between r and s and between t and u
            % use cWW pairs to avoid pairs between i and i+1 making junctions

      if Verbose > 1,
        fprintf('Found nested interactions between %s and %s and between %s and %s\n', File.NT(r).Number, File.NT(s).Number, File.NT(t).Number, File.NT(u).Number);
      end

      pMakeNodesJunction                   % identify make models for junctions
      return                               % nothing left to do!

    else                                   % not a junction

      % ---------------------------------- Identify basepair or cluster

      aaa = min([length(File.NT) a+cdepth floor((a+B)/2)]);
      BBB = max([1 B-cdepth ceil((a+B)/2)]);
      LS = (a+1):aaa;                         % left strand
      RS = BBB:(B-1);                         % right strand

      if HasMotif(a),   % ---------------------- Insert motif model when needed
 
        pMakeNodesMotif
        pMakeNodesProbeForInsertions          % add Initial node if needed

      elseif ((H(a,B) > 0) && sum(sum(G(a,[LS RS]))) == 0 ...
                           && sum(sum(G([LS RS],B))) == 0), 
                      % a and B interact, but not also with other nearby bases
        pMakeNodesBasepair                       % add basepair with insertions

      else     % a and B also interact with nearby bases - use a cluster node

        % [a aaa BBB B]
        % zShowInteractionTable(File,unique([a:aaa BBB:B]));

        pMakeNodesCluster

        % pause

      end                                          % basepair or cluster

      % ------------------- check for truncation and hairpin

      if (EndLoop == 0),
        if ismember(a,Truncate) || ismember(a-1,Truncate) || isempty(File.NT(a).Base),

          pMakeNodesTruncate

        elseif (a == B) || ((sum(sum(G(a:B,a:B))) == 0)), % time for a hairpin
          if (TertiaryFreeNode > 0) && isempty(Truncate),
                     % extensible region, not a truncated model
            pMakeNodesExtraBasepairs;       % add extra basepairs if extensible
          end

          pMakeNodesHairpin
        else

          pMakeNodesProbeForInsertions       % add insertions if needed

        end                                  % hairpin or insertions
      end                                    % if EndLoop == 0
    end                                      % junction and junction cluster
end                                       % while (EndLoop == 0) & (a <= N),

% ---------------------------------- Poisson distribution for lengths -------

function [d] = subPoisson(m)

n = max(3,2*m);

d = exp(-m) * (m .^ (0:n)) ./ factorial(0:n);

d = d / sum(d);                     % normalize
