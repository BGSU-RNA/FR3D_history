% probe for insertions on left and right strands

aa = a;                                 % initial values of these
BB = B;

LeftIns = 0;                            % number of insertions on the left
 
aaa = min([length(File.NT) a+cdepth floor((a+B)/2)]);
LS  = (a+1):aaa;                        % left strand
BBB = max([1 B-cdepth floor((a+B)/2)+1]);  % how far down right strand to look
RS  = BBB:B;                            % right strand

% File.NT(a).Number

while sum(abs(G(a,[LS RS]))) == 0 && sum(H(a,aa:BB)) == 0,
                                        % no interaction w/in this loop
  if Verbose > 0,
    fprintf('    Insertion %4s      %s\n', File.NT(a).Number, File.NT(a).Base);
  end
  LeftIns = LeftIns + 1;           % increase mean number of insertions
  a = a + 1;                               % next base on left
  aaa = min([length(File.NT) a+cdepth floor((a+B)/2)]);
  BBB = max([1 B-cdepth floor((a+B)/2)+1]);% how far down right strand to look
  LS  = (a+1):aaa;                        % left strand
  RS  = BBB:B;                            % right strand
end

% File.NT(a).Number

RightIns = 0;

aaa = min([length(File.NT) a+cdepth floor((a+B)/2)]);
LS  = a:aaa;                            % left strand
BBB = max([1 B-cdepth floor((a+B)/2)+1]);  % how far down right strand to look
RS  = BBB:(B-1);                        % right strand

% File.NT(B).Number
% G([LS RS],B)
% sum(H(B,:))

while sum((G([LS RS],B))) == 0 && sum(H(B,aa:BB)) == 0,
                                        % no interaction w/in this loop
  if Verbose > 0,
    fprintf('    Insertion      %4s %s\n', File.NT(B).Number, File.NT(B).Base);
  end
  RightIns = RightIns + 1;         % increase mean number of insertions
  B = B - 1;                               % next base on right
  aaa = min([length(File.NT) a+cdepth floor((a+B)/2)]);
  LS  = a:aaa;                             % left strand
  BBB = max([1 B-cdepth floor((a+B)/2)+1]);
  RS  = BBB:(B-1);                         % right strand
end

% File.NT(B).Number
 
if strcmp(Node(n).type,'Basepair'),   % add insertions to basepair
  Node(n).lpar = (0.01+LeftIns)  * [ones(16,1); 0];
  Node(n).rpar = (0.01+RightIns) * [ones(16,1); 0];
  Node(n).leftLengthDist  = subPoisson(0.01+LeftIns);
  Node(n).leftLetterDist  = [0.25 0.25 0.25 0.25];
  Node(n).rightLengthDist = subPoisson(0.01+RightIns); 
  Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];
  Node(n).LeftLetter  = [Node(n).LeftLetter cat(2,File.NT((Node(n).LeftIndex+1):(a-1)).Base)];
  Node(n).RightLetter = [cat(2,File.NT((B+1):(Node(n).RightIndex-1)).Base) Node(n).RightLetter];

elseif strcmp(Node(n).type,'Initial'),
  Node(n).lpar            = 0.01+LeftIns;       % left insertion parameter
  Node(n).rpar            = 0.01+RightIns;      % right insertion parameter
  Node(n).leftLengthDist  = subPoisson(Node(n).lpar);
  Node(n).leftLetterDist  = [0.25 0.25 0.25 0.25];
  Node(n).rightLengthDist = subPoisson(Node(n).rpar); 
  Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];

  Node(n).LeftLetter  = cat(2,File.NT(Node(n).LeftIndex:(a-1)).Base);
  Node(n).RightLetter = cat(2,File.NT((B+1):Node(n).RightIndex).Base);

elseif strcmp(Node(n).type,'Cluster') && (LeftIns > 0 || RightIns > 0),
  n = n + 1;                          % initial node after cluster
  Node(n).type = 'Initial';
  Node(n).nextnode  = n+1;            % index of next node in tree
  Node(n).lpar = LeftIns;
  Node(n).rpar = RightIns;
  Node(n).LeftIndex  = aa;
  Node(n).RightIndex = BB;
  Node(n).LeftLetter  = cat(2,File.NT(aa:(aa+LeftIns-1)).Base);
  Node(n).RightLetter = cat(2,File.NT((BB-RightIns+1):BB).Base);
  Node(n).leftLengthDist  = subPoisson(Node(n).lpar);
  Node(n).leftLetterDist  = [0.25 0.25 0.25 0.25];
  Node(n).rightLengthDist = subPoisson(Node(n).rpar); 
  Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];
end

if strcmp(Node(n).type,'Initial'),
  Node(n).Comment = [' // Initial node ' File.NT(aa).Base File.NT(aa).Number ' - ' File.NT(BB).Base File.NT(BB).Number ' from junction ' Node(n).id];

  if Verbose > 0,
    fprintf('%3d Initial   %4s (%d insertion) and %4s (%d insertion)\n', n, File.NT(Node(n).LeftIndex).Number, LeftIns, File.NT(Node(n).RightIndex).Number, RightIns);
  end
end

