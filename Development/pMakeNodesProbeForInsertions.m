% probe for insertions on left and right strands

aa = a;                                 % initial values of these
BB = B;

% ['Top ' File.NT(a).Number ' ' File.NT(B).Number]

LeftIns  = 0;                            % number of insertions on the left
RightIns = 0;                            % number of insertions on the right

% ---------------------------------------- Easy insertions on the left

while sum(sum(G(a,a:B))) == 0,          % a makes no interactions w/in loop
  if Verbose > 0,
    fprintf('    Insertion %4s      %s\n', File.NT(a).Number, File.NT(a).Base);
  end
  a = a + 1;
  LeftIns = LeftIns + 1;
end

% ---------------------------------------- Easy insertions on the right

while sum(sum(G(a:B,B))) == 0,          % B makes no interactions w/in loop
  if Verbose > 0,
    fprintf('    Insertion      %4s %s\n', File.NT(B).Number, File.NT(B).Base);
  end
  RightIns = RightIns + 1;         % increase mean number of insertions
  B = B - 1;                               % next base on right
end

% ['Done with easy insertions ' File.NT(a).Number ' ' File.NT(B).Number]

% ---------------------------------------- a interacts, but skip anyway?

aaa = min([length(File.NT) a+cdepth floor((a+B)/2)]);
LS  = (a+1):aaa;                        % left strand not including a
BBB = max([1 B-cdepth floor((a+B)/2)+1]);  % how far down right strand to look
RS  = BBB:B;                            % right strand including B
aB  = a:B;                              % indices from a to B
na  = find(H(a,aB));                       % nested interaction(s) a makes
nB  = find(H(B,aB));                       % nested interaction(s) B makes

if ~isempty(na),
%  ['a interacts with ' File.NT(aB(max(na))).Number]
end

if isempty(na),                             % a makes no nested interactions
  movea = 1;
elseif sum(sum(H((a+1):aB(max(na)),(a+1):(aB(max(na)))))) > 0, % 
  movea = 0;
%  disp('Found that the interaction a makes contains a stem');
elseif ~isempty(nB),                        % B makes a nested interaction
  if min(abs(nonzeros(G(a,aB(na))))) < min(abs(nonzeros(G(B,aB(nB))))),
    movea = 0;
%    disp('Both a and B make nested interactions, but a makes a more important one');
  end
elseif any(abs(G(a,aB(na)))==1),
  movea = 0;
else
  movea = 1;
end

%movea
%full(nonzeros(abs(G(a,aB))))
%full(nonzeros(abs(G(B,aB))))

while sum(abs(G(a,[LS RS]))) == 0 && movea == 1,
                                        % no interaction w/in this loop
  if Verbose > 0,
    fprintf('    Insertion %4s      %s\n', File.NT(a).Number, File.NT(a).Base);
  end
  LeftIns = LeftIns + 1;           % increase mean number of insertions
  a = a + 1;                              % next base on left
  aaa = min([length(File.NT) a+cdepth floor((a+B)/2)]);
  BBB = max([1 B-cdepth floor((a+B)/2)+1]);% how far down right strand to look
  LS  = (a+1):aaa;                        % left strand
  RS  = BBB:B;                            % right strand
  aB  = a:B;                            % indices from a to B
  na  = find(H(a,aB));                       % nested interaction(s) a makes
  nB  = find(H(B,aB));                       % nested interaction(s) B makes

  if isempty(na),                             % a makes no nested interactions
    movea = 1;
  elseif sum(sum(H((a+1):aB(max(na)),(a+1):(aB(max(na)))))) > 0, % 
    movea = 0;
%    disp('Found that the interaction a makes contains a stem');
  elseif ~isempty(nB),                        % B makes a nested interaction
    if min(abs(nonzeros(G(a,aB(na))))) < min(abs(nonzeros(G(B,aB(nB))))),
      movea = 0;
%      disp('Both a and B make nested interactions, but a makes a more important one');
    end
  elseif any(abs(G(a,aB(na)))==1),
    movea = 0;
  else
    movea = 1;
  end

%movea
%abs(G(a,aB))
%abs(G(B,aB))

end

%['Moved a ' File.NT(a).Number ' ' File.NT(B).Number]

% ---------------------------------------- Easy insertions on the right

while sum(sum(G(a:B,B))) == 0,          % B makes no interactions w/in loop
  if Verbose > 0,
    fprintf('    Insertion      %4s %s\n', File.NT(B).Number, File.NT(B).Base);
  end
  RightIns = RightIns + 1;         % increase mean number of insertions
  B = B - 1;                               % next base on right
end

% --------------------------------------- B interacts, but skip anyway?

aaa = min([length(File.NT) a+cdepth floor((a+B)/2)]);
LS  = a:aaa;                            % left strand including a
BBB = max([1 B-cdepth floor((a+B)/2)+1]);  % how far down right strand to look
RS  = BBB:(B-1);                        % right strand not including B
aB  = a:B;                              % indices from a to B
na  = find(H(a,aB));                       % nested interaction(s) a makes
nB  = find(H(B,aB));                       % nested interaction(s) B makes

if ~isempty(nB),
%  ['B interacts with ' File.NT(aB(max(nB))).Number]
end

if isempty(nB),                             % a makes no nested interactions
  moveB = 1;
elseif sum(sum(H(aB(min(nB)):(B-1),aB(min(nB)):(B-1)))) > 0, % 
  moveB = 0;
%  disp('Found that the interaction B makes contains a stem');
elseif ~isempty(nB),                        % B makes a nested interaction
  if min(abs(nonzeros(G(B,aB(nB))))) < min(abs(nonzeros(G(a,aB(na))))),
    moveB = 0;
%    disp('Both a and B make nested interactions, but B makes a more important one');
  end
elseif any(abs(G(B,aB(nB)))==1),
  moveB = 0;
else
  moveB = 1;
end

while sum((G([LS RS],B))) == 0 && moveB == 1,
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
  aB  = a:B;                            % indices from a to B
  na  = find(H(a,aB));                       % nested interaction(s) a makes
  nB  = find(H(B,aB));                       % nested interaction(s) B makes

  if ~isempty(nB),
%    ['B interacts with ' File.NT(aB(max(nB))).Number]
  end

  if isempty(nB),                             % a makes no nested interactions
    moveB = 1;
  elseif sum(sum(H(aB(min(nB)):(B-1),aB(min(nB)):(B-1)))) > 0, % 
    moveB = 0;
%    disp('Found that the interaction B makes contains a stem');
  elseif ~isempty(nB),                        % B makes a nested interaction
    if min(abs(nonzeros(G(B,aB(nB))))) < min(abs(nonzeros(G(a,aB(na))))),
      moveB = 0;
%      disp('Both a and B make nested interactions, but B makes a more important one');
    end
  elseif any(abs(G(B,aB(nB)))==1),
    moveB = 0;
  else
    moveB = 1;
  end

end

%['Moved B ' File.NT(a).Number ' ' File.NT(B).Number]

 
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


