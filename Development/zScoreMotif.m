% zScoreMotif takes a set of nucleotides from a FR3D search, tallies up the observed frequency of the different sequences, and scores them according to various pairwise scoring schemes


% load 2008-06-06_13_44_49-Internal_loops_flanked_by_nested_cWW_pairs.mat
% load 2010-05-18_15_42_07-Sarcin_7_mixed_2aw4.mat

method = 2;

Filenames = Search.Filenames;

Verbose = 1;

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData(Filenames,2,[],Verbose);   % load PDB data
else
  [File,SIndex] = zAddNTData(Filenames,2,File,Verbose); %add PDB data if needed
end                       % SIndex tells which elements of File to search

File = File(SIndex);

Cand = Search.Candidates;

[L,N] = size(Cand);        % number of candidates
N = N - 1;                 % number of nucleotides

File = zAttachAlignment(File);              % attach alignment data

AllLett = [];
for c = 1:L,
  Lett = [];
  for n = 1:N,
    Lett = [Lett File.NT(Cand(c,n)).FASTA];
  end
  AllLett = [AllLett; Lett];
end

i = find(sum(AllLett == 'N',2) == 0);
AllLett = AllLett(i,:);

i = find(sum(AllLett == '-',2) == 0);
AllLett = AllLett(i,:);

[seqs,counts] = zUnique(AllLett)

S = length(seqs(:,1));

codes = 1*(seqs=='A') + 2*(seqs=='C') + 3*(seqs=='G') + 4*(seqs=='U');

% find consensus interaction list

Edge = pConsensusInteractions(Search);

% score motif according to basepairs

score = ones(S,1);          % initial scores

for p = 1:N,
  for q = 1:N,                              % loop through interactions
    if abs(Edge(p,q)) > 0 && abs(Edge(p,q)) < 13,  % basepair
      IS = pIsoScore(Edge(p,q),codes(1,p),codes(1,q),method); % use most common seq
      for s = 1:S,                          % loop through sequences
        score(s) = score(s) * IS(codes(s,p),codes(s,q));
      end
    end
  end
end

loglog(counts,score,'*')
hold on
for s = 1:S,
  text(counts(s),score(s),seqs(s,:),'fontsize',6);
end

% ------------------------------------------- score all possible sequences

allcodes = zeros(4^N,N);                     % all sequences of length N
c = ones(1,N);                               % start with all 1's
allcodes(1,:) = c;                           % first row
r = 2;                                       % current row
a = N;                                       % start with last digit

while a > 0,
  if c(a) < 4,
    c(a) = c(a) + 1;                         % move to next sequence
    allcodes(r,:) = c;                       % store this one
    r = r + 1;
  else
    while a >= 1 && c(a) == 4,
      c(a) = 1;                                % roll this one over
      a = a - 1;
    end
    if a > 0,
      c(a) = c(a) + 1;  
      allcodes(r,:) = c;
      r = r + 1;
      a = N;
    end
  end
end

allscores = ones(4^N,1);                    % place to store all scores

for p = 1:N,
  for q = 1:N,                              % loop through interactions
    if abs(Edge(p,q)) > 0 && abs(Edge(p,q)) < 13,  % basepair
      IS = pIsoScore(Edge(p,q),codes(1,p),codes(1,q),method); % use most common seq
      for s = 1:(4^N),                      % loop through all sequences
        allscores(s) = allscores(s) * IS(allcodes(s,p),allcodes(s,q));
      end
    end
  end
end

[y,i] = sort(-allscores);

allscores = allscores(i);
allcodes  = allcodes(i,:);
Lett = 'ACGU';
allseqs   = Lett(allcodes);

for i = 1:100,
  plot(0.1,allscores(i),'r*');
  text(0.1,allscores(i),allseqs(i,:),'fontsize',6);
  fprintf('%4d %s %16.12f\n', i, allseqs(i,:), allscores(i));
%  plot(0.1,allscores(end-i),'r*');
%  text(0.1,allscores(end-i),allseqs(end-i,:));
end

% score motif according to basepairs and BPh


% score motif according to pairs, BPh, and stacking






break
    
clear Text

r = 1;
while r <= L/2,
  m = min(Cand(2*r,[1 2]));                 % the smaller 
  M = max(Cand(2*r,[1 2]));
  i = m:M;

  mm = min(Cand(2*r,[3 4]));
  MM = max(Cand(2*r,[3 4]));
  j = mm:MM;

  f = Cand(2*r,N+1);

  [a,b,c] = find(File(f).Edge(i,j));
  c = abs(c);
  c = c .* (c > 0) .* (c < 13);             % "true" basepairs
  c = c(find(c));
  c = sort(c);

  if length(c) == 2,
    [a,b,c] = find(File(f).Edge(i,j));
    c = abs(c);
    c = c .* (c > 100) .* (c < 113);       % "near" basepairs
    c = c(find(c));
    c = sort(c);
    c = [1; 1; c];
  end

  if length(c) > 0,
    Text{r} = '';
    for g = 1:length(c),
      Text{r} = [Text{r} sprintf(' %s', zEdgeText(c(g)))];
    end
  else
    Text{r} = ' No interaction';
  end

%  fprintf('%s\n',Text{r});


  if Verbose > 1,
    zShowInteractionTable(File(f),[mm:MM m:M]);

    VP.Sugar = 1;
    clf
    zDisplayNT(File(f), [mm:MM m:M],VP)
    disp(Text{r})
    drawnow
    pause
  end

  r = r + 1;

end

[t,p] = sort(Text);
Text = Text(p);

[b,i,j] = unique(Text);

% i points to the unique rows of Text
% j maps 1:L to 1:length(b)

fprintf('Of the %d internal loops identified, there were %d unique sequence signatures.\n', L, length(i));

for r = 1:length(i),
  fprintf('%5d instances of %s\n', length(find(j==r)), b{r});
end
