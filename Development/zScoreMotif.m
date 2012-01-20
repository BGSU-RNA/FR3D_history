% zScoreMotif takes a set of nucleotides from a FR3D search, tallies up the observed frequency of the different sequences, and scores them according to various pairwise scoring schemes


% load 2008-06-06_13_44_49-Internal_loops_flanked_by_nested_cWW_pairs.mat
% load 2010-05-18_15_42_07-Sarcin_7_mixed_2aw4.mat

method = 2;                                  % pIsoScore method to use

Verbose = 1;

Filenames = Search.Filenames;

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData(Filenames,0,[],Verbose);   % load PDB data
  File = zAttachAlignment(File);              % attach alignment data
else
  [File,SIndex] = zAddNTData(Filenames,0,File,Verbose); %add PDB data if needed
end                       % SIndex tells which elements of File to search

File = File(SIndex);

Cand = Search.Candidates;

[L,N] = size(Cand);        % number of candidates
N = N - 1;                 % number of nucleotides

% ------------------------------------------- collect all sequence variants

AllLett = [];
for c = 1:L,
  Lett = [];
  for n = 1:N,
    Lett = [Lett File.NT(Cand(c,n)).FASTA];
  end
  AllLett = [AllLett; Lett];
end

% ------------------------------------------- remove non-ACGU letters

i = find(sum(AllLett == 'N',2) == 0);
AllLett = AllLett(i,:);

i = find(sum(AllLett == '-',2) == 0);
AllLett = AllLett(i,:);

% ------------------------------------------- find unique sequences

[seqs,counts] = zUnique(AllLett)

S = length(seqs(:,1));                      % number of unique sequences

codes = 1*(seqs=='A') + 2*(seqs=='C') + 3*(seqs=='G') + 4*(seqs=='U');

% --------------------------------------- find consensus interaction list

[Edge,BPh] = pConsensusInteractions(Search);

% ------------------------------------------- generate all possible sequences

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

% ------------------------------------------ score all sequences

best = codes(1,:);                      % most common sequence

[score,ScoringMethodName]     = zScoreMotifSequences(codes,best,Edge,BPh,method);
[allscores,ScoringMethodName] = zScoreMotifSequences(allcodes,best,Edge,BPh,method);


% ------------------------------------------ display scores and sequences

Lett = 'ACGU';

for k = 1:length(allscores(1,:)),          % loop through scoring methods

  figure(k)
  clf

  loglog(counts,score(:,k),'*')
  hold on
  for s = 1:S,
    text(counts(s)*1.1,score(s,k),seqs(s,:),'fontsize',6);
  end

  [y,i] = sort(-allscores(:,k));           % sort by decreasing score
  allscores = allscores(i,:);              % re-order scores
  allcodes  = allcodes(i,:);               % re-order codes

  M = 100;                                   % number of sequences to display
  i = 1:(M + S);

  d = zDistance(allcodes(i,:),codes);
  i = find(min(d,[],2) > 0);                 % unobserved codes

  i = i(1:M);

  fprintf('Scoring method: %s\n', ScoringMethodName{k});

  for j = 1:M,
    plot(0.1,allscores(i(j),k),'r*');
    text(0.1*1.1,allscores(i(j),k),Lett(allcodes(i(j),:)),'fontsize',6);
    fprintf('%4d %s %16.12f\n', i(j), Lett(allcodes(i(j),:)), allscores(i(j),k));
%    plot(0.1,allscores(end-i(j)),'r*');
%    text(0.1,allscores(end-i(j)),allseqs(end-i(j),:));
  end

  title(['Motifs scored by ' ScoringMethodName{k}]);

end

% score motif according to basepairs and BPh


% score motif according to pairs, BPh, and stacking



