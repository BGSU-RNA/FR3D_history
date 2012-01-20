% zScoreMotif takes a set of nucleotides from a FR3D search, tallies up the observed frequency of the different sequences, and scores them according to various pairwise scoring schemes

method = 2;                                  % pIsoScore method to use
Verbose = 1;

if ~exist('Stacks'),
  load 2010-05-30_10_07_51-s35_NR_4A.mat
  Stacks{1} = Search;
  load 2010-05-29_10_24_40-s33_NR_4A_exclude_redundant.mat
  Stacks{2} = Search;
  load 2010-05-29_10_39_07-s55_NR_4A.mat
  Stacks{3} = Search;
end

for motif = 1:1,

switch motif
case 0,                                      % don't load a new motif

case 1,
  load 2010-06-07_16_14_35-Sarcin_7_mixed_2qan_2qbe.mat
  MotifName = 'Sarcin 7-nucleotide';
case 2,
  load 2010-05-22_23_06_50-Sarcin_9_mixed_2aw4_2avy.mat
  MotifName = 'Sarcin 9-nucleotide';
case 3,
%  load 2010-05-22_23_38_17-C_loop_core_with_flanking_cWW_centroid_tuned_2aw4_2avy.mat
  load 2010-05-23_10_17_08-C_loop_core_with_flanking_cWW_2aw4_2avy_ACAUAU.mat
  MotifName = 'C-loop';
case 4,
  load 2010-05-22_23_54_01-Kink_turn_highly_constrained.mat
  MotifName = 'Kink-turn';
case 5,
  load 2010-05-25_12_25_48-IL_tSH_tWH_tHS_2aw4_2avy
  Cand = Search.Candidates
  Cand = Cand(:,[2 3 4 7 8 9 11]);      % remove flanking cWW pairs
  Search.Candidates = Cand;
  MotifName = 'tSH-tWH-tHS IL';
case 6,
  load 2010-05-23_00_28_00-Helix_10_nucleotides_2avy_2aw4.mat
  MotifName = '8-nucleotide helix';
end




Cand = Search.Candidates;

[L,N] = size(Cand);        % number of candidates
N = N - 1;                 % number of nucleotides

Filenames = Search.Filenames;

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData(Filenames,0,[],Verbose);   % load PDB data
else
  [File,SIndex] = zAddNTData(Filenames,0,File,Verbose); %add PDB data if needed
end                       % SIndex tells which elements of File to search

File = File(SIndex);

if ~isfield(File(1).NT(1),'FASTA')
  File = zAttachAlignment(File);              % attach alignment data
end

% ------------------------------------------- collect all sequence variants

AllLett = [];
clear Lett
for c = 1:L,                                % loop through candidates
  f = Cand(c,N+1);                          % file number
  Let = [];
  for n = 1:N,
    Let = [Let File(f).NT(Cand(c,n)).FASTA];
  end

  AllLett = [AllLett; Let];

  [seqs,counts] = zUnique(Let);

  codes = 1*(seqs=='A') + 2*(seqs=='C') + 3*(seqs=='G') + 4*(seqs=='U');
  i = find(min(codes,[],2) > 0);

  Let = seqs(i,:);

  Counts{c} = counts(i);
  Lett{c} = Let;

  if Verbose > 0,
    fprintf('File %s instance ', File(f).Filename);
    for n = 1:N,
      fprintf('%s%s ', File(f).NT(Cand(c,n)).Base, File(f).NT(Cand(c,n)).Number);
    end
    fprintf('\n');

    for i = 1:length(Counts{c}),
      fprintf('  %s occurs %4d times\n', Lett{c}(i,:), Counts{c}(i));
    end

  end
end

% ------------------------------------------- find unique sequences

[seqs,counts] = zUnique(AllLett);
codes = 1*(seqs=='A') + 2*(seqs=='C') + 3*(seqs=='G') + 4*(seqs=='U');

% ------------------------------------------- remove non-ACGU letters

i = find(min(codes,[],2) > 0);
seqs = seqs(i,:);
codes = codes(i,:);
counts = counts(i);
S = length(seqs(:,1));                       % number of unique sequences

seqs
counts

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
      c(a) = 1;                              % roll this one over
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

if 0 > 1,
  NRList = 'Nonredundant_2009-05-14_list';
  if ~exist('AllFile'),                           % if no molecule data is loaded,
    [AllFile,SIndex] = zAddNTData(NRList,0,[],Verbose);   % load PDB data
  else
    [AllFile,SIndex] = zAddNTData(NRList,0,AllFile,Verbose); %add PDB data if needed
  end                       % SIndex tells which elements of File to search
end

f = Search.Candidates(1,N+1);
clear Motif
for i = 1:N,
  Motif.NT(i) = Search.File(f).NT(Cand(1,i));
end

best = codes(1,:);                      % most common sequence

[allscores,ScoringMethodName] = zScoreMotifSequences(allcodes,Search,Edge,BPh,method,Stacks);

[score,ScoringMethodName] = zScoreMotifSequences(codes,Search,Edge,BPh,method,Stacks);

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

  M = 100;                                 % number of sequences to display

  fprintf('Scoring method: %s\n', ScoringMethodName{k});

  g = 0;                                  % number of observed sequences found
  a = 1;                                  % current sequence
  while g < S && a < M,                   % go until all observed are listed
    fprintf('Sequence %4d %s Score %16.16f ', a, Lett(allcodes(a,:)), allscores(a,k));
    y = find(zDistance(allcodes(a,:),codes) == 0);
    if ~isempty(y),
      g = g + 1;
      fprintf(' observed %5d times\n', counts(y));
    else
      fprintf('\n');
    end
    a = a + 1;
  end

  i = 1:(M + S);
  d = zDistance(allcodes(i,:),codes);
  i = find(min(d,[],2) > 0);               % unobserved codes to plot
  i = i(1:M);

  for j = 1:M,
    plot(0.1*1.02,allscores(i(j),k),'r*');
    text(0.1*1.12,allscores(i(j),k),Lett(allcodes(i(j),:)),'fontsize',6);
%    plot(0.1,allscores(end-i(j)),'r*');
%    text(0.1,allscores(end-i(j)),allseqs(end-i(j),:));
  end

  title(MotifName);
  xlabel(['Motifs scored by ' ScoringMethodName{k}]);

  orient tall

  saveas(gcf,['Motif scoring\' MotifName '_BPhGeom_Depth' num2str(k) '_BP' num2str(method) '.pdf'],'pdf');
  saveas(gcf,['Motif scoring\' MotifName '_BPhGeom_Depth' num2str(k) '_BP' num2str(method) '.png'],'png');


end


end
