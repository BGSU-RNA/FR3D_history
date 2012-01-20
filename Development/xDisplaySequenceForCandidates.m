
g = 2;

switch g
case 1,
  load LIB00005_IL_tSH-tHH-cSH-tWH-tHS_sarcin-ricin.mat
  NofI = 1:13;
case 2,
  load 2009-01-28_16_32_47-GC_1BPh_swap.mat
  NofI = [1 2];                     % nucleotides of interest for comparison
end

if exist('File'),

else
  [File,FIndex] = zAddNTData(Search.CandidateFilenames,2);
  File = File(FIndex);
  File = zAttachAlignment(File,1);          
end

Cand = Search.Candidates;


[L,N] = size(Cand);
N = N - 1;

for c = 1:L,
 f = Search.Candidates(c,N+1);               % file number
 s = [];

 if isfield(File(f).NT(Cand(c,1)),'FASTA'),

  for j = 1:N,
    s(j) = length(File(f).NT(Cand(c,j)).FASTA);
  end
  if min(s) == max(s) && min(s) > 0,
    seq = [];
    for j = NofI,
      seq = [seq File(f).NT(Cand(c,j)).FASTA];
    end

    
    fprintf('Candidate %d is',c);
    for j = 1:N,
      fprintf(' %s%4s', File(f).NT(Cand(c,j)).Base,File(f).NT(Cand(c,j)).Number);
    end
    fprintf(' from %s chain %s\n', File(f).Filename, File(f).NT(Cand(c,1)).Chain);
    fprintf('Candidate %d has relevant sequence %s\n', c, cat(2,File(f).NT(Cand(c,NofI)).Base));
    fprintf('Summary of sequence variants found in alignment:\n');



    [b,t] = zUniqueRows(seq);

    L = min(10,length(t));

    for r = 1:L,
      fprintf('Sequence %s occurs %4d times\n', b(r,:), t(r));
    end

    if L < length(t),
      fprintf('The remaining variants account for %8.2f%% of the sequences\n', 100*(1-sum(t(1:L))/sum(t)));
    end

    fprintf('\n');
  end
 end
end
