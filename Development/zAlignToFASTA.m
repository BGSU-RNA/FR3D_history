% zAlignToFASTA determines a vector FastaCol which maps the indices in File corresponding to Chain to the corresponding column numbers in the Aligned field of the data in FASTA.  If Entry is specified or non-zero, zAlignToFASTA uses the indicated entry of the FASTA file as the key.  Otherwise, it tries aligning to each sequence successively until a good match is found.

% The presumption is that an alignment contains information from only one chain

% File = zAddNTData('2avy');
% FASTA = zReadFASTA('Alignments\16S_Bacterial_Stombaugh_et_al_Sup_Mat_S2.fasta');
% Entry = 1;

% F = zAlignToFASTA(File,'A',FASTA,28,2);


function [File] = zAlignToFASTA(File,Chain,FASTA,Entry,Verbose)

if nargin < 5,
  Verbose = 0;
end

if nargin < 4,
  Entry = 0;
elseif Entry == 0,
  Verbose = 2;
end

C = cat(2,File.NT.Chain);
a = find(C == Chain);

StructureSequence = cat(2,File.NT(a).Base);

if Entry == 0,
  if Verbose > 0,
    fprintf('Aligning sequence from structure to each entry in FASTA file\n');
  end
  for e = 1:length(FASTA),
    p = 0.99;
    d = 2;
%   [matches,align1,align2,s1,s2] = dNeedlemanWunsch(StructureSequence,FASTA(e).Sequence,p,d);

    s = FASTA(e).Sequence;
    s = strrep(s,'?','N');

    [matches,align1,align2,s1,s2] = zNeedlemanWunsch(StructureSequence,s);

    if Verbose > 1,
      fprintf('Sequence %4d has %4d matches and characters %8s; from %s\n', e, matches, unique(s), FASTA(e).Header);
    end

    m(e) = matches;
  end

  [y,Entry] = max(m);

  if Verbose > 1,
    [y,i] = sort(-m);
    fprintf('Top choices:\n');
    for j = 1:10,
      fprintf('Sequence %4d has %4d matches and characters %8s; from %s\n', i(j), m(i(j)), unique(FASTA(i(j)).Sequence), FASTA(i(j)).Header);
    end
  end      

  if Verbose > 0,
    fprintf('Using sequence %d, which has %d matches out of %d bases.\n', Entry, m(Entry), length(StructureSequence));
  end

  [n,t,r] = xlsread('Alignments\StructureToAlignmentMap.xls');
  for i = 1:length(r(:,1)),
    if strcmp(upper(File.Filename),upper(r{i,1})) && (Chain == r{i,2}),
      r{i,4} = Entry;
    end
  end
  xlswrite('Alignments\StructureToAlignmentMap.xls',r);
  if Verbose > 0,
    fprintf('Added this choice to structure to alignment map.\n');
  end
end

  p = 0.99;
  d = 2;
%  [matches,align1,align2,s1,s2] = dNeedlemanWunsch(StructureSequence,FASTA(Entry).Sequence,p,d);

  [matches,align1,align2,s1,s2] = zNeedlemanWunsch(StructureSequence,FASTA(Entry).Sequence);


  k = double(FASTA(Entry).Aligned ~= '-');  % 1 if not a gap
  m = cumsum(k);
  
  x = zInvertFunction(align1);              % map 
  y = [align2 (length(FASTA(Entry).Sequence)+1)]; 
  z = [zInvertFunction(m) length(FASTA(Entry).Aligned)]; % maps elements of seq to columns

%align1

%length(x)
%length(y)
%length(z)


  FastaCol = z(y(x));  

  Alignment = cat(1,FASTA.Aligned);

%FastaCol
%min(FastaCol)
%max(FastaCol)

  for i = a,
%[i i-a(1)+1]
    File.NT(i).FASTA = Alignment(:,FastaCol(i-a(1)+1));
    File.NT(i).FASTACol = FastaCol(i-a(1)+1);
  end

  if Verbose > 1,
    fprintf('Structure: %s\n', StructureSequence);
    fprintf('Alignment: %s\n', FASTA(Entry).Aligned(FastaCol));

    clear c d
    for i = 1:length(File.NT),
      if length(File.NT(i).FASTA > 0),
        c(1,i) = File.NT(i).Base;
        d(1,i) = File.NT(i).FASTA(Entry);
      end
    end
    fprintf('\n');
    fprintf('Structure: %s\n', c);
    fprintf('Alignment: %s\n', d);
  end

return



% check:

F = File(2);



for i = 1:length(FASTA),
  fprintf('%d %s\n', i, FASTA(i).Header);
end