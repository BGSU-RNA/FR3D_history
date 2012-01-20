% xAnnotateWithKnownMotifs(File) finds all instances of motifs from the motif
% library in the folder MotifLibrary in the given file(s)
 
function [File] = xAnnotateWithKnownMotifs(Files,Verbose)

if nargin < 2,
  Verbose = 0;
end

Motif = dir([pwd filesep 'MotifLibrary' filesep 'LIB*.mat']);

if strcmp(class(Files),'char'),
  Filename = Files;
  Files = zAddNTData(Filename,0,[],Verbose);
end

UsingLibrary = 1;                                % so FR3D just searches

for f = 1:length(Files),
  File = Files(f);
  FN   = upper(File.Filename);
  Filenames = {File.Filename};

  LText{1} = ['<a href = "index.html">Return to FR3D home page for ' FN '</a><br>'];
  LText{2} = ['<a href = "' FN '_interactions.html">List of all pairwise interactions in ' FN '</a><br>'];
  LText{3} = ['<a href = "' FN '_basepairs.html">List of basepair interactions in ' FN '</a><br>'];
  LText{4} = ['<a href = "' FN '_stacking.html">List of stacking interactions in ' FN '</a><br>'];
%  LText{4} = ['<a href = "' FN '_basephosphate.html">List of base-phosphate interactions in ' FN '</a><br>'];
  LText{5} = ['<a href = "' FN '_motifs.html">List of motifs found in ' FN '</a><br>'];
  LText{6} = ['<a href="http://www.rcsb.org/pdb/explore/explore.do?structureId=' FN '">PDB entry for ' FN '</a><br>'];
  LText{7} = ['<a href="../">Return to list of analyzed structures</a><br>'];
  LText{8} = ['<a href="../../basepairs">Basepair catalog</a><br>'];
  LText{9} = ['<a href="../../MotifLibrary/index.html">FR3D motif library</a><br>'];
  LText{10} = ['<a href="../../index.html">FR3D home page</a><br>'];
  LText{11} = ['<a href="http://rna.bgsu.edu">BGSU RNA group home page</a><br><br>'];

  DN = [pwd filesep 'Web' filesep 'AnalyzedStructures' filesep FN];

  if ~(exist(DN) == 7),        % if directory doesn't yet exist
    mkdir(DN);
  end

  ffid = fopen([DN filesep FN '_motifs.html'],'w');

  fprintf(ffid,'<html>\n<title>%s motif list\n</title>\n', FN);
  fprintf(ffid,'<body>\n');

  fprintf(ffid,'<h1>Motif list for %s</h1>\n', FN);

  for L = 1:length(LText),
    if L ~= 5,
      fprintf(ffid,'%s\n', LText{L});
    end
  end

  for m = 1:length(Motif),

   if (Motif(m).name(9) == '_') && (Motif(m).name(4) ~= 'U'),

    load(['MotifLibrary' filesep Motif(m).name]);

    MotifNumber = Motif(m).name(1:8);
    MotifName   = Motif(m).name(10:end);

    if Verbose > 0,
      fprintf('\nSearching %s for %s\n', FN, Motif(m).name);
    end

    Query = Search.Query;

    clear Search

    xFR3DSearch

    [s,t] = size(Candidates);

    if s > 0,
      fprintf(ffid,'<a name=%s>\n',MotifNumber);
      fprintf(ffid,'<h2><a href="%s">%s</a></h2>\n',['http://rna.bgsu.edu/FR3D/MotifLibrary/' MotifNumber '/index.html'], strrep(Motif(m).name,'.mat',''));

      [y,i] = sort(Candidates(:,1));               % sort by 1st nucleotide #
      Search.Candidates = Candidates(i,:);         % re-order candidates
      Search.Discrepancy = Search.Discrepancy(i);

      fprintf(ffid,'<pre>\n');
      Text = xListCandidates(Search,Inf,6,{MotifNumber});
      for c = 4:length(Text),
        fprintf(ffid,'%s\n',Text{c});
      end
      fprintf(ffid,'</pre>\n');
    end
   end
  end
  fclose(ffid);
end
