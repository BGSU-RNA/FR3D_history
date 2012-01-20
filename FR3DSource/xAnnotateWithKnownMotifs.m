% xAnnotateWithKnownMotifs(File) finds all instances of motifs from the motif
% library in the folder MotifLibrary in the given file(s)
 
function [File] = xAnnotateWithKnownMotifs(Files)

Motif = dir([pwd filesep 'MotifLibrary' filesep 'LIB*.mat']);
% Motif = [Motif; dir(['MotifLibrary' filesep 'LIB*.MAT'])];         
% for the Mac, but then you need to eliminate redundant listings for the PC

if strcmp(class(Files),'char'),
  Filename = Files;
  Files = zAddNTData(Filename,2);
end

UsingLibrary = 1;                                % so FR3D just searches

%for m = 1:length(Motif),
%  load(['MotifLibrary' filesep Motif(m).name]);
%  MSearch(m) = Search;
%end


for f = 1:length(Files),
  File = Files(f);
  FN = upper(File.Filename);
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

  ffid = fopen([pwd filesep 'Web' filesep 'AnalyzedStructures' filesep FN filesep FN '_motifs.html'],'w');

  fprintf(ffid,'<html>\n<title>%s motif list\n</title>\n', FN);
  fprintf(ffid,'<body>\n');


  fprintf(ffid,'<h1>Motif list for %s</h1>\n', FN);

  for L = 1:length(LText),
    if L ~= 5,
      fprintf(ffid,'%s\n', LText{L});
    end
  end

  for m = 1:length(Motif),

    load(['MotifLibrary' filesep Motif(m).name]);

    MotifNumber = Motif(m).name(1:8);

%    fprintf('\nSearching %s for %s\n', FN, Search.Query.Description);

    Query = Search.Query;

    clear Search

    xFR3DSearch

    [s,t] = size(Candidates);

    fprintf(ffid,'<a name=%s>\n',MotifNumber);
    fprintf(ffid,'<a href="%s"><h2>%s</h2></a>\n',['../../MotifLibrary/' MotifNumber '/index.html'], Query.Description);

    if s > 0,
      [y,i] = sort(Candidates(:,1));               % sort by 1st nucleotide #
      Search.Candidates = Candidates(i,:);         % re-order candidates
      Search.Discrepancy = Search.Discrepancy(i);

      fprintf(ffid,'<pre>\n');
      Text = xListCandidates(Search,Inf,6,{MotifNumber});
      for c = 1:length(Text),
        fprintf(ffid,'%s\n',Text{c});
      end
      fprintf(ffid,'</pre>\n');
    else
      fprintf(ffid,'No motifs of this type found\n');
    end
  end

  fclose(ffid);
end
