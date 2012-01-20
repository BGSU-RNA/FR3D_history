% xFindMotifsFromLibrary finds all instances of motifs from the motif
% library in the folder MotifLibrary in the given file(s)
 
%function [File] = xFindMotifsFromLibrary(Files)

Motif = dir(['MotifLibrary' filesep 'LIB*.mat']);
% Motif = [Motif; dir(['MotifLibrary' filesep 'LIB*.MAT'])];         
% for the Mac, but then you need to eliminate redundant listings for the PC

if strcmp(class(Files),'char'),
  Filename = Files;
  Files = zAddNTData(Filename,2);
end

UsingLibrary = 1;                                % so FR3D just searches

for m = 1:length(Motif),
  load(['MotifLibrary' filesep Motif(m).name]);
  MSearch(m) = Search;
end


for f = 1:length(Files),
  File = Files(f);
  Filenames = {File.Filename};

  if ~(exist(['MotifList' filesep File.Filename]) == 7),        % if directory doesn't yet exist
    mkdir(['MotifList' filesep File.Filename]);
  end

  ffid = fopen(['MotifList' filesep File.Filename filesep 'index.html'],'w');
  fprintf(ffid,'<html>\n<title>%s motif list\n</title>\n', File.Filename);
  fprintf(ffid,'<body>\n');
  fprintf(ffid,'<h1>Motif list for %s</h1>\n', File.Filename);
  fprintf(ffid,'<a href="interactions.html">List of pairwise interactions found in %s</a><br>', File.Filename);
  fprintf(ffid,'<a href="../">Return to list of PDB files</a><br>');

  for m = 1:length(Motif),
    MotifNumber = Motif(m).name(1:8);
    Search = MSearch(m);

%    fprintf('\nSearching %s for %s\n', File.Filename, Search.Query.Description);

    Query = Search.Query;

    clear Search

    FR3D

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