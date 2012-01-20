% xWriteHTMLForMotifs writes an html file for each motif in the motif library
% and lists a centroid and instances

function [void] = xWriteHTMLForMotifs

Motif = dir(['MotifLibrary' filesep 'LIB*.mat']);
% Motif = [Motif; dir(['MotifLibrary' filesep 'LIB*.MAT'])];         
% for the Mac, but then you need to eliminate redundant listings for the PC

% ----------- open index file for motif library

ifid = fopen(['MotifLibrary' filesep 'index.html'],'w');
fprintf(ifid,'<html>\n<title>NDB/FR3D motif library</title>\n');
fprintf(ifid,'<body>\n');
fprintf(ifid,'<h1>NDB/FR3D motif library</h1>\n');
fprintf(ifid,'<a href="../MotifList/">Return to list of PDB files</a><br>');
for m = 1:length(Motif),
  load(['MotifLibrary' filesep Motif(m).name]);
  MotifNumber = Motif(m).name(1:8);

  fprintf('Writing results for motif %s: %s\n', MotifNumber, Search.Query.Description);
  % --------------------------- Add to index file for motif library

  fprintf(ifid,'<a href="%s">%s: %s</a><br>\n', [MotifNumber filesep 'index.html'], MotifNumber, Search.Query.Description);

  % --------------------------- Write index file and instances for this motif




  if ~(exist(['MotifLibrary' filesep MotifNumber]) == 7), % if directory doesn't yet exist
    mkdir(['MotifLibrary' filesep MotifNumber]);
  end

  fid = fopen(['MotifLibrary' filesep MotifNumber filesep 'index.html'],'w');
  fprintf(fid,'<html>\n<title>%s motif: %s\n</title>\n', MotifNumber, Search.Query.Description);
  fprintf(fid,'<body>\n');
  fprintf(fid,'<h1>Motif %s: %s</h1>\n', MotifNumber, Search.Query.Description);

  fprintf(fid,'<a href="%s">List of instances</a><br>\n',['instances.html']);

  fprintf(fid,'<a href="%s">FR3D search parameters</a><br>\n',['FR3D_search_params.html']);

  fprintf(fid,'<a href="%s">Alignment and statistical information</a><br>\n',['alignment.html']);

  fprintf(fid,'<img src = "%s">\n',['annotation.png']);

  fprintf(fid,'<img src = "%s">\n',['centroid.png']);

  fprintf(fid,'<a href="%s">Return to motif library index</a><br>\n',['../index.html']);

  fclose(fid);

  % ----------------------------- write list of instances of this motif

  fid = fopen(['MotifLibrary' filesep MotifNumber filesep 'instances.html'],'w');
  fprintf(fid,'<h1>Instances of %s: %s</h1>\n',MotifNumber,Search.Query.Description);

  if length(Search.Candidates(:,1)) > 0,
    fprintf(fid,'<pre>\n');
    Text = xListCandidates(Search,Inf,6,{MotifNumber});
    for c = 1:length(Text),
      fprintf(fid,'%s\n',Text{c});
    end
    fprintf(fid,'</pre>\n');
  else
    fprintf(fid,'No motifs of this type found\n');
  end

  fclose(fid);

  % ----------------------------- write alignment of instances of this motif

  fid = fopen(['MotifLibrary' filesep MotifNumber filesep 'alignment.html'],'w');
  fprintf(fid,'<h1>Alignment of instances of %s: %s</h1>\n',MotifNumber,Search.Query.Description);

  if length(Search.Candidates(:,1)) > 0,
    fprintf(fid,'<pre>\n');
    Text = xAlignCandidates(Search.File,Search,1,{MotifNumber});
    for c = 1:length(Text),
      fprintf(fid,'%s\n',Text{c});
    end
    fprintf(fid,'</pre>\n');
  else
    fprintf(fid,'No motifs of this type found\n');
  end

  fclose(fid);

end

fclose(ifid);