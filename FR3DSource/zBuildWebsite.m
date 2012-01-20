
load PDBInfo

Names = t(:,1);                        % names of files from PDB/NDB

%for f = 1:length(Names),  

%1354

for f = 1391:length(Names),  
 t = cputime;

 File = zAddNTData(Names{f},0,[],2);              % load RNA data

 t(end+1) = cputime;

 if ~isempty(File.NT),

  fprintf('Writing HTML files for %s, file %d of %d\n', Names{f}, f, length(Names));

  zAnalyzedFilesHTML(File);                        % make HTML

  t(end+1) = cputime;

  % ----------------------------------------------- Create circular diagram

  mypath = [pwd filesep 'Web' filesep 'AnalyzedStructures'];
  mypath = [mypath filesep File.Filename filesep];

  clf
  zCircularDiagram(File,1);
  saveas(gcf,[mypath File.Filename '_circular_diagram.png'],'png');
  [X,map] = imread([mypath File.Filename '_circular_diagram.png']);
  Y = X(30:830,210:1030,:);
  imwrite(Y,[mypath File.Filename '_circular_diagram.png']);

  t(end+1) = cputime;

  clf
  zCircularDiagram(File,0.1);
  saveas(gcf,[mypath File.Filename '_circular_diagram.pdf'],'pdf');

  t(end+1) = cputime;

  % ----------------------------------------------- Write PDB file

  zWritePDB(File,[mypath File.Filename '_RNA.pdb']);

  t(end+1) = cputime;

  % ----------------------------------------------- Annotate with known motifs

  xAnnotateWithKnownMotifs(File);           % find and list motifs

  t(end+1) = cputime;

  fprintf('Time taken:');
  fprintf(' %6.2f', diff(t));
  fprintf(' seconds \n');
 end
end
