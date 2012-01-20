
StartFile = '409d';
%StartFile = '';

[NamesLists,Names] = mGetPDBFilenames;

%Names = zReadPDBList('NonRedundant_2008_02_21_list');

s = find(ismember(upper(Names),upper(StartFile)));

if isempty(s)
  s = 1;
end

for f = s(1):length(Names),  
  File = zAddNTData(Names{f},0,[],2);              % load RNA data

  zAnalyzedFilesHTML(File);                 % make HTML

  % ----------------------------------------------- Create circular diagram

  mypath = [pwd filesep 'Web' filesep 'AnalyzedStructures'];
  mypath = [mypath filesep File.Filename filesep];

  clf
  zCircularDiagram(File,1);
  saveas(gcf,[mypath File.Filename '_circular_diagram.png'],'png');
  [X,map] = imread([mypath File.Filename '_circular_diagram.png']);
  Y = X(30:830,210:1030,:);
  imwrite(Y,[mypath File.Filename '_circular_diagram.png']);

  clf
  zCircularDiagram(File,0.1);
  saveas(gcf,[mypath File.Filename '_circular_diagram.pdf'],'pdf');

  % ----------------------------------------------- Write PDB file

  zWritePDB(File,[mypath File.Filename '_RNA.pdb']);

  % ----------------------------------------------- Annotate with known motifs

  xAnnotateWithKnownMotifs(File);           % find and list motifs
end
