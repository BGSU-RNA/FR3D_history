
WhatToWrite = [1 1 1 1];           % which pieces to write right now
WhatToWrite = [0 1 0 0];           % which pieces to write right now

% WhatToWrite(1)  - html pages
% WhatToWrite(2)  - circular diagram
% WhatToWrite(3)  - pdb file
% WhatToWrite(4)  - annotate with known motifs


load PDBInfo

Names = {'2AW4','2J01','1S72','2AVY','1J5E','3BWP','1ZZN'};
Names = {'1HR2'};
Names = {'3BWP'};
Names = zReadPDBList('Nonredundant_2008_02_21_list');
Names = t(:,1);                        % names of files from PDB/NDB

%for f = 415:length(Names),  
for f = 1:415,  
 t = cputime;

 File = zAddNTData(Names{f},0,[],2);              % load RNA data

 t(end+1) = cputime;

 if ~isempty(File.NT),

  fprintf('Writing HTML files for %s, file %d of %d\n', Names{f}, f, length(Names));

  if WhatToWrite(1) > 0,
    zAnalyzedFilesHTML(File);                        % make HTML
  end

  t(end+1) = cputime;

  % ----------------------------------------------- Create circular diagram

  if WhatToWrite(2) > 0,
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

  end

  t(end+1) = cputime;

  % ----------------------------------------------- Write PDB file

  if WhatToWrite(3) > 0,
  zWritePDB(File,[mypath File.Filename '_RNA.pdb']);
  end
  t(end+1) = cputime;

  % ----------------------------------------------- Annotate with known motifs

  if WhatToWrite(4) > 0,
  xAnnotateWithKnownMotifs(File,1,1);           % find and list motifs
  end

  t(end+1) = cputime;

  fprintf('Time taken:');
  fprintf(' %6.2f', diff(t));
  fprintf(' seconds \n');

 end
end

break

% ------------------------------------------------- Write index files

DN = [pwd filesep 'Web' filesep 'AnalyzedStructures' filesep FN];

if ~(exist(DN) == 7),        % if directory doesn't yet exist
  mkdir(DN);
end




fid = fopen([DN 'short_list.html'],'w');

