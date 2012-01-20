
StartFile = '299D';

[NamesLists,Names] = mGetPDBFilenames;

s = 1;

for i=1:length(Names),
  if strcmp(StartFile,Names{i}),
    s = i;
  end
end


for f = s:length(Names),
  File = zAddNTData(Names{f},2);            % load RNA data
  zAnalyzedFilesHTML(File);                 % make HTML
  xAnnotateWithKnownMotifs(File);           % find and list motifs
end
