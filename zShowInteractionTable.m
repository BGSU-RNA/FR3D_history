% zShowInteractionTable(File,NTList) displays a table of interactions
% among the bases listed in Indices

function [void] = zShowInteractionTable(File,NTList)

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(NTList),'char'),
  NTList = {NTList};
end

if strcmp(class(NTList),'cell'),
  Indices = zIndexLookup(File,NTList);
else
  Indices = NTList;
end

fprintf('      ');
for j=2:length(Indices),
  fprintf('%6s',[File.NT(Indices(j)).Base File.NT(Indices(j)).Number]);
end
fprintf('\n');

for i=1:(length(Indices)-1),
  fprintf('%6s',[File.NT(Indices(i)).Base File.NT(Indices(i)).Number]);
  for j=2:length(Indices),
    if j > i,
      fprintf('%6.1f', File.Inter(Indices(i),Indices(j)))
    else
      fprintf('      ');
    end
  end
  fprintf('\n');
end
   