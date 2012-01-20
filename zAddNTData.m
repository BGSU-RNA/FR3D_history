% zAddNTData(Filenames,ReadCode,File) reads RNA structure data files, if
% necessary, so that all molecules listed in Filenames are present in File.
% The parameter Index is a non-redundant list of indices of File
% corresponding to names in Filenames.

function [File,Index] = zAddNTData(Filenames,ReadCode,File)

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};                % make into a cell array
end

% ----------------------------------------- Read PDB lists, if any

FullList = [];

for j=1:length(Filenames),
  FullList = [FullList; zReadPDBList(Filenames{j})];
end

if length(FullList) > 0,

for f = 1:length(FullList),
  FL{f} = lower(FullList{f});
  FirstOccurence(f) = isempty(strmatch(FL{f},FL(1:(f-1)),'exact'));
end

FullList = FullList(find(FirstOccurence));

if nargin < 2,
  ReadCode = 0;
end

if nargin < 3,                            % no data already loaded
  File = zGetNTData(FullList{1},ReadCode);
  FullList = FullList(2:end);             % remove first filename
  Index(1) = 1;                           
end

for j = 1:length(File),
  LoadedFiles{j} = lower(File(j).Filename);
end

F = length(File);

for f = 1:length(FullList),                       % loop through PDB list
  i = strmatch(lower(FullList{f}), LoadedFiles, 'exact');
  if isempty(i),                                  % if PDB not loaded,
    File(F+1) = zGetNTData(FullList{f},ReadCode); %   load it
    F = length(File);
    LoadedFiles = [LoadedFiles FullList{f}];
    Index(f) = F;                                 %   point to it
  else                                            % but if PDB has been loaded
    Index(f) = i(1);                              %   point to first instance
    if ReadCode > 0,                              % reanalyze
      File(Index(f)) = zGetNTData(File(Index(f)).Filename,ReadCode);
    end
  end
end

else

fprintf('No files specified to read in %s\n', Filenames{1});

end
