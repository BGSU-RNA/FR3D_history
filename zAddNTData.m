% zAddNTData(Filenames,ReadCode,File) reads RNA structure data files, if
% necessary, so that all molecules listed in Filenames are present in File.
% The parameter Index is a non-redundant list of indices of File
% corresponding to names in Filenames.

function [File,Index] = zAddNTData(Filenames,SizeCode,File)

if nargin < 2,
  SizeCode = 1;                           % default is to read full files
end

if SizeCode == 0,
  SizeCode = 1;
end

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};                % make into a cell array
end

% ----------------------------------------- Read PDB lists, if any

FullList = [];

for j=1:length(Filenames),
  FullList = [FullList; zReadPDBList(Filenames{j})];
end

% ----------------------------------------- Read PDB files

if length(FullList) > 0,

for f = 1:length(FullList),
  FL{f} = lower(FullList{f});
  FirstOccurence(f) = isempty(strmatch(FL{f},FL(1:(f-1)),'exact'));
end

FullList = FullList(find(FirstOccurence)); % remove multiple instances

JustRead = 0;

if nargin < 3,                            % no data already loaded
  File = zGetNTData(FullList{1},0,SizeCode);
  JustRead = 1;
  Index(1) = 1;
end

for j = 1:length(File),
  LoadedFiles{j} = lower(File(j).Filename);
end

F = length(File);

for f = 1:length(FullList),                       % loop through PDB list
  i = strmatch(lower(FullList{f}), LoadedFiles, 'exact');
  if isempty(i),                                  % if PDB not loaded,
    File(F+1) = zGetNTData(FullList{f},0,SizeCode); %   load it
    F = length(File);
    LoadedFiles = [LoadedFiles FullList{f}];
    Index(f) = F;                                 %   point to it
  else                                            % but if PDB has been loaded
    Index(f) = i(1);                              %   point to first instance
    if length(File(i(1)).NT) == 0,
      File(Index(f)) = zGetNTData(File(Index(f)).Filename,0,SizeCode);
    end
  end
end

else

fprintf('No files specified to read in %s\n', Filenames{1});

end
