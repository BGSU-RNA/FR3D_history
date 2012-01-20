% zAddNTData(Filenames,ReadCode,File) reads RNA structure data files, if
% necessary, so that all molecules listed in Filenames are present in File.
% The parameter Index is a non-redundant list of indices of File
% corresponding to names in Filenames.

% SizeCode = 0, 1  : load full .mat files
% SizeCode = 2     : load _small.mat files
% SizeCode = 3     : load, but do not append to File (for reclassification)
% SizeCode = 4     : load _small.mat files, do not compute distances

function [File,Index] = zAddNTData(Filenames,SizeCode,File,PDBStart)

if nargin < 2,
  SizeCode = 1;                           % default is to read full files
end

LoadedFiles = [];
F = 0;

if nargin == 3,
  F = length(File);
  for j = 1:length(File),
    LoadedFiles{j} = lower(File(j).Filename);
  end
  if isempty(File),
    clear File
  end
end

if SizeCode == 0,
  SizeCode = 1;
end

if SizeCode == 3,
  File = [];
end

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};                % make into a cell array
end

% ----------------------------------------- Read PDB lists, if any

FullList = [];

for j=1:length(Filenames),
  FullList = [FullList; zReadPDBList(Filenames{j})];
end

% ----------------------------------------- Skip some files

if nargin == 4,
  if strcmp(PDBStart,'back') == 1,
    FullList = FullList(end:-1:1);
  else
    keep = [];
    for j=1:length(FullList),
      if issorted([PDBStart(1:4); FullList{j}(1:4)],'rows'),
        keep = [keep j];
      end
    end
    FullList = FullList(keep);
  end
end

% ----------------------------------------- Read PDB files

if length(FullList) > 0,

for f = 1:length(FullList),                       % loop through PDB list
  if ~isempty(FullList{f}),
    i = strmatch(lower(FullList{f}), LoadedFiles, 'exact');
    if isempty(i),                                  % if PDB not loaded,
      NewF = zGetNTData(FullList{f},0,SizeCode); %   load it
      if SizeCode ~= 3,
        File(F+1) = NewF;
      end
      clear NewF;
      F = length(File);
      k = length(LoadedFiles);
      LoadedFiles{k+1} = FullList{f};
      Index(f) = F;                           %   point to it
    else                                      % but if PDB has been loaded
      Index(f) = i(1);                        %   point to first instance
      if length(File(i(1)).NT) == 0,
        NewF = zGetNTData(File(Index(f)).Filename,0,SizeCode);
        if SizeCode ~= 3,
          File(Index(f)) = NewF;
          clear NewF;
        end
      end
    end
  end
end

% -----------------------------------  allow for File(Index)
F = length(File);

for i = 1:length(Index),             
  if Index(i) == 0,
    File(F+1).Filename = 'Fictitious';
    File(F+1).NumNT = 0;               % create a fictitious file
    Index(i) = F+1;                    % point to the fictitious file
  end
end

else

fprintf('No files specified to read in %s\n', Filenames{1});

end
