% zAddNTData(Filenames,ReadCode,File) reads RNA structure data files, if
% necessary, so that all molecules listed in Filenames are present in File

function [File,Index] = zAddNTData(Filenames,ReadCode,File)

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};
end

if nargin < 3,                            % no data already loaded
  File = zGetNTData(Filenames,ReadCode);
  Index = 1:length(Filenames);
else
  for j = 1:length(File),
    LoadedFiles{j} = lower(File(j).Filename);
  end

  F = length(File);

  for f = 1:length(Filenames),
    i = strmatch(lower(Filenames{f}), LoadedFiles);
    if isempty(i),
      File(F+1) = zGetNTData(Filenames{f},ReadCode);
      F = length(File);
      Index(f) = F;
    else
      Index(f) = i(1);
    end
  end
end
