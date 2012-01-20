% zGetNTData reads data files, depending on ReadCode
%
% If ReadCode = 0, it looks for Filename.mat and reads it if it exists.
%    If no hand classification data is present, it looks for Filename.class.
% If ReadCode = 1, it looks for Filename.mat, reads it, and re-does 
%    the classification of interacting pairs
% If ReadCode = 2, it looks for Filename.mat, reads it, and 
%    re-reads the hand classification file Filename.hand.
% If ReadCode = 3, it looks for Filename.mat, reads it, re-does the
%    classification of pairs and reads the hand classification file
% If ReadCode = 4, it reads Filename.pdb, analyzes each nucleotide, reads the
%    hand classification file, and classifies interacting pairs

function [Files] = zGetNTData(Filenames,ReadCode)

CurrentVersion = 3.1;                       % version number of class limits

if nargin < 2,
  ReadCode = 0;
end

path(path,pwd);

if ~(exist('PDBFiles') == 7),        % if directory doesn't yet exist
  mkdir('PDBFiles');
end
path(path,[pwd filesep 'PDBFiles']);

if ~(exist('PrecomputedData') == 7),        % if directory doesn't yet exist
  mkdir('PrecomputedData');
end
path(path,[pwd filesep 'PrecomputedData']);

if ~(exist('SearchSaveFiles') == 7),        % if directory doesn't yet exist
  mkdir('SearchSaveFiles');
end
path(path,[pwd filesep 'SearchSaveFiles']);

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};
end

for f=1:length(Filenames),
  Filename = Filenames{f};
  FILENAME = upper(Filename);
  filename = lower(Filename);
  if (ReadCode < 4) & (exist(strcat(Filename,'.mat'),'file') > 0),
      load(strcat(Filename,'.mat'),'File','-mat');
      fprintf('Loaded %s\n', [Filename '.mat']);
      ClassifyCode = 0;
  elseif (ReadCode < 4) & (exist(strcat(FILENAME,'.MAT'),'file') > 0),  % helps on a Mac
      load(strcat(FILENAME,'.MAT'),'File','-mat');
      fprintf('Loaded %s\n', [FILENAME '.MAT']);
      ClassifyCode = 0;
  elseif (ReadCode < 4) & (exist(strcat(filename,'.mat'),'file') > 0),  % helps on a Mac
      load(strcat(filename,'.mat'),'File','-mat');
      fprintf('Loaded %s\n', [FILENAME '.MAT']);
      ClassifyCode = 0;
  else
      File = zReadandAnalyze(Filename);
      ClassifyCode = 1;
  end

  if isfield(File,'BI'),
    File = rmfield(File,'BI');
  end

  if isfield(File,'BermanClass'),
    File = rmfield(File,'BermanClass');
  end

  if (ReadCode == 2) | (ReadCode == 3) | (ReadCode == 4),
    File = zReadHandFile(File);
  end

  Overlap = 0;

  if ~isfield(File,'ClassVersion'),
    File.ClassVersion = 0;
  end

  if length(File.NT) > 0,                    % if it has nucleotides,

    c = cat(1,File.NT(1:File.NumNT).Center);
    File.Distance = zMutualDistance(c,35); 

    if (ReadCode == 1) | (ReadCode == 3) | (ReadCode == 4) | ... 
      (ClassifyCode == 1) | (File.ClassVersion < CurrentVersion),
      File.Edge = sparse(File.NumNT,File.NumNT);

      d = sort(nonzeros(File.Distance));
      if d(min(10,length(d))) < 1,
        fprintf('%s has overlapping nucleotides and should be avoided as such\n',File.Filename);
        Overlap = 1;
      else
        File = zClassifyPairs(File);
        File = zUpdateDistanceToExemplars(File);
        File.ClassVersion = CurrentVersion;
        ClassifyCode = 1;
      end
    end

    if ~isfield(File.NT(1),'Syn'),
      SynList = mSynList(File);
      for k=1:length(File.NT),
        File.NT(k).Syn = SynList(k);
      end
      ClassifyCode = 1;
    end

%File.Header = zExtractAtomsPDB(Filename,'##TempPDB');

    if ~isfield(File,'Header'),
      File.Header.ModelStart = [];
      File.Header.ExpData    = '';
      File.Header.Resolution = '';
      ClassifyCode = 1;
    end

    File = zGetPDBInfo(File);          % get resolution and other info

    File = orderfields(File);

    if Overlap == 0,
      if ((ReadCode > 0) | (ClassifyCode > 0)) & (File.NumNT > 0),
        zSaveNTData(File);
      end

      Files(f) = File;
    end
  else
    Files(f) = File;
  end
end
