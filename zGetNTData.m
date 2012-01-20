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

  if (ReadCode < 4) & (exist(strcat(Filename,'.mat'),'file') > 0),
      load(strcat(Filename,'.mat'),'File','-mat');
      File.Distance = tril(File.Distance) + tril(File.Distance)';  
                                          % only saved lower triangular part
      fprintf('Loaded %s\n', [Filename '.mat']);
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

  if (ReadCode == 1) | (ReadCode == 3) | (ReadCode == 4) | ... 
    (ClassifyCode == 1) | (length(fieldnames(File)) < 11) | ...
    (max(max(File.Inter)) < 100),
    File.Edge = sparse(File.NumNT,File.NumNT);
    File = zClassifyPairs(File);
    File = zUpdateDistanceToExemplars(File);
    ClassifyCode = 1;
  end

  if ~isfield(File.NT(1),'Syn'),
    SynList = mSynList(File);
    for k=1:length(File.NT),
      File.NT(k).Syn = SynList(k);
    end
    ClassifyCode = 1;
 end

  File = orderfields(File);

  if ((ReadCode > 0) | (ClassifyCode > 0)) & (File.NumNT > 0),
    File.Modified = 0;
    zSaveNTData(File);
  end

  Files(f) = File;
end
