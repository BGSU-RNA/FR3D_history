% zGetNTData(Filenames,ReadCode,SizeCode,Verbose) reads data files, depending on ReadCode
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

% If SizeCode = 1, it returns a full version of the file(s)
% If SizeCode = 2, it returns a small version, suitable for motif searching

function [Files] = zGetNTData(Filenames,ReadCode,SizeCode,Verbose)

CurrentVersion = 5.0;                       % version number of class limits

if nargin < 2,
  ReadCode = 0;
end

if nargin < 3,
  SizeCode = 1;
end

if nargin < 4,
  Verbose = 0;
end

path(path,pwd);

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};
end

for f=1:length(Filenames),
  Filename = Filenames{f};
  FILENAME = upper(Filename);
  filename = lower(Filename);

  ClassifyCode = 0;

  if (SizeCode == 1) || (ReadCode > 0),
    ReadFull = 1;
  else
    ReadFull = 0;
  end

  if ReadCode == 4,                     % re-read the PDB file
    File = zReadandAnalyze(Filename,Verbose);   % might not work on a Mac
    ClassifyCode = 1;
  elseif SizeCode >= 2,                 % try to load a small version
    if (exist(strcat(Filename,'_small.mat'),'file') > 0),
      load(strcat(Filename,'_small.mat'),'File','-mat');
      if Verbose > 0,
        fprintf('Loaded %s\n', [Filename '_small.mat']);
      end
    elseif (exist(strcat(FILENAME,'_SMALL.MAT'),'file') > 0),  % helps on a Mac
      load(strcat(FILENAME,'_SMALL.MAT'),'File','-mat');
      if Verbose > 0,
        fprintf('Loaded  %s\n', [FILENAME '_SMALL.MAT']);
      end
    elseif (exist(strcat(filename,'_small.mat'),'file') > 0),  % helps on a Mac
      load(strcat(filename,'_small.mat'),'File','-mat');
      if Verbose > 0,
        fprintf('Loaded   %s\n', [filename '_small.MAT']);
      end
    else
      ReadFull = 1;
    end

    if (ReadFull == 0),
      if File.ClassVersion < CurrentVersion,
         ReadFull = 1;                    % read the full version to classify
      end
    end
  end

  if (ReadFull == 1),                     % try to load the full version
    if (exist(strcat(Filename,'.mat'),'file') > 0),
      load(strcat(Filename,'.mat'),'File','-mat');
      if Verbose > 0,
        fprintf('Loaded %s\n', [Filename '.mat']);
      end
    elseif (exist(strcat(FILENAME,'.MAT'),'file') > 0),  % helps on a Mac
      load(strcat(FILENAME,'.MAT'),'File','-mat');
      if Verbose > 0,
        fprintf('Loaded  %s\n', [FILENAME '.MAT']);
      end
    elseif (exist(strcat(filename,'.mat'),'file') > 0),  % helps on a Mac
      load(strcat(filename,'.mat'),'File','-mat');
      if Verbose > 0,
        fprintf('Loaded   %s\n', [filename '.mat']);
      end
    else
      File = zReadandAnalyze(Filename,Verbose);
      ClassifyCode = 1;
    end
  end

  if ReadFull == 1,
    File.SizeCode = 1;
    if (ReadCode == 2) | (ReadCode == 3) | (ReadCode == 4),
      File = zReadHandFile(File);
    end
    File = zGetPDBInfo(File);          % get resolution and other info
  else
    File.SizeCode = 2;
  end

  if ~isfield(File,'ClassVersion'),
    File.ClassVersion = 0;
  end

  if ~isfield(File,'BasePhosphate'),
    File.BasePhosphate = 0;
  end

  if isfield(File,'Inter'),
    File = rmfield(File,'Inter');
  end

  if isfield(File,'Header'),
    File = rmfield(File,'Header');
  end

  Overlap = 0;

  File.Distance = [];                        % only compute this if needed

  if length(File.NT) > 1,                    % if it has nucleotides,
    if length(File.NT(1).Sugar(:,1)) < 13,
      File = zStoreO3(File);
    end

    c = cat(1,File.NT(1:File.NumNT).Center); % nucleotide centers

    if (ReadCode == 1) | (ReadCode == 3) | (ReadCode == 4) | ... 
      (ClassifyCode == 1) | (File.ClassVersion < CurrentVersion),
      File.Edge = sparse(File.NumNT,File.NumNT);

      File.Distance = zMutualDistance(c,11); % compute distances < 20 Angstroms
      d = sort(nonzeros(File.Distance));

      if ~isempty(d),

       if d(min(10,length(d))) < 1,
         fprintf('%s has overlapping nucleotides and should be avoided as such\n',File.Filename);
         Overlap = 1;
         File.Pair = [];
       else
         t = cputime;
         File = zClassifyPairs(File,Verbose);

         if Verbose > 0,
           fprintf(' Base-phosphate interactions ...');
         end

         File = zPhosphateInteractions(File);
         File.ClassVersion = CurrentVersion;
         ClassifyCode = 1;

         if Verbose > 0,
           fprintf(' Interaction range ... ');
         end

         File = zInteractionRange(File,Verbose);

         if Verbose > 0,
           fprintf('\nClassification took %4.2f minutes\n', (cputime-t)/60);
         end
       end
      end
    end

    if ~isfield(File.NT(1),'Syn'),
      SynList = mSynList(File);
      for k=1:length(File.NT),
        File.NT(k).Syn = SynList(k);
      end
      ClassifyCode = 1;
    end
  else
    File.Pair     = [];
  end

  if ~isfield(File,'Range'),
    File = zInteractionRange(File,Verbose);
    ClassifyCode = 1;
  end

  File = orderfields(File);

  Saved = 0;

  if (ReadCode > 0) || (ClassifyCode > 0),
    zSaveNTData(File,Verbose);
    Saved = 1;
  end

  if SizeCode == 2,
    File = zSmallVersion(File);

    if (ReadFull == 1) && (Saved == 0),
      zSaveNTData(File,Verbose);
    end
  end

  Files(f) = File;

end
