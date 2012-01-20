% pMakeModelsFromLibrary reads the motif library and makes SCFG/MRF models for each motif

% ----------------------------------------- Read file names from the library

Filenames = dir(['SearchSaveFiles' filesep 'LIB*']);

keep = [];

for m = 1:length(Filenames),
  if strcmp(Filenames(m).name(1:3),'LIB') && (Filenames(m).name(9) == '_') && (Filenames(m).name(10) == 'I'),             % internal loops only for now
    keep(m) = 1;
    Filenames(m).name = strrep(Filenames(m).name,'.mat','');
  end
end

Filenames = Filenames(find(keep));

% ----------------------------------------- Write names of models

Types = {'HL','IL','JL'};

for t = 1:length(Types),
  fid = fopen(['models' filesep Types{t} '_Models.txt'],'w');
  for m = 1:length(Filenames),
    if strcmp(Types{t},Filenames(m).name(10:11)),
      fprintf(fid,'%s\n',[Filenames(m).name '.txt']);
    end
  end
  fclose(fid);
end

% ----------------------------------------- Load each seach, make a model

for m = 1:length(Filenames),
  MN = Filenames(m).name;
  FN = ['SearchSaveFiles' filesep MN '.mat'];
  load(FN,'Search','-mat')                             % Load search data

  fprintf('Making a JAR3D model for %s\n', MN);

  % --------------------------------------- Write sequences in FASTA format
  Text = xFASTACandidates(Search.File,Search,1,MN(1:8));

%fprintf('%s\n',MN);

  fid = fopen(['sequences' filesep MN '.fasta'],'w');
  for t = 1:length(Text),
    fprintf(fid,'%s\n',Text{t});

%fprintf('%s\n',Text{t});

  end
  fclose(fid);

%fprintf('\n');

  % --------------------------------------- Make model and write it

  Node = pMakeModelFromSearchSaveFile(Search,MN(10:11),1);

  if strcmp(MN(10:11),'HL'),
    pWriteJavaNodeFile(Search.Query,Node,4,[MN '.txt']);
  elseif strcmp(MN(10:11),'IL'),
    pWriteJavaNodeFile(Search.Query,Node,5,[MN '.txt']);
  elseif strcmp(MN(10:11),'JL'),
    pWriteJavaNodeFile(Search.Query,Node,5,[MN '.txt']);
  end

  if m < length(Filenames)
%    pause
  end
end