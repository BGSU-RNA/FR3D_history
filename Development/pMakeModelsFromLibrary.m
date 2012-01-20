% pMakeModelsFromLibrary reads the motif library and makes SCFG/MRF models for each motif

% ----------------------------------------- Read file names from the library

Focus = 'H';                              % hairpins only
Focus = 'I';                              % internal loops only

Types = {'HL','IL','JL'};                 % types of models we have

Filenames = dir(['SearchSaveFiles' filesep 'LIB*']);
keep = [];                               % of all models, which to keep

switch Focus,

case 'H',
  typ = 1;
  for m = 1:length(Filenames),
    if strcmp(Filenames(m).name(1:4),'LIB7') && (Filenames(m).name(9) == '_') && (Filenames(m).name(10) == 'H'),
      keep(m) = 1;
      Filenames(m).name = strrep(Filenames(m).name,'.mat','');
    end 
  end

case 'I',
  typ = 2;
  for m = 1:length(Filenames),
    if strcmp(Filenames(m).name(1:3),'LIB') && (Filenames(m).name(9) == '_') && (Filenames(m).name(10) == 'I'), 
      keep(m) = 1;
      Filenames(m).name = strrep(Filenames(m).name,'.mat','');
    end 
  end

case 'J',
  typ = 3;
  for m = 1:length(Filenames),
    if strcmp(Filenames(m).name(1:3),'LIB') && (Filenames(m).name(9) == '_') && (Filenames(m).name(10) == 'J'), 
      keep(m) = 1;
      Filenames(m).name = strrep(Filenames(m).name,'.mat','');
    end 
  end

end

Filenames = Filenames(find(keep));

% ----------------------------------------- Write names of models

fid = fopen(['models' filesep Types{typ} '_Models.txt'],'w');
for m = 1:length(Filenames),
  if strcmp(Types{typ},Filenames(m).name(10:11)),
    fprintf(fid,'%s\n',[Filenames(m).name '.txt']);
  end
end
fclose(fid);

% ----------------------------------------- Load each seach, make a model

for m = 1:length(Filenames),
  MN = Filenames(m).name;
  FN = ['SearchSaveFiles' filesep MN '.mat'];
  load(FN,'Search','-mat')                             % Load search data

  fprintf('Making a JAR3D model for %s\n', MN);

  % --------------------------------------- Write sequences in FASTA format
  Text = xFASTACandidates(Search.File,Search,1,MN(1:8));

  fid = fopen(['sequences' filesep MN '.fasta'],'w');
  for t = 1:length(Text),
    fprintf(fid,'%s\n',Text{t});
  end
  fclose(fid);

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