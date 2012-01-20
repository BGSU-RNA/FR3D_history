% pMakeModelsFromLibrary reads the motif library and makes SCFG/MRF models for each motif

loopType = 'HL';                          % hairpins only
loopType = 'IL';                          % internal loops only

Param = [1 4];                            % Verbose, pair substitution method

% ----------------------------------------- Set variables

Verbose = Param(1);
Focus = loopType(1);
Types = {'HL','IL','JL'};                 % types of models we have

% ----------------------------------------- Read file names from the library

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

% ----------------------------------------- Write names of models for JAR3D

fid = fopen(['models' filesep Types{typ} '_Models.txt'],'w');
for m = 1:length(Filenames),
  if strcmp(Types{typ},Filenames(m).name(10:11)),
    fprintf(fid,'%s\n',[Filenames(m).name '.txt']);
  end
end
fclose(fid);

% ----------------------------------------- 

% ----------------------------------------- Load each search, make a model

for m = 1:length(Filenames),
  MN = Filenames(m).name;
  FN = ['SearchSaveFiles' filesep MN '.mat'];
  load(FN,'Search','-mat')                             % Load search data

  if Verbose > 0,
    fprintf('Making a JAR3D model for %s\n', MN);
  end

  % --------------------------------------- Write sequences in FASTA format
  Text = xFASTACandidates(Search.File,Search,0,MN(1:8));

  fid = fopen(['sequences' filesep MN '.fasta'],'w');
  for t = 1:length(Text),
    fprintf(fid,'%s\n',Text{t});
  end
  fclose(fid);

  % --------------------------------------- Make model and write it

  Node = pMakeModelFromSearchSaveFile(Search,Param);

  if ~strcmp(MN(4:5),'99'),                 % don't do this for helices!
    for n = 1:length(Node),
      if fix(abs(Node(n).Edge)) == 1,       % cWW basepair
        Node(n).SubsProb = ones(4,4)/16;    % noninformative
      end
    end
  end

  if strcmp(MN(1:8),'LIB00012'),
    Node(3).SubsProb = [1 1 1 1; 1 1 10 1; 1 1 1 1; 1 1 1 1]/25;
  end

  if strcmp(MN(10:11),'HL'),
    pWriteJavaNodeFile(Search.Query,Node,4,[MN '.txt']);
  elseif strcmp(MN(10:11),'IL'),
    pWriteJavaNodeFile(Search.Query,Node,5,[MN '.txt']);
  elseif strcmp(MN(10:11),'JL'),
    pWriteJavaNodeFile(Search.Query,Node,5,[MN '.txt']);
  end
end

% ----------------------------------------- Run JAR3D on these models
break

%clear java; JAR3D_path;  loopType = 'IL'; clc; S = JAR3DMatlab.MotifTest(pwd,loopType);

%clc; JAR3D_path;  loopType = 'IL'; S = JAR3DMatlab.MotifTest(pwd,loopType);

JAR3D_path;  loopType = 'IL'; S = JAR3DMatlab.MotifTest(pwd,loopType);

[s,t] = size(S);
switch loopType,
case 'IL'
  S = S(:,1:(2*s));
case 'HL'
  S = S(:,1:s);
end

% ----------------------------------------- Display results graphically


loopType = 'IL';

pDisplayModelScores

