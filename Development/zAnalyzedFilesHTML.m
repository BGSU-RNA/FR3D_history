% zAnalyzedFilesHTML produces an .html file listing basepairing and base
% stacking interactions for each molecule in File

function [void] = zAnalyzedFilesHTML(File)

zBackboneCodes;                           % load conformation codes

for f = 1:length(File),
  FN = upper(File(f).Filename);

  LText{1} = ['<a href = "index.html">Return to FR3D home page for ' FN '</a><br>'];
  LText{2} = ['<a href = "' FN '_interactions.html">List of all pairwise interactions in ' FN '</a><br>'];
  LText{3} = ['<a href = "' FN '_basepairs.html">List of basepair interactions in ' FN '</a><br>'];
  LText{4} = ['<a href = "' FN '_stacking.html">List of stacking interactions in ' FN '</a><br>'];
  LText{5} = ['<a href = "' FN '_base_phosphate.html">List of base-phosphate interactions in ' FN '</a><br>'];
  LText{6} = ['<a href = "' FN '_backbone_connectivity.html">List of backbone connectivity relations in ' FN '</a><br>'];
  LText{7} = ['<a href = "' FN '_backbone_conformation.html">List of backbone conformations found in ' FN '</a><br>'];
  LText{8} = ['<a href = "' FN '_motifs.html">List of motifs found in ' FN '</a><br>'];
  LText{9} = ['<a href="http://www.rcsb.org/pdb/explore/explore.do?structureId=' FN '">PDB entry for ' FN '</a><br>'];
  LText{10} = ['<a href="../">Return to list of analyzed structures</a><br>'];
  LText{11} = ['<a href="../../basepairs">FR3D Basepair catalog</a><br>'];
  LText{12} = ['<a href="../../BasePhosphates">FR3D Base-phosphate interactions</a><br>'];
  LText{13} = ['<a href="../../MotifLibrary/index.html">FR3D motif library</a><br>'];
  LText{14} = ['<a href="../../index.html">FR3D home page</a><br>'];
  LText{15} = ['<a href="http://rna.bgsu.edu">BGSU RNA group home page</a><br><br>'];

  HText{1} = ['FR3D classification version ' num2str(File(f).ClassVersion) ' ' date];

  HText{2} = '<p>Basepairing follows the paper Leontis, Stombaugh, Westhof Nucleic Acids Research <b>30</b> No. 16, August 15, 2002.  See the <a href="../../basepairs">FR3D Basepair catalog</a>.  Basepairs are either <i>cis</i> or <i>trans</i>, denoted c and t below.  Each base can use one of three edges, the Waston-Crick edge (W), the Hoogsteen edge (H) or the Sugar edge (S).  Basepairs listed below indicate which base is using which edge.  For example, a line reading A108 G130 tSH would mean that A108 is using its Sugar Edge and G130 is using its Hoogsteen edge.  In the case that both bases use the same edge, a capital letter indicates which base uses the edge in the dominant way.  For perfectly symmetric basepairs such as AU cWW, the capital and lowercase letters are irrelevant.  Bifurcated basepairs are indicated by the text bif.  It does not indicate which base uses which edge.';

  HText{3} = '<p>Base stacking is divided into three subcategories, according to the faces used by each base.  The faces are named this way:  Imagine a base making a Watson-Crick basepair in a helical context.  The side of the base that faces toward the 3'' end of the strand is called the 3 face, while the side that faces the 5'' end is called the 5 face.  The stacking one finds in a helix is thus refered to as s35, meaning that the first base uses the 3 face and the second uses the 5 face.  Two other types of stacking (s33 and s55) occur in other contexts.';

  HText{4} = '<p>Each pair is listed twice.  For instance, if A108 G130 tSH is listed, then so is G130 A108 tHS.  The order in which the edges are listed still corresponds to the base which is using that edge.  Similarly with stacking.  The chain is indicated in parentheses.';

  HText{5} = '<p>Starting from hairpin loops, cWW interactions are classified as being nested if they do not cross any previously established nested cWW pairs.  The last column indicates, for all interactions, the number of nested cWW pairs crossed by the interaction (which may be visualized on the circular diagram of pairwise interactions).';

  HText{6} = '<p>Classification of basepairs and base stacking is done using the program FR3D developed by the <A href="http://rna.bgsu.edu">Bowling Green State University RNA group</a>, including Neocles Leontis, Craig L. Zirbel, Jesse Stombaugh, Ali Mokdad, and Michael Sarver.';

  HText{7} = '<p>Base-phosphate interactions are described in a forthcoming paper by the BGSU RNA group.  See the catalog of <a href="../../BasePhosphates">FR3D Base-phosphate interactions</a>.  Backbone connectivity relations are denoted c35 or c53.  c35 indicates that the first nucleotide uses its O3'' atom to connect to the phosphorus of the second nucleotide and then to the O5'' atom of the second nucleotide.  c53 is used when the order of the nucleotides is reversed.';

  HText{8} = '<p>Backbone conformations are two-character codes which apply to the portion of the backbone shared by the two listed nucleotides.  They are only listed once, with the nucleotides in increasing order.  Nucleotides missing a base in the PDB file are omitted by FR3D.  Conformations are calculated by Dangle and Suitename from the <A href="http://kinemage.biochem.duke.edu/index.php">Richardson lab</A> at Duke University.'; 

  HText{9} = '<p>Please write to Craig L. Zirbel at <img src="http://www-math.bgsu.edu/z/cz.gif" align="absbottom"> with questions or comments.<hr>';

  HText{10} = '<pre>';

% --------------------------------------------- Produce interaction list

  c = 1;                                 % counter for pairs

  IText{1} = '';
  InterType = [];

  E = File(f).Edge;
  E = E .* (abs(E)>0) .* (abs(E) < 30);     % basepairing and stacking

  BPh = File(f).BasePhosphate;
  BPh = BPh .* (BPh < 20) .* (BPh > 0);     % base phosphates

  BC = File(f).Covalent;                    % covalent connections

  for i = 1:File(f).NumNT,
    N1 = File(f).NT(i);

    % ------------------------------------- Find basepairing, stacking
    j = find(E(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = sprintf('%4d', File(f).Crossing(i,j(k)));
      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %7s - %s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zEdgeText(File(f).Edge(i,j(k)),0,N1.Code,N2.Code),r);

      InterType(c) = abs(File(f).Edge(i,j(k)));

      c = c + 1;
    end

    % ------------------------------------- Find base phosphate interactions
    j = find(BPh(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = sprintf('%4d', File(f).Range(i,j(k)));

      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %7s - %s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zBasePhosphateText(BPh(i,j(k)),1), r);

      if i == j(k),
        InterType(c) = 200.1;                 % self interaction
      else
        InterType(c) = 200;                   % non-self interaction
      end

      c = c + 1;
    end

    % ---------------------------------- Find backbone connectivity relations
    j = find(BC(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = sprintf('%4d', 0);

      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %6s  - %s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zBackboneContinuityText(BC(i,j(k))),r);

      InterType(c) = 300;                 % code for later

      c = c + 1;
    end

    % ---------------------------------- Find backbone conformation
    j = find(File(f).Backbone(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = sprintf('%4d', 0);

      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %6s  - %s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, Codes{File(f).Backbone(i,j(k))},r);

      InterType(c) = 400;                 % code for later

      c = c + 1;
    end
  end

% -------------------------------------------------- Write html files

  warning off

  mypath = [pwd filesep 'Web' filesep 'AnalyzedStructures'];
  mkdir(mypath, FN);
  warning on

  mypath = [mypath filesep FN filesep];

% -------------------------------------------------- Write index.html file

  fid = fopen([mypath 'index.html'],'w');

  fprintf(fid,'<html>\n<head>\n<title>%s analysis by FR3D</title>\n<script src="../../jmol/Jmol.js" type="text/javascript"></script></head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D analysis of ' FN '.pdb</h1>']);

  for L = 2:length(LText),
    fprintf(fid,'%s\n', LText{L});
  end

  Diagram = ['<a href="' FN '_circular_diagram.pdf"> <img src="' FN '_circular_diagram.png" alt="Click for high resolution pdf" width = "615" > </a>'];

  DiagramText = '<p>Circular basepair diagram in Nussinov style.  Click the diagram for a high-resolution PDF version.';

%  Dark blue chords indicate nested Watson-Crick basepairs, cyan indicates nested non-Watson-Crick basepairs, red indicates non-nested Watson-Crick basepairs, green indicates non-nested non-Watson-Crick basepairs, and yellow indicates long-range stacking interactions.

  fprintf(fid,'<table border="1">\n<tr>\n<th>Circular interaction diagram</th>\n<th>RNA 3D structure</th>\n</tr>\n');

  fprintf(fid,'<tr>\n<td>\n%s\n</td>\n',Diagram);

  fprintf(fid,'<td>\n<script type=''text/javascript''>\n jmolInitialize(''../../jmol'');\n jmolApplet(600, ''load %s_RNA.pdb;spacefill off'',''3D structure'');\n </script>\n</td>\n</tr>\n', FN);

  fprintf(fid,'<tr>\n<td>\n%s\n</td>\n', DiagramText);

  fprintf(fid,'<td><a href="%s_RNA.pdb">Click here</a> to download the RNA coordinate file.\n</td>\n</tr>\n</table>\n',FN);

  fprintf(fid,'<b>Resolution: </b>%7.1f<br>\n', File(f).Info.Resolution);
  fprintf(fid,'<b>Descriptor: </b>%s<br>\n', File(f).Info.Descriptor);
  fprintf(fid,'<b>Experimental technique: </b>%s<br>\n', File(f).Info.ExpTechnique);
  fprintf(fid,'<b>Release Date: </b>%s<br>\n', File(f).Info.ReleaseDate);
  fprintf(fid,'<b>Author: </b>%s<br>\n', File(f).Info.Author);
  if ~isempty(File(f).Info.Source),
    fprintf(fid,'<b>Biological source: </b>%s<br>\n', File(f).Info.Source);
  end

  fprintf(fid,'</html>\n');

  fclose(fid);

% ----------------------------------------------- Write FN_interactions file  

  fid = fopen([mypath FN '_interactions.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s interactions from FR3D</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D list of pairwise interactions in ' FN '</h1>']);

  for L = 1:length(LText),
    if L ~= 2,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  for i = 1:length(IText),
    fprintf(fid,'%5d %s\n',i,IText{i});
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

% ----------------------------------------------- Write FN_basepairs file  

  fid = fopen([mypath FN '_basepairs.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s basepairs from FR3D</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D list of basepair interactions in ' FN '</h1>']);

  for L = 1:length(LText),
    if L ~= 3,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  k = find(InterType < 15);

  for i = 1:length(IText(k)),
    fprintf(fid,'%5d %s\n',i,IText{k(i)});
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

% ----------------------------------------------- Write FN_stacking file  

  fid = fopen([mypath FN '_stacking.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s stacking interactions from FR3D</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D list of stacking interactions in ' FN '</h1>']);

  for L = 1:length(LText),
    if L ~= 4,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  k = find((InterType > 19) .* (InterType < 100));

  for i = 1:length(IText(k)),
    fprintf(fid,'%5d %s\n',i,IText{k(i)});
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

% --------------------------------------------- Write FN_base_phosphate file  

  fid = fopen([mypath FN '_base_phosphate.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s base-phosphate interactions from FR3D</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D list of base-phosphate interactions in ' FN '</h1>']);

  for L = 1:length(LText),
    if L ~= 5,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:(length(HText)-2),
    fprintf(fid,'%s\n',HText{i});
  end

  fprintf(fid,'<p>Non-self interactions listed first, if any, then self interactions.\n');

  for i = (length(HText)-1):length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  k = find(InterType == 200);               % non-self interactions

  c = 1;

  for i = 1:length(IText(k)),
    fprintf(fid,'%5d %s\n',c,IText{k(i)});
    c = c + 1;
  end

  k = find(InterType == 200.1);               % self interactions

  for i = 1:length(IText(k)),
    fprintf(fid,'%5d %s\n',c,IText{k(i)});
    c = c + 1;
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

% --------------------------------------- Write FN_backbone_connectivity file  

  fid = fopen([mypath FN '_backbone_connectivity.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s backbone connectivity relations from FR3D</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>FR3D backbone connectivity relations in ' FN '</h1>']);

  for L = 1:length(LText),
    if L ~= 6,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  k = find(InterType == 300);

  for i = 1:length(IText(k)),
    fprintf(fid,'%5d %s\n',i,IText{k(i)});
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

% --------------------------------------- Write FN_backbone_conformation file  

  fid = fopen([mypath FN '_backbone_conformation.html'],'w'); % open for writing

  fprintf(fid,'<html>\n<head>\n<title>%s backbone conformation relations from Richardson lab programs</title>\n</head>\n<body>',FN);
  fprintf(fid,'%s\n',['<h1>Backbone conformation relations in ' FN '</h1>']);
  fprintf(fid,'As computed by Dangle and Suitename from the <A href="http://kinemage.biochem.duke.edu/index.php">Richardson lab</A> at Duke University</p>\n');

  for L = 1:length(LText),
    if L ~= 7,
      fprintf(fid,'%s\n', LText{L});
    end
  end

  for i = 1:length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  k = find(InterType == 400);

  for i = 1:length(IText(k)),
    fprintf(fid,'%5d %s\n',i,IText{k(i)});
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);


  clear IText HText 

end


