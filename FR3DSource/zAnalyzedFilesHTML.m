% zAnalyzedFilesHTML produces an .html file listing basepairing and base
% stacking interactions for each molecule in File

function [void] = zAnalyzedFilesHTML(File)

for f = 1:length(File),
  FN = upper(File(f).Filename);

  fprintf('Writing HTML files for %s, file %d\n', FN, f);

  LText{1} = ['<a href = "index.html">Return to FR3D home page for ' FN '</a><br>'];
  LText{2} = ['<a href = "' FN '_interactions.html">List of all pairwise interactions in ' FN '</a><br>'];
  LText{3} = ['<a href = "' FN '_basepairs.html">List of basepair interactions in ' FN '</a><br>'];
  LText{4} = ['<a href = "' FN '_stacking.html">List of stacking interactions in ' FN '</a><br>'];
%  LText{4} = ['<a href = "' FN '_basephosphate.html">List of base-phosphate interactions in ' FN '</a><br>'];
  LText{5} = ['<a href = "' FN '_motifs.html">List of motifs found in ' FN '</a><br>'];
  LText{6} = ['<a href="http://www.rcsb.org/pdb/explore/explore.do?structureId=' FN '">PDB entry for ' FN '</a><br>'];
  LText{7} = ['<a href="../">Return to list of analyzed structures</a><br>'];
  LText{8} = ['<a href="../../basepairs">Basepair catalog</a><br>'];
  LText{9} = ['<a href="../../MotifLibrary/index.html">FR3D motif library</a><br>'];
  LText{10} = ['<a href="../../index.html">FR3D home page</a><br>'];
  LText{11} = ['<a href="http://rna.bgsu.edu">BGSU RNA group home page</a><br><br>'];

  HText{1} = ['FR3D classification version ' num2str(File(f).ClassVersion) ' ' date];

  HText{2} = '<p>Basepairing follows the paper Leontis, Stombaugh, Westhof Nucleic Acids Research <b>30</b> No. 16, August 15, 2002.  Basepairs are either <i>cis</i> or <i>trans</i>, denoted c and t below.  Each base can use one of three edges, the Waston-Crick edge (W), the Hoogsteen edge (H) or the Sugar edge (S).  Basepairs listed below indicate which base is using which edge.  For example, a line reading A108 G130 tSH would mean that A108 is using its Sugar Edge and G130 is using its Hoogsteen edge.  In the case that both bases use the same edge, a capital letter indicates which base uses the edge in the dominant way.  For perfectly symmetric basepairs such as AU cWW, the capital and lowercase letters are irrelevant.  Bifurcated basepairs are indicated by the text bif.  It does not indicate which base uses which edge.';

  HText{3} = '<p>Base stacking is divided into three subcategories, according to the faces used by each base.  The faces are named this way:  Imagine a base making a Watson-Crick basepair in a helical context.  The side of the base that faces toward the 3'' end of the strand is called the 3 face, while the side that faces the 5'' end is called the 5 face.  The stacking one finds in a helix is thus refered to as s35, meaning that the first base uses the 3 face and the second uses the 5 face.  Two other types of stacking (s33 and s55) occur in other contexts.';

  HText{4} = '<p>Each pair is listed twice.  For instance, if A108 G130 tSH is listed, then so is G130 A108 tHS.  The order in which the edges are listed still corresponds to the base which is using that edge.  Similarly with stacking.  The chain is indicated in parentheses.';

  HText{5} = '<p>Classification of basepairs and base stacking is done using the program FR3D developed by the <A href="http://rna.bgsu.edu">Bowling Green State University RNA group</a>, including Neocles Leontis, Craig L. Zirbel, Jesse Stombaugh, Ali Mokdad, and Michael Sarver.';

  HText{6} = '<p>Please write to Craig L. Zirbel at <img src="http://www-math.bgsu.edu/z/cz.gif" align="absbottom"> with questions or comments.<hr>';

  HText{7} = '<pre>';

% --------------------------------------------- Produce interaction list

  c = 1;                                 % counter for pairs

  IText{1} = '';
  InterType = [];

  E = File(f).Edge;
  E = E .* (abs(E)>0) .* (abs(E) < 30); 
  for i = 1:File(f).NumNT,
    N1 = File(f).NT(i);
    j = find(E(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = '';
      if (File(f).Range(i,j(k)) == 0) && abs(E(i,j(k))) < 15,
        r = 'Nested';
      elseif (File(f).Range(i,j(k)) <= 10) && abs(E(i,j(k))) < 15,
        r = 'Local';
      elseif File(f).Range(i,j(k)) > 10,
        r = 'Long-range';
      end
      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %5s - %10s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zEdgeText(File(f).Edge(i,j(k)),0),r);

      InterType(c) = abs(File(f).Edge(i,j(k)));

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

  DiagramText = '<p>Circular basepair diagram in Nussinov style.  Click the diagram for a high-resolution PDF version.  Dark blue chords indicate nested Watson-Crick basepairs, cyan indicates nested non-Watson-Crick basepairs, red indicates non-nested Watson-Crick basepairs, green indicates non-nested non-Watson-Crick basepairs, and yellow indicates long-range stacking interactions.';

  fprintf(fid,'<table border="1">\n<tr>\n<th>Circular interaction diagram</th>\n<th>RNA 3D structure</th>\n</tr>\n');

  fprintf(fid,'<tr>\n<td>\n%s\n</td>\n',Diagram);

  fprintf(fid,'<td>\n<script type=''text/javascript''>\n jmolInitialize(''../../jmol'');\n jmolApplet(600, ''load %s_RNA.pdb;spacefill off'',''3D structure'');\n </script>\n</td>\n</tr>\n', FN);

  fprintf(fid,'<tr>\n<td>\n%s\n</td>\n', DiagramText);

  fprintf(fid,'<td><a href="%s_RNA.pdb>Click here</a> to download the RNA coordinate file.\n</td>\n</tr>\n</table>\n',FN);

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

  k = find(InterType >19);

  for i = 1:length(IText(k)),
    fprintf(fid,'%5d %s\n',i,IText{k(i)});
  end

  fprintf(fid,'</pre>\n</html>\n');

  fclose(fid);

% ----------------------------------------------- Create circular diagram

%  figure(1)
  clf
  E = fix(abs(File(f).Edge));
  B = E .* (E > 0) .* (E < 24);                 % pairs and stacks
  S = File(f).Range;
  zNussinovPlot(File(f), triu(B), (B==1).*(S==0) + 2*(B>1).*(B<14).*(S==0) + 3*(B==1).*(S>0) + 4*(B > 1).*(B < 14) .*(S>0) + 5*(B > 20) .* (B < 25) .* (S > 10),1);

  saveas(gcf,[mypath FN '_circular_diagram.png'],'png');
  [X,map] = imread([mypath FN '_circular_diagram.png']);
  Y = X(30:830,210:1030,:);
  imwrite(Y,[mypath FN '_circular_diagram.png']);

  clf
  zNussinovPlot(File(f), triu(B), (B==1).*(S==0) + 2*(B>1).*(B<14).*(S==0) + 3*(B==1).*(S>0) + 4*(B > 1).*(B < 14) .*(S>0) + 5*(B > 20) .* (B < 25) .* (S > 10),0.1);
  saveas(gcf,[mypath FN '_circular_diagram.pdf'],'pdf');

% ----------------------------------------------- Write PDB file

  zWritePDB(File(f),[mypath FN '_RNA.pdb']);

  clear IText HText 

end

