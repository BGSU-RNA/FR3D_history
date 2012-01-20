% zListInteractions produces an .html file listing basepairing and base
% stacking interactions for each molecule in File

function [void] = zListInteractions(File)

for f = 1:length(File),

  fprintf('Writing interaction list for %s, file %d\n', File(f).Filename, f);

  HText{1} = ['<html><title>' File(f).Filename ' basepairing and base stacking</title>'];
  HText{2} = ['<h1>' File(f).Filename ' basepairing and base stacking</h1><a href = "index.html">List of motifs found in ' File(f).Filename '</a><br><a href="../">Return to list of PDB files</a><br><a href="../../MotifLibrary/index.html">Return to motif library</a><br>'];
  HText{3} = ['FR3D classification version ' num2str(File(f).ClassVersion) ' ' date];
  HText{4} = '<p>Basepairing follows the paper Leontis, Stombaugh, Westhof Nucleic Acids Research <b>30</b> No. 16, August 15, 2002.  Basepairs are either <i>cis</i> or <i>trans</i>, denoted c and t below.  Each base can use one of three edges, the Waston-Crick edge (W), the Hoogsteen edge (H) or the Sugar edge (S).  Basepairs listed below indicate which base is using which edge.  For example, a line reading A108 G130 tSH would mean that A108 is using its Sugar Edge and G130 is using its Hoogsteen edge.  In the case that both bases use the same edge, a capital letter indicates which base uses the edge in the dominant way.  For perfectly symmetric basepairs such as AU cWW, the capital and lowercase letters are irrelevant.  Bifurcated basepairs are indicated by the text bif.  It does not indicate which base uses which edge.';

  HText{5} = '<p>Base stacking is divided into three subcategories, according to the faces used by each base.  The faces are named this way:  Imagine a base making a Watson-Crick basepair in a helical context.  The side of the base that faces toward the 3'' end of the strand is called the 3 face, while the side that faces the 5'' end is called the 5 face.  The stacking one finds in a helix is thus refered to as s35, meaning that the first base uses the 3 face and the second uses the 5 face.  Two other types of stacking (s33 and s55) occur in other contexts.';

  HText{6} = '<p>Each pair is listed twice.  For instance, if A108 G130 tSH is listed, then so is G130 A108 tHS.  The order in which the edges are listed still corresponds to the base which is using that edge.  Similarly with stacking.  The chain is indicated in parentheses.';

  HText{7} = '<p>Classification of basepairs and base stacking is done using the program FR3D developed by the <A href="http://rna.bgsu.edu">Bowling Green State University RNA group</a>, including Neocles Leontis, Craig L. Zirbel, Jesse Stombaugh, Ali Mokdad, and Michael Sarver.';

  HText{8} = '<p>Please write to Craig L. Zirbel at <img src="http://www-math.bgsu.edu/z/cz.gif" align="absbottom"> with questions or comments.<hr>';

  HText{9} = ['<a href="' File(f).Filename '_Circular_diagram.pdf"> <img src="' File(f).Filename '_Circular_diagram.png" alt="Click for high resolution pdf" width = "50%" > </a>'];

  HText{10} = '<p>Circular basepair diagram in Nussinov style.  Dark blue chords indicate nested Watson-Crick basepairs, cyan indicates nested non-Watson-Crick basepairs, red indicates non-nested Watson-Crick basepairs, green indicates non-nested non-Watson-Crick basepairs, and yellow indicates long-range stacking interactions.';

  HText{11} = '<pre>';

% --------------------------------------------- Produce interaction list

  c = 1;                                 % counter for pairs

  Text{1} = '';

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
      Text{c} = sprintf('%5d  %s%4s(%s) - %s%4s(%s) - %5s - %10s', c, N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zEdgeText(File(f).Edge(i,j(k)),0),r);
      c = c + 1;
    end
  end

% -------------------------------------------------- Write HTML file

  warning off
  mkdir([pwd filesep 'Interactions'], File(f).Filename);
  warning on

  fid = fopen(['Interactions' filesep File(f).Filename filesep File(f).Filename '_interactions.html'],'w'); % open for writing

%['Interactions' filesep File(f).Filename filesep File(f).Filename '_interactions.html']

  for i = 1:length(HText),
    fprintf(fid,'%s\n',HText{i});
  end

  for i = 1:length(Text),
    fprintf(fid,'%s\n',Text{i});
  end

  fprintf(fid,'</pre></html>\n');

  fclose(fid);




  figure(4)
  clf
  E = fix(abs(File(f).Edge));
  B = E .* (E > 0) .* (E < 24);                 % pairs and stacks
  S = File(f).Range;
  zNussinovPlot(triu(B), (B==1).*(S==0) + 2*(B>1).*(B<14).*(S==0) + 3*(B==1).*(S>0) + 4*(B > 1).*(B < 14) .*(S>0) + 5*(B > 20) .* (B < 25) .* (S > 10),1);
  saveas(gcf,['Interactions' filesep File(f).Filename filesep File(f).Filename '_Circular_diagram.png'],'png');

  clf
  zNussinovPlot(triu(B), (B==1).*(S==0) + 2*(B>1).*(B<14).*(S==0) + 3*(B==1).*(S>0) + 4*(B > 1).*(B < 14) .*(S>0) + 5*(B > 20) .* (B < 25) .* (S > 10),0.1);
  saveas(gcf,['Interactions' filesep File(f).Filename filesep File(f).Filename '_Circular_diagram.pdf'],'pdf');


  clear Text, HText

end
