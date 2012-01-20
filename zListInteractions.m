% zListInteractions produces an .html file listing basepairing and base
% stacking interactions for each molecule in File

function [void] = zListInteractions(File)


for f = 1:length(File),
  clear Text
  Text{1} = ['<html><title>' File(f).Filename ' basepairing and base stacking</title>'];
  Text{2} = ['<h1>' File(f).Filename ' basepairing and base stacking</h1>'];
  Text{3} = ['FR3D classification version ' num2str(File(f).ClassVersion) ' ' date];
  Text{4} = '<p>Basepairing follows the paper Leontis, Stombaugh, Westhof Nucleic Acids Research <b>30</b> No. 16, August 15, 2002.  Basepairs are either <i>cis</i> or <i>trans</i>, denoted c and t below.  Each base can use one of three edges, the Waston-Crick edge (W), the Hoogsteen edge (H) or the Sugar edge (S).  Basepairs listed below indicate which base is using which edge.  For example, a line reading A108 G130 tSH would mean that A108 is using its Sugar Edge and G130 is using its Hoogsteen edge.  In the case that both bases use the same edge, a capital letter indicates which base uses the edge in the dominant way.  For perfectly symmetric basepairs such as AU cWW, the capital and lowercase letters are irrelevant.  Bifurcated basepairs are indicated by the text bif.  It does not indicate which base uses which edge.';

  Text{5} = '<p>Base stacking is divided into three subcategories, according to the faces used by each base.  The faces are named this way:  Imagine a base making a Watson-Crick basepair in a helical context.  The side of the base that faces toward the 3'' end of the strand is called the 3 face, while the side that faces the 5'' end is called the 5 face.  The stacking one finds in a helix is thus refered to as s35, meaning that the first base uses the 3 face and the second uses the 5 face.  Two other types of stacking (s33 and s55) occur in other contexts.';

  Text{6} = '<p>Each pair is listed twice.  For instance, if A108 G130 tSH is listed, then so is G130 A108 tHS.  The order in which the edges are listed still corresponds to the base which is using that edge.  Similarly with stacking.  The chain is indicated in parentheses.';

  Text{7} = '<p>Classification of basepairs and base stacking is done using the program FR3D developed by the <A href="http://rna.bgsu.edu">Bowling Green State University RNA group</a>, including Neocles Leontis, Craig L. Zirbel, Jesse Stombaugh, Ali Mokdad, and Michael Sarver.';

  Text{8} = '<p>Please write to Craig L. Zirbel at <img src="http://www-math.bgsu.edu/z/cz.gif" align="absbottom"> with questions or comments.<hr>';

  Text{9} = '<pre>';

  t = length(Text);
  c = 1;                                 % counter for pairs

  E = File(f).Edge;
  E = E .* (abs(E)>0) .* (abs(E) < 30); 
  for i = 1:File(f).NumNT,
    N1 = File(f).NT(i);
    j = find(E(i,:));
    for k = 1:length(j),
      N2 = File(f).NT(j(k));
      Text{c+t} = sprintf('%5d  %s%4s(%s) - %s%4s(%s) - %4s', c, N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zEdgeText(File(f).Edge(i,j(k))));
      c = c + 1;
    end
  end

  t = length(Text);
  Text{t+1} = '</pre></html>';

  fid = fopen(['Interactions' filesep File(f).Filename '_Interaction_List.html'],'w'); % open for writing
  for i = 1:length(Text),
    fprintf(fid,'%s\n',Text{i});
  end

  fclose(fid);

end

