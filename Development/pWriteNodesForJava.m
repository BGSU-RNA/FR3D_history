% pWriteNodesForJave(File,Nodes) writes Java code describing each node

function [void] = pWriteNodesForJave(File,Node)

fid = fopen('Parser/NodeFile.java','w');

Text{1} = '';
t = 1;

for n=1:length(Node),
  switch Node(n).type
    case 'Initial'
      if n == 1,
        Text{t} = [Text{t} sprintf('first = new InitialNode(null')];
      else
        Text{t} = [Text{t} sprintf('current = new InitialNode(current')];
      end

      Text{t} = [Text{t} sprintf(', new double[]')];
      Text{t} = [Text{t} subWrite(fid, Node(n).leftLengthDist)];

      Text{t} = [Text{t} sprintf(', new double[]')];
      Text{t} = [Text{t} subWrite(fid, Node(n).leftLetterDist)];

      Text{t} = [Text{t} sprintf(', new double[]')];
      Text{t} = [Text{t} subWrite(fid, Node(n).rightLengthDist)];

      Text{t} = [Text{t} sprintf(', new double[]')];
      Text{t} = [Text{t} subWrite(fid, Node(n).rightLetterDist)];

      Text{t} = [Text{t} sprintf(');\n')];

      if n == 1,
        Text{t} = [Text{t} sprintf('current = first;\n')];
      end

    case 'Basepair'
      Text{t} = [Text{t} sprintf('current = new BasepairNode(current, ')];
      Text{t} = [Text{t} sprintf('%0.6f', Node(n).Delete)];

      tt = Node(n).PIns(1:16);
      tt = tt/sum(tt);

      Text{t} = [Text{t} sprintf(', new double[]')];
      Text{t} = [Text{t} subWrite(fid, tt)];

      Text{t} = [Text{t} sprintf(', new double[]')];
      Text{t} = [Text{t} subWrite(fid, Node(n).leftLengthDist)];

      Text{t} = [Text{t} sprintf(', new double[]')];
      Text{t} = [Text{t} subWrite(fid, Node(n).leftLetterDist)];

      Text{t} = [Text{t} sprintf(', new double[]')];
      Text{t} = [Text{t} subWrite(fid, Node(n).rightLengthDist)];

      Text{t} = [Text{t} sprintf(', new double[]')];
      Text{t} = [Text{t} subWrite(fid, Node(n).rightLetterDist)];
 
      Text{t} = [Text{t} sprintf(');')];

      LI = Node(n).LeftNTIndex;
      RI = Node(n).RightNTIndex;

%[length(LI) length(RI)]

      Text{t} = [Text{t} sprintf(' // %s%s - %s%s %s\n', File.NT(LI).Base, File.NT(LI).Number, File.NT(RI).Base, File.NT(RI).Number, zEdgeText(File.Edge(LI,RI)))];

    case 'Cluster'

      % initialize Cluster node

      Text{t} = [Text{t} sprintf('current = new ClusterNode(current, ')];
      Text{t} = [Text{t} sprintf('%0.6f, ', Node(n).Delete)];
      Text{t} = [Text{t} sprintf('%d, ', length(Node(n).Left))];
      Text{t} = [Text{t} sprintf('%d);\n', length(Node(n).Right))];

      % add interactions between bases

      Indices = [Node(n).LeftIndex(Node(n).Left) Node(n).RightIndex(Node(n).Right)];
      for i = 1:length(Node(n).IBases(:,1)),
        Text{t} = [Text{t} sprintf('((ClusterNode)current).addInteraction(')];
        Text{t} = [Text{t} sprintf('%3d, %3d, new double[][]', Node(n).IBases(i,1), Node(n).IBases(i,2))];
        p = exp(Node(n).Score(:,:,i));                % exponentiate
        p = p / sum(sum(p));                          % normalize
        Text{t} = [Text{t} subWriteArray(fid, p)];
        Text{t} = [Text{t} sprintf(');')];
        LI = Indices(Node(n).IBases(i,1));
        RI = Indices(Node(n).IBases(i,2));
        Text{t} = [Text{t} sprintf(' // %s%s - %s%s %s', File.NT(LI).Base, File.NT(LI).Number, File.NT(RI).Base, File.NT(RI).Number, zEdgeText(File.Edge(LI,RI)))];
        Text{t} = [Text{t} sprintf('\n')];
      end

      % add insertions between bases

      % later!

      Text{t} = [Text{t} sprintf('((ClusterNode)current).normalize();\n')];

    case 'Junction'
      Text{t} = [Text{t} sprintf('current = new JunctionNode(current,')];
      Text{t} = [Text{t} sprintf('%d',length(Node(n).nextnode))];
      Text{t} = [Text{t} sprintf(');\n')];

    case 'Hairpin'
      Text{t} = [Text{t} sprintf('current = new HairpinNode(current, "%s");\n', Node(n).subtype)];

  end
end

for i = 1:length(Text)
  fprintf('%s',Text{i});
end


%fclose(fid)

function [newText] = subWrite(fid, a)

  newText = '';
  newText = [newText sprintf('{')];
  for i = 1:length(a)-1,
    newText = [newText sprintf('%0.6f,', a(i))];
  end
  newText = [newText sprintf('%0.6f}', a(end))];

function [newText] = subWriteArray(fid, a)

  newText = '';
  newText = [newText sprintf('{')];
  for i = 1:4,
    newText = [newText sprintf('{')];
    for j = 1:3,
      newText = [newText sprintf('%0.6f,', a(i,j))];
    end
    newText = [newText sprintf('%0.6f}',a(i,4))];
    if i < 4,
      newText = [newText sprintf(',')];
    end
  end
  newText = [newText sprintf('}')];

