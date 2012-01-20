% pMakeNodesDiagnostics runs some simple diagnostics of the model

% ------------------------------ Compare basepairs in the model structure

FF = File(1);
Node(1).Edge(FF.NumNT,FF.NumNT) = 0;
FF.Edge = sparse(Node(1).Edge + Node(1).Edge');
FF.Filename = [File(1).Filename '_JAR3D'];
FF.BasePhosphate = 0*FF.BasePhosphate;
zSummarizeInteractions([File(1) FF],1);

% ------------------------------ List nested basepairs not used in model

[i,j,k] = find(triu(FF.Edge == 0) .* (abs(File(1).Edge) < 13) .* (abs(File(1).Edge) > 0) .* (File(1).Crossing == 0));
for a = 1:length(i),
  fprintf('Nested basepair %s%4s - %s%4s %4s is missed by pMakeNodes\n', File(1).NT(i(a)).Base, File(1).NT(i(a)).Number, File(1).NT(j(a)).Base, File(1).NT(j(a)).Number, zEdgeText(File(1).Edge(i(a),j(a))));
end

% ------------------------------ List basepairs that extend across junctions

kk = Node(1).JunctionDeletion;
kk = [kk abs(kk)];
[y,i] = sortrows(kk,2);
for i = 1:length(kk),
  fprintf('Removed a %4s basepair across a junction.\n', zEdgeText(kk(i,1)));
end
