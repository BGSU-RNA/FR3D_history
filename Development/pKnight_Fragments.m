
File = zAddNTData('2avy');

% remove a cWW-cWW triple

i1 = zIndexLookup(File(1),{'454'});
i2 = zIndexLookup(File(1),{'478'});
File(1).Edge(i1,i2) = 0;
File(1).Edge(i2,i1) = 0;

% remove a triple that is absent in 1j5e, just for comparison

i1 = zIndexLookup(File(1),{'69'});
i2 = zIndexLookup(File(1),{'100'});
File(1).Edge(i1,i2) = 0;
File(1).Edge(i2,i1) = 0;

if 0 > 1,
  F = File;
else
  F = xAnnotateWithKnownMotifs(File(1),1);
end
Node = pMakeNodes(F,1,'107','338');
Node = pShiftNodeIndices(Node,Node(1).LeftIndex-1);
Node(1).leftLengthDist = ones(1,40)/40;
pWriteJavaNodeFile(F,Node,4,'16S_2avy_107-338.txt');
fprintf('Running JAR3D\n');
JAR3D('Knight-16S-soil-segments.fasta','16S_2avy_107-338.txt',10,2);
