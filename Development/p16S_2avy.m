% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

if exist('File'),
  [File,FIndex] = zAddNTData('2avy',0,File);
else
  [File,FIndex] = zAddNTData('2avy',0);
end

F = File(FIndex);

nMin = 1;
nMax = length(File(FIndex).NT);

F = File(FIndex);

%F.Edge = F.Edge .* (F.Range <= 10);

% zSecondaryStructure(F,nMin,nMax);

Node = pMakeNodes(F,nMin,nMax);

%Node = pShiftNodeIndices(Node,nMin);

S = pTheoreticalAlignment(Node,1);

M = strrep(S{2},'-','');
M = strrep(M,'(','');
M = strrep(M,')','');
M = strrep(M,'<','');
M = strrep(M,'>','');
M = strrep(M,'{','');
M = strrep(M,'}','');

fprintf('Diagnostic:  Compare bases from structure (first line) with bases in experimental alignment:\n');
fprintf('%s seq: %s\n', File(FIndex).Filename, cat(2,File(FIndex).NT(nMin:nMax).Base));
fprintf('pTheoret: %s\n', M);
fprintf('\n');

fprintf('Diagnostic:  Compare experimental alignment with Java alignment.\n');
fprintf('Experimental alignment for %s:\n', File(FIndex).Filename);
fprintf('%s\n', S{1});
fprintf('%s\n', S{2});
fprintf('\n');

pWriteJavaNodeFile(File(FIndex),Node,4,'16S_from_2avy.txt');

figure(1)
clf
zCircularDiagram(File,0.2,[1 1 1 1 0 0 0]);
saveas(gcf,'2avy_circular_basepairs.pdf','pdf');

figure(2)
clf
FF = File;
Node(1).Edge(File.NumNT,File.NumNT) = 0;
FF.Edge = sparse(Node(1).Edge + Node(1).Edge');
zCircularDiagram(FF,0.2,[1 1 1 1 0 0 0]);
saveas(gcf,'2avy_circular_SCFG_basepairs.pdf','pdf');
