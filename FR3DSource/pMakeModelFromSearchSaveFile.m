% pMakeModelFromSearchSaveFile(Search) creates an SCFG/MRF Node variable corresponding to the model in Search

function [Node] = pMakeModelFromSearchSaveFile(Search)

if strcmp(class(Search),'char'),
  load(['SearchSaveFiles' filesep Search,'Search','-mat');
end

Node = pMakeNodes(Search.Query,1);

a = Node(1).LeftIndex;
b = Node(1).RightIndex-length(Node);

for n = 1:length(Node),
  Node(n).LeftIndex  = Node(n).LeftIndex - a + 1;
  Node(n).RightIndex = Node(n).RightIndex - b + 1;
  Node(n).RightIndex = max(Node(n).RightIndex,max(Node(n).LeftIndex)+1);
  Node(n).MiddleIndex = Node(n).MiddleIndex - a + 1;
end

% the next program uses a strange filename extension ::::::::

pWriteJavaNodeFile(Search.Query,Node,5);
