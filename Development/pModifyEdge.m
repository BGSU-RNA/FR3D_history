
function [F] = pModifyEdge(F,N1,N2,EdgeText)

i1 = zIndexLookup(F,{N1});
i2 = zIndexLookup(F,{N2});

F.Edge(i1,i2) =  zEdgeFromEdgeText(EdgeText);
F.Edge(i2,i1) = -zEdgeFromEdgeText(EdgeText);

