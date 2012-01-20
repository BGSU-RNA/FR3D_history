
function Colors = zColorNucleotides(File,NTList,Color,Colors)

Indices = zIndexLookup(File,NTList);

for k = Indices,
  Colors(k,:) = Color;
end
