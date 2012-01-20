a = dir(['PDBFiles' filesep '*.pdb']);

for i=1:length(a),
  name = a(i).name;
  if strcmp('pdb',name(1:3)),
    fprintf('ren %s %s\n',name,name(4:end));
  end
end
