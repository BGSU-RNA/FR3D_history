
File = zAddNTData('NonRedundant_list',1);      % reclassify
zFindExemplars
File = zUpdateDistanceToExemplars(File);

zExemplarTable(1,3.5)
zExemplarTable(2,3.5)
zExemplarTable(3,3.5)
zExemplarTable(4,3.5)
zExemplarTable(5,3.5)
zExemplarTable(6,3.5)
zExemplarTable(7,3.5)
zExemplarTable(8,3.5)
zExemplarTable(9,3.5)
zExemplarTable(10,3.5)
zExemplarTable(11,3.5)
zExemplarTable(12,3.5)
zExemplarTable(13,3.5)
zExemplarTable([1:13],3.5)
zExemplarTable([1 5],3.5)
zExemplarTable(1:6,3.5)
zExemplarTable(7:12,3.5)

break


for f = 1:length(File),
  if isfield(File(f).Header,'Expdata'),
    disp(File(f).Header.Expdata)
  end
end

for f = 1:length(File),
  if isfield(File(f).Header,'Resolution'),
    disp(File(f).Header.Resolution)
  end
end


