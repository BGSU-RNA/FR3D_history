
Filenames = zReadPDBList('Nonredundant_2009-05-14_list',1);

for f = 1:length(Filenames),
  copyfile(['PrecomputedData\' Filenames{f} '.mat'],['NRData\' Filenames{f} '.mat'])
end
