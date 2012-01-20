% zListNonredundantSet lists PDB ID's and chains that are meant to be non-redundant

Verbose = 1;

NRList = 'Nonredundant_4A_2010-05-19_list';

Filenames = zReadPDBList(NRList,1);

File = zAddNTData(NRList,0,[],1);

fid = fopen([pwd filesep 'Web' filesep 'AnalyzedStructures' filesep NRList '.txt'],'w');

Vers = num2str(File(1).ClassVersion);
fprintf(fid,'# PDB_ID_FR3D_Version_%s\n',Vers);

fprintf(fid,'PDB_ID\tChain(s)\n');

zBestChains

fclose(fid);



















