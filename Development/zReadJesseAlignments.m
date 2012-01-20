
function [File,i1,i2] = zReadJesseAlignment(Molecule)

switch Molecule

case '16S',  % ----------------------------------------------- 16S alignment

  [n,t,r] = xlsread('Stombaugh_et_al_Sup_Mat_S9.xls','16S_3D_Struct_Aln');

  cols = [5 7 18 20];
  startrow = 3;
  Chain1 = 'A';
  Chain2 = 'A';

  File1 = zAddNTData('2avy');
  File2 = zAddNTData('1j5e');

case '5S',  % ----------------------------------------------- 5S alignment

  [n,t,r] = xlsread('Stombaugh_et_al_Sup_Mat_S9.xls','5S_3D_Struct_Aln');

  cols = [5 7 18 20];
  startrow = 3;
  Chain1 = 'A';
  Chain2 = 'B';

  File1 = zAddNTData('2aw4');
  File2 = zAddNTData('2j01');

case '23S',  % ----------------------------------------------- 23S alignment

  [n,t,r] = xlsread('Stombaugh_et_al_Sup_Mat_S9.xls','23S_3D_Struct_Aln');

  cols = [5 7 18 20];
  startrow = 3;
  Chain1 = 'B';
  Chain2 = 'A';

  File1 = zAddNTData('2aw4');
  File2 = zAddNTData('2j01');

end


[i1,i2] = LookUpIndices(File1,File2,Chain1,Chain2,r,cols,startrow);

rWriteAlignmentBPComparison(File1,1:length(File1.NT),File2,1:length(File2.NT),i1,i2);


% ---------------------------------------------- LookUpIndices routine

function [i1,i2] = LookUpIndices(File1,File2,Chain1,Chain2,r,cols,startrow);

for a = startrow:length(r(:,1)),
  if ~isempty(r{a,cols(1)}),                         % row has an NT number
    i = zIndexLookup(File1,r{a,cols(1)},Chain1);
    j = zIndexLookup(File1,r{a,cols(2)},Chain1);
    k = zIndexLookup(File2,r{a,cols(3)},Chain2);
    m = zIndexLookup(File2,r{a,cols(4)},Chain2);
    indices = [indices; [i k]; [j m]];               % append these matches
  end
end

indices = unique(indices);                           % remove redundant aligns

i1 = indices(:,1);
i2 = indices(:,2);
