
% File = zAddNTData({'2avy','1j5e'});
for f = 1:length(File),
  c = cat(1,File(f).NT(1:File(f).NumNT).Center); % nucleotide centers
  File(f).Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
end

Verbose = 1;

% -------------------------------------- Load Ryan's 3D to 3D alignment

load 16SalignmentInd

a = 1;
Al(a).ModelStructure = ind2';
Al(a).InferStructure = ind1';
Al(a).Name           = 'Ryan 3D to 3D';

% -------------------------------------- Load Jesse's 3D to 3D alignment

[File1,i1,File2,i2] = zReadJesseAlignments('16S');

a = 2;
Al(a).ModelStructure = i1';
Al(a).InferStructure = i2';
Al(a).Name           = 'Jesse 3D to 3D';

% ------------------------------ JAR3D alignment of 1j5e seq to 2avy structure

% >E. coli sequence from 2avy -1661.9402926282287 

ES = 'UGAA---GAGU<UUGAUCAUG>GCUC--AGAUUGAACGC--{UGGCGGCAGG-{CCUAA{CA{CAUGCAA--GUC{GAACG{GUAACAGGAAGAAGC<UU>-GCUUCUUUGCUGACG}AGU}GGCGGACG-GGUGAGU-AAU-GUCUGGGAA-ACUGCCUGAUGGA--GGGGGAUAACUACUGG<AA>ACGGUAGCUAAUACCGC--AUAACGUCG<CA>AGACCAA-A-GAGGGGGA--CCU<UC>GGGCCUCUU-GCCAUCGGAUGUGCCCAGAU--GGG{AUUAGCUAGUAGGUGGG<GUAACGG>CUCACCUAGGCGACGAU}CCC-UA-GCUGGUCUG<AG>AGGAUGACCAGC-CACA----CUGGAACU<GAG>ACACGGUCCAG--ACUCCU<AC>GGGAGG--CAGCAG}UG}GGG-AAU}--AUUGC{ACAAUGGGCG<CA>AGCCUGAU}GCAG-CCAUGCCGC}G-U-GUAUG{AAGAAGGCCU<UC>GGGUUGUAA}AGUAC--{UUUCAGCGGGGAGGAAG{GGAGUAAAG<UUAAUAC>CUUUGCUCAU}UGACGUUACCCGCAGAAG}-AA---GC{ACCGGCUAACUCCG{UGCCAGC<AGCC>GCGGUAAUA}CGGAG-GGU}GC-AAGCGUUAAUCGG-AAUU--ACUG{GGCGUAAAGC--GCACGCAGGCG-GUUUGUU{AAG{UCAGAUGUGAAAUCCCCGGGC<UCA>ACCUGGGAACUGCAUCUGAUA}CU}GGCAAGC--{UUGAGUCUCGUAGA{GGGGGGU-AGAAUUCCAGGUGUAGCGG<UGAA>AUGCGUAGAGAUCUGGAGGAAUACCG-GUGGC<GAAG>GCG-GCCCCCU}GGACGAAGACUGAC}-GCUCAGGUGCG-AAAGCGUGGGGAGCAAACAGG<AUUAGAUAC>CCUGGUAGUCCACGCCGU--AAACGA{UGUCGACU-UGGAGGUUGUGCC<CUUGA>GGCGUGGCUUCCGG-AGC<UAAC>GCG-UUAAGUCGAC}-CGCC}U--GGGGA{GUACGGCCG<CA>AGGUUAA}AACUCAAAUGAAU----{UG-ACGGGGGCCCGC{ACAA-GCGGUG-GAGCAUGUGGUU{UAAUU-CGA<UGCA>ACGCGAAGAACCUUAC-CUGGUCU-{UGACAUCCACGGAA--GUUUUCA<GAGA>UGAGAAU--GU-GCC<UUCG>GGA---ACCGUGAGAC}--{AGGUGCU{GCAUG{GCUGUCGUCA-GCUCGU-GUU<GUGA>AAUGU-UGGGU<UAAG>UCCCG-CAACGAGCGCAAC-CCUUA-UCCUUUGUUGCCAG{CGG<UCCG>GCCG}GGAACUCAAAGGA--{GACUGCCAGUG<AUA>A-ACUGGAGGAAGG}-UGGGG--AUGACGUCAAGUCAU}CAU}GGC-CCUUAC}-GACCAG--G}GCUACACACGUGCUA--C-AAUGGCGCAUACAAAGAGAA{GCGACCUC<GCGA>GAGCAAGCG}GACCUCAUAAAGUGCGUCGUAGUC-CGGAUUGGAGUCU<GCA>ACUCGACUCCAUGAAG-UCG-GAAUCGCU-A{GUAAUCGUGGAU<CAGA>A-UGCCACGGUGA}A-UACGUUCC}CGGGCCUUGUACACA}--CCGCCCGUCACACCAUGGGAGUGGGUUGCAAAAGAAGUAGG{UAGCUUAACCU<UC>GGGAGGGCG}CUUACCACUUUGUGAUUCAUGACUGGGGUGAAGUCGUAACAAG-GU-AACCGUAGG<GGAA>CCUGCGGUU--G-----GAUCA';

% >T.th. sequence from 1j5e -1798.11196029636 

TS = 'UGGA---GAGU<UUGAUCCUG>GCUC--AGGGUGAACGC--{UGGCGGCGUG-{CCUAA{GA{CAUGC----AAG{U-CGU{GCG-GGCCGCGGGGU<UU>UACUCCGUGGUCAGCG}GCG}GAC----G-GGUGAGU-AAC-GCGUGGGUGACCUACCCGGAAGA--GGGGGACAACCCGGGG<AA>ACUCGGGCUAAUCCCCC--AUGUGGACC<CG>CCCCUUG---GGGUGUGUCCAAA<GG>GCUUUGCCC-GCUUCCGGAUGGGCCCGCGU--CCC{AUCAGCUAGUUGGUGGG<GUAAUGG>CCCACCAAGGCGACGAC}GGG-UA-GCCGGUCUG<AG>AGGAUGGCCGGC-CACA----GGGGCACU<GAG>ACACGGGCCCC--ACUCCU<AC>GGGAGG--CAGCAG}UU}AGGAAUC}--UU-CC{GCAAUGGGCG<CA>AGCCUGAC}GGAG-CGACGCCGC}U-U-GGAGG{AAGAAGCCCU<UC>GGGGUGUAA}ACUCC--{UG-AACCCGGGAC----{GAAACCC-C<CGACGAG>GGG-ACUGAC}----GGUACCGGGGUAAU}--A---GC{GCCGGCCAACUCCG{UGCCAGC<AGCC>GCGGUAAUA}CGGAG-GGC}GC-GAGCGUUACCCGG-AUUC--ACUG{GGCGUAAAGG--GCGUGUAGGCG-GCCUGGG{GCG{UCCCAUGUGAAAGACCACGGC<UCA>ACCGUGGGGGAGCGUGGGAUA}CG}CUCAGGC--{UAGACGGUGGGAGA{GGGUGGU-GGAAUUCCCGGAGUAGCGG<UGAA>AUGCGCAGAUACCGGGAGGAACGCCG-AUGGC<GAAG>GCA-GCCACCU}GGUCCACCCGUGAC}-GCUGAGGCGCG-AAAGCGUGGGGAGCAAACCGG<AUUAGAUAC>CCGGGUAGUCCACGCCCU--AAAC--{GAUGCGC--GCUAGGUCUCUGG<GUCUC>CUGGGGGCCGAAGC-UAA<CGCG>UUA----AGCGCGC}-CGCC}U--GGGGA{GUACGGCCG<CA>AGGCUGA}AACUCAAAGGAAU----{UG-ACGGGGGCCCGC{ACAA-GCGGUG-GAGCAUGUGGUU{UAAUU-CGA<AGCA>ACGCGAAGAACCUUAC-CAGGCCU-{UGACAUGCUAGGGAA-CCCGGGU<GAAA>GCCUGGGGUGC-CCC<GCGA>GGG-GAGCCCUAGCAC}--{AGGUGCU{GCAUG{GCCGUCGUCA-GCUCGU-GCC<GUGA>GGUGU-UGGGU<UAAG>UCCCG-CAACGAGCGCAAC-CCCCG-CCGUUAGUUGCCAG{CGG<UUCG>GCCG}GGCACUCUAACGG--{GACUGCCCGC-<GAA>-AGCGGGAGGAAGG}--AGGG-GACGACGUCUGGUCAG}CAU}GGC-CCUUAC}-GGCCUG--G}GCGACACACGUGCUA--C-AAUGCCCACUACAAAGCGAU{GCCACCCG<GCAA>CGGGGAGCU}AAUCGCAAAAAGGUGGGCCCAGUU-CGGAUUGGGGUCU<GCA>ACCCGACCCCAUGAAG-CCG-GAAUCGCU-A{GUAAUCGCGGAU<CAGC>CAUGCCGCGGUGA}A-UACGUUCC}CGGGCCUUGUACACA}--CCGCCCGUCACGCCAUGGGAGCGGGCUCUACCCGAAGUCGC{CGGGA---GCC<UA>CGG--GCAG}GCGCCGAGGGUAGGGCCCGUGACUGGGGCGAAGUCGUAACAAG-GU-AGCUGUACC<GGAA>GGUGCGGCU-GGAUCACUUUCU';

[i1,i2] = pGetAlignedIndices(ES,TS);

a = 3;
Al(a).ModelStructure = i1;
Al(a).InferStructure = i2;
Al(a).Name           = 'JAR3D Version 1';

% --------------------------------------- Needleman-Wunsch alignment

[matches,align1,align2,s1,s2] = zNeedlemanWunsch(cat(2,File(1).NT.Base),cat(2,File(2).NT.Base));

a = 4;
Al(a).ModelStructure = align1;
Al(a).InferStructure = align2;
Al(a).Name           = 'Needleman-Wunsch';

% --------------------------------------- Calculations for each alignment

for a = 4:length(Al),
  Al(a).Matrix = sparse(Al(a).ModelStructure,Al(a).InferStructure,ones(1,length(Al(a).ModelStructure)));
  Al(a).Matrix(length(File(1).NT),length(File(2).NT)) = 0;
  if Verbose > 1,
    Al(a).Tally = zAlignmentDiagram(File,Al(a).ModelStructure,Al(a).InferStructure,1);
  else
    Al(a).Tally = zAlignmentDiagram(File,Al(a).ModelStructure,Al(a).InferStructure,0);
  end
end

% --------------------------------------- Superimpose local neighborhoods

for a = 1:length(Al),
  Al(a).Discrep = zHistogramDiscrepanciesInAlignment(File,Al(a).ModelStructure,Al(a).InferStructure);
end

figure(10)
clf

m = 3;

for a = 1:length(Al),
  subplot(2,2,a)
  n = hist(min(Al(a).Discrep,m),-0.025+(0:0.05:m));
  hist(min(Al(a).Discrep,m),-0.025+(0:0.05:m))
  title(Al(a).Name);
  axis([0 m 0 max(n)*1.1])
end

saveas(gcf,'Histogram Discrepancies in 16S Alignments.pdf','pdf')

% --------------------------------------- Compare to Alignment 1

for a = 1:length(Al),
  Agree  = sum(sum(Al(a).Matrix .* Al(1).Matrix == 1));
  Missed = sum(sum(Al(1).Matrix > Al(a).Matrix));
  Extra  = sum(sum(Al(a).Matrix > Al(1).Matrix));

  T{ 1,a+1} = Al(a).Name;
  T{ 2,a+1} = length(Al(a).ModelStructure);
  for v = 1:6,
    T{v+2,a+1} = Al(a).Tally(v);
  end
  T{ 9,a+1} = Agree;
  T{10,a+1} = Missed;
  T{11,a+1} = Extra;
  T{12,a+1} = mean(Al(a).Discrep);
  T{13,a+1} = median(Al(a).Discrep);
end

T{1,1} = [File(1).Filename ' as the model, ' File(2).Filename ' as the unknown structure'];
T{2,1} = 'Number aligned';
T{3,1} = ['Nested cWW correctly inferred ' File(2).Filename];
T{4,1} = ['Nested non-cWW correctly inferred in ' File(2).Filename];
T{5,1} = ['Non-nested cWW correctly inferred in ' File(2).Filename];
T{6,1} = ['Non-nested non-cWW correctly inferred in ' File(2).Filename];
T{7,1} = ['Stacking correctly inferred in ' File(2).Filename];
T{8,1} = ['Base-phosphate correctly inferred in ' File(2).Filename];
T{9,1} = ['Number of correspondences agreeing with ' Al(1).Name];
T{10,1} = ['Number of correspondences missing, compared to ' Al(1).Name];
T{11,1} = ['Number of correspondences extra, compared to ' Al(1).Name];
T{12,1} = 'Mean geometric discrepancy at 8 Angstroms';
T{13,1} = 'Median geometric discrepancy at 8 Angstroms';

xlswrite('16S_alignment_comparison.xls',T);

T

% ---------------------------------------- Visually compare to Alignment 1

if Verbose > 1,

for a = 2:length(Al),
  zCompareAlignment(File,Al(1).ModelStructure,Al(1).InferStructure,Al(a).ModelStructure,Al(a).InferStructure,Al(1).Name,Al(a).Name);
  Titl = ['Agreement between ' Al(a).Name ' and ' Al(1).Name];
%  title(Titl);
  saveas(gcf,[strrep(Titl,' ','_') '.pdf'],'pdf');
end



for a = 2:length(Al),
  [i1,i2,s] = find(Al(a).Matrix .* Al(1).Matrix == 1);
  figure(1)
  clf
  zShowAlignment(File,i1,i2);
  Titl = ['Agreement between ' Al(a).Name ' and ' Al(1).Name];
  title(Titl);
  saveas(gcf,[strrep(Titl,' ','_') '.pdf'],'pdf');

  [i1,i2,s] = find(Al(1).Matrix > Al(a).Matrix);
  figure(2)
  clf
  zShowAlignment(File,i1,i2);
  Titl = ['Correspondences missing in ' Al(a).Name ' compared to ' Al(1).Name];
  title(Titl);
  saveas(gcf,[strrep(Titl,' ','_') '.pdf'],'pdf');

  [i1,i2,s] = find(Al(a).Matrix > Al(1).Matrix);
  figure(3)
  clf
  zShowAlignment(File,i1,i2);
  Titl = ['Extra correspondences in ' Al(a).Name ' compared to ' Al(1).Name];
  title(Titl);
  saveas(gcf,[strrep(Titl,' ','_') '.pdf'],'pdf');
end

end