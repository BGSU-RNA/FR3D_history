% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

File = zAddNTData({'2avy','1j5e'});

Node = pMakeNodes(File(1));
pWriteJavaNodeFile(File(1),Node,4,'16S_from_2avy.txt');

F = xAnnotateWithKnownMotifs(File(1),1);
NodeM = pMakeNodes(F);
pWriteJavaNodeFile(F,NodeM,4,'16S_from_2avy_with_motifs.txt');

JAR3D('16S_sequence_from_2avy_1j5e.fasta','16S_from_2AVY_with_motifs.txt');

% ------------------------------------------------ Start diagnostics

S = pTheoreticalAlignment(Node,1);

M = strrep(S{2},'-','');
M = strrep(M,'(','');
M = strrep(M,')','');
M = strrep(M,'<','');
M = strrep(M,'>','');
M = strrep(M,'{','');
M = strrep(M,'}','');

fprintf('Diagnostic:  Compare bases from structure (first line) with bases in experimental alignment:\n');
fprintf('%s seq: %s\n', File(1).Filename, cat(2,File(1).NT(nMin:nMax).Base));
fprintf('pTheoret: %s\n', M);
fprintf('\n');

fprintf('Diagnostic:  Compare experimental alignment with Java alignment.\n');

clear Sequence
Sequence(1).Organism = 'E. coli sequence from 2avy';
Sequence(1).Fasta    = cat(2,File(1).NT.Base);
pWriteFASTA(Sequence,['sequences' filesep '16S_sequence_from_2avy.fasta']);

JAR3D('16S_sequence_from_2avy.fasta','16S_from_2AVY.txt');

fprintf('%s\n', S{2});
fprintf('%s\n', S{1});
fprintf('Experimental alignment for %s is above.\n', File(1).Filename);
fprintf('\n');

fprintf('Diagnostic:  Compare alignment of T.th. sequence to E.coli structure.\n');

Sequence(2).Organism = 'T.th. sequence from 1j5e';
Sequence(2).Fasta    = '';
S                    = cat(2,File(2).NT.Base);
for i = 1:length(S),
  Sequence(2).Fasta = [Sequence(2).Fasta S(i)];
  if mod(i,60) == 0,
    Sequence(2).Fasta = [Sequence(2).Fasta '-'];
  end
end
        
pWriteFASTA(Sequence,['sequences' filesep '16S_sequence_from_2avy_1j5e.fasta']);

JAR3D('16S_sequence_from_2avy_1j5e.fasta','16S_from_2AVY.txt');

% The following diagnostic was lifted from pInferVsActual3DStructure.m

JH = '[---[[[((((<FFFFFFFFF>))))][(((((-((((([[{I-IIIIII(([{I--I({I({II(((--[[((({I-III{I(--(((((((((((<FF>-))))))))-))))II}III})))----][((((---[(--[((((((((--(((--(((((((([[(-((((((((((((((<FF>)))))))))))-)))-)][((---((((<FF>))))-))]][((((((----(((<FF>)))))))))]))))))))-)))))))))))][((({I-III((--((((((((<FFFFFFF>))))))))))II--III})))])][(((((((((<FF>)))))---))))]))))]][[((((((((<FFF>)-)--))))))][((((((<FF>))))))]]--)))I})I})II---I}][((-(({I--I((((((<FF>)))))-)I}))))]))IIIIIII}][([((((({I--I[(((((<FF>)))))]III})))))][{II-I[(((((((((-(-{I((((((-(<FFFFFFF>)))-))))II})--))-)))))))-]III}]-)]][(({III----I[((((({I---III<FFFF>III---I-I})))))]III}))]--)))))))))-)]---][[(---{III---II(-[[((((((((((-[((((((({II({I(((((((-----((((((((<FFF>))))))))---)))))))III})I})))))))][{I-I(((((((((--{I(((((([((-((((((((((((((-(<FFFF>)-))---))))))))))))-))---][(--((<FFFF>)))]))))))I}-)))))-))))III}]))))))))))][(--((((((---(((((((((<FFFFFFFFF>)))))))))))))))--)]][-----{IIIIIII([(((((((((((((<FFFFF>)))))))))))))][(((<FFFF>)))]--)IIIIIII}])III})][((((({II(-(((((<FF>))))))I})))))----]--]][[[{II[(((-(((((((-{I(-([(((((([(((((((((((-{I-I[-[(((<FFFF>)--))----------][((((((-[{II--II((((((---[(((((((<FFFF>)))))))---][(((<FFFF>)))]---)))))-)II}][{III[((({I-III{I((((((---[(((((([(((<FFFF>)))-][(((((<FFFF>)))))]--))))))----][((((-[(((((((--(((--{III<FFFF>IIII})))---)))))))][{III-I((((((<FFF>)-)))))IIII-II}]-))))]---)--)))-))I--I}III})))]III-II}]))))))]]I}---))-)))))))))][([(--(((((((((--((((--{I(((-(((<FFFF>)))-)))II}---))))--)))))))))----)][(((((((((((((<FFF>)---))))))))))))]--)]-))))))][({II((((((((((<FFFF>)--)))))))))I})])--)I--I}))))))))))]II-I}][(-(-(((-(--((((((((((((((((((((((((-((((({II((((--(((<FF>)))-))))I})))))))))))))))))))))))))))))---)-))-)-)--)]-][(((((((((<FFFF>)))))))))]-]---------]';

% >E. coli sequence from 2avy -1661.9402926282287 

ES = 'UGAA---GAGU<UUGAUCAUG>GCUC--AGAUUGAACGC--{UGGCGGCAGG-{CCUAA{CA{CAUGCAA--GUC{GAACG{GUAACAGGAAGAAGC<UU>-GCUUCUUUGCUGACG}AGU}GGCGGACG-GGUGAGU-AAU-GUCUGGGAA-ACUGCCUGAUGGA--GGGGGAUAACUACUGG<AA>ACGGUAGCUAAUACCGC--AUAACGUCG<CA>AGACCAA-A-GAGGGGGA--CCU<UC>GGGCCUCUU-GCCAUCGGAUGUGCCCAGAU--GGG{AUUAGCUAGUAGGUGGG<GUAACGG>CUCACCUAGGCGACGAU}CCC-UA-GCUGGUCUG<AG>AGGAUGACCAGC-CACA----CUGGAACU<GAG>ACACGGUCCAG--ACUCCU<AC>GGGAGG--CAGCAG}UG}GGG-AAU}--AUUGC{ACAAUGGGCG<CA>AGCCUGAU}GCAG-CCAUGCCGC}G-U-GUAUG{AAGAAGGCCU<UC>GGGUUGUAA}AGUAC--{UUUCAGCGGGGAGGAAG{GGAGUAAAG<UUAAUAC>CUUUGCUCAU}UGACGUUACCCGCAGAAG}-AA---GC{ACCGGCUAACUCCG{UGCCAGC<AGCC>GCGGUAAUA}CGGAG-GGU}GC-AAGCGUUAAUCGG-AAUU--ACUG{GGCGUAAAGC--GCACGCAGGCG-GUUUGUU{AAG{UCAGAUGUGAAAUCCCCGGGC<UCA>ACCUGGGAACUGCAUCUGAUA}CU}GGCAAGC--{UUGAGUCUCGUAGA{GGGGGGU-AGAAUUCCAGGUGUAGCGG<UGAA>AUGCGUAGAGAUCUGGAGGAAUACCG-GUGGC<GAAG>GCG-GCCCCCU}GGACGAAGACUGAC}-GCUCAGGUGCG-AAAGCGUGGGGAGCAAACAGG<AUUAGAUAC>CCUGGUAGUCCACGCCGU--AAACGA{UGUCGACU-UGGAGGUUGUGCC<CUUGA>GGCGUGGCUUCCGG-AGC<UAAC>GCG-UUAAGUCGAC}-CGCC}U--GGGGA{GUACGGCCG<CA>AGGUUAA}AACUCAAAUGAAU----{UG-ACGGGGGCCCGC{ACAA-GCGGUG-GAGCAUGUGGUU{UAAUU-CGA<UGCA>ACGCGAAGAACCUUAC-CUGGUCU-{UGACAUCCACGGAA--GUUUUCA<GAGA>UGAGAAU--GU-GCC<UUCG>GGA---ACCGUGAGAC}--{AGGUGCU{GCAUG{GCUGUCGUCA-GCUCGU-GUU<GUGA>AAUGU-UGGGU<UAAG>UCCCG-CAACGAGCGCAAC-CCUUA-UCCUUUGUUGCCAG{CGG<UCCG>GCCG}GGAACUCAAAGGA--{GACUGCCAGUG<AUA>A-ACUGGAGGAAGG}-UGGGG--AUGACGUCAAGUCAU}CAU}GGC-CCUUAC}-GACCAG--G}GCUACACACGUGCUA--C-AAUGGCGCAUACAAAGAGAA{GCGACCUC<GCGA>GAGCAAGCG}GACCUCAUAAAGUGCGUCGUAGUC-CGGAUUGGAGUCU<GCA>ACUCGACUCCAUGAAG-UCG-GAAUCGCU-A{GUAAUCGUGGAU<CAGA>A-UGCCACGGUGA}A-UACGUUCC}CGGGCCUUGUACACA}--CCGCCCGUCACACCAUGGGAGUGGGUUGCAAAAGAAGUAGG{UAGCUUAACCU<UC>GGGAGGGCG}CUUACCACUUUGUGAUUCAUGACUGGGGUGAAGUCGUAACAAG-GU-AACCGUAGG<GGAA>CCUGCGGUU--G-----GAUCA';

% >T.th. sequence from 1j5e -1798.11196029636 

TS = 'UGGA---GAGU<UUGAUCCUG>GCUC--AGGGUGAACGC--{UGGCGGCGUG-{CCUAA{GA{CAUGC----AAG{U-CGU{GCG-GGCCGCGGGGU<UU>UACUCCGUGGUCAGCG}GCG}GAC----G-GGUGAGU-AAC-GCGUGGGUGACCUACCCGGAAGA--GGGGGACAACCCGGGG<AA>ACUCGGGCUAAUCCCCC--AUGUGGACC<CG>CCCCUUG---GGGUGUGUCCAAA<GG>GCUUUGCCC-GCUUCCGGAUGGGCCCGCGU--CCC{AUCAGCUAGUUGGUGGG<GUAAUGG>CCCACCAAGGCGACGAC}GGG-UA-GCCGGUCUG<AG>AGGAUGGCCGGC-CACA----GGGGCACU<GAG>ACACGGGCCCC--ACUCCU<AC>GGGAGG--CAGCAG}UU}AGGAAUC}--UU-CC{GCAAUGGGCG<CA>AGCCUGAC}GGAG-CGACGCCGC}U-U-GGAGG{AAGAAGCCCU<UC>GGGGUGUAA}ACUCC--{UG-AACCCGGGAC----{GAAACCC-C<CGACGAG>GGG-ACUGAC}----GGUACCGGGGUAAU}--A---GC{GCCGGCCAACUCCG{UGCCAGC<AGCC>GCGGUAAUA}CGGAG-GGC}GC-GAGCGUUACCCGG-AUUC--ACUG{GGCGUAAAGG--GCGUGUAGGCG-GCCUGGG{GCG{UCCCAUGUGAAAGACCACGGC<UCA>ACCGUGGGGGAGCGUGGGAUA}CG}CUCAGGC--{UAGACGGUGGGAGA{GGGUGGU-GGAAUUCCCGGAGUAGCGG<UGAA>AUGCGCAGAUACCGGGAGGAACGCCG-AUGGC<GAAG>GCA-GCCACCU}GGUCCACCCGUGAC}-GCUGAGGCGCG-AAAGCGUGGGGAGCAAACCGG<AUUAGAUAC>CCGGGUAGUCCACGCCCU--AAAC--{GAUGCGC--GCUAGGUCUCUGG<GUCUC>CUGGGGGCCGAAGC-UAA<CGCG>UUA----AGCGCGC}-CGCC}U--GGGGA{GUACGGCCG<CA>AGGCUGA}AACUCAAAGGAAU----{UG-ACGGGGGCCCGC{ACAA-GCGGUG-GAGCAUGUGGUU{UAAUU-CGA<AGCA>ACGCGAAGAACCUUAC-CAGGCCU-{UGACAUGCUAGGGAA-CCCGGGU<GAAA>GCCUGGGGUGC-CCC<GCGA>GGG-GAGCCCUAGCAC}--{AGGUGCU{GCAUG{GCCGUCGUCA-GCUCGU-GCC<GUGA>GGUGU-UGGGU<UAAG>UCCCG-CAACGAGCGCAAC-CCCCG-CCGUUAGUUGCCAG{CGG<UUCG>GCCG}GGCACUCUAACGG--{GACUGCCCGC-<GAA>-AGCGGGAGGAAGG}--AGGG-GACGACGUCUGGUCAG}CAU}GGC-CCUUAC}-GGCCUG--G}GCGACACACGUGCUA--C-AAUGCCCACUACAAAGCGAU{GCCACCCG<GCAA>CGGGGAGCU}AAUCGCAAAAAGGUGGGCCCAGUU-CGGAUUGGGGUCU<GCA>ACCCGACCCCAUGAAG-CCG-GAAUCGCU-A{GUAAUCGCGGAU<CAGC>CAUGCCGCGGUGA}A-UACGUUCC}CGGGCCUUGUACACA}--CCGCCCGUCACGCCAUGGGAGCGGGCUCUACCCGAAGUCGC{CGGGA---GCC<UA>CGG--GCAG}GCGCCGAGGGUAGGGCCCGUGACUGGGGCGAAGUCGUAACAAG-GU-AGCUGUACC<GGAA>GGUGCGGCU-GGAUCACUUUCU';

% ------------------------------------------- Determine alignment from JAR3D

[i1,i2] = pGetAlignedIndices(ES,TS);

zAlignmentDiagram(File,i1,i2);

a = 1;
b = 1;

j = 1;            % position in Java sequence
e = 1;            % position in experimental sequence

[Ecti,Eitc] = pColumnToIndex(ES);
[Tcti,Titc] = pColumnToIndex(TS);

Text{1} = '';
Text{2} = '';
Text{3} = '';
Text{4} = '';

for i = 1:length(Titc),
  Text{1} = [Text{1} ES(Eitc(i))];
  Text{2} = [Text{2} TS(Titc(i))];
end

fprintf('%s E. coli sequence\n', Text{1});
fprintf('%s T.th. sequence  \n', Text{2});

break

% ------------------------------------------------ Circular diagrams and counts

figure(7)
clf
zCircularDiagram(File,0.2,[1 1 1 1 0 0 0]);
saveas(gcf,'2avy_circular_basepairs.pdf','pdf');

figure(8)
clf
FF = File;
Node(1).Edge(File.NumNT,File.NumNT) = 0;
FF.Edge = sparse(Node(1).Edge + Node(1).Edge');
zCircularDiagram(FF,0.2,[1 1 1 1 0 0 0]);
saveas(gcf,'2avy_circular_SCFG_basepairs.pdf','pdf');



