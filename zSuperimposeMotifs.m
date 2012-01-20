
Motif1.Filename = 'rr0033_23S';
Motif1.Bases = {'469' '470' '471' '472'}';
Motif1.AllBases    = {'469' '470' '471' '472'}';

Motif2.Filename = 'rr0033_23S';
Motif2.Bases = {'1629' '1630' '1631' '1632'};
Motif2.AllBases    = {'1629' '1630' '1631' '1632'};

Motif1 = zGetMotifData(Motif1);
Motif2 = zGetMotifData(Motif2);

% ----------------- Mike superimposes the two motifs, calls them Rotate1, 2

Rotate1 = Motif1; % temporary
Rotate2 = sRotateMotif(Motif1, Motif2);

% ----------------- Write the new pdb files

zWriteMotifPDB(Rotate1);
zWriteMotifPDB(Rotate2);

