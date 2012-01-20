% zAlignmentDiagram takes an alignment of the NTs in File(1) and File(2) and produces circular diagrams showing which NTs are aligned, which basepairs are aligned, etc.

% File = zAddNTData({'2avy','1j5e'});

% function [void] = zAlignmentDiagram(File,Aligned1,Aligned2)

Aligned1 = [3 4 5 7 8 10];
Aligned2 = [3 4 5 9 10 15];

Aligned1 = [10:500 length(File(1).NT)];
Aligned2 = [10:500 length(File(2).NT)];

Aligned1 = [ind2' length(File(1).NT)];
Aligned2 = [ind1' length(File(2).NT)];

% R1 is a function from the indices of File(1) to the row of the alignment

R1 = 1:(Aligned1(1)-1);
R2 = R1(end) + (1:(Aligned2(1)-1));

for a = 1:(length(Aligned1)-1)
  r  = R2(end)+1;
  R1 = [R1 r];
  R2 = [R2 r];

  R1 = [R1 r+(1:(Aligned1(a+1)-1-Aligned1(a)))];
  R2 = [R2 R1(end)+(1:(Aligned2(a+1)-1-Aligned2(a)))];
end

r  = R2(end)+1;
R1 = [R1 r];
R2 = [R2 r];

%%%%%%%%%%% need to add in any NTs after the last aligned NTs!!!!!

% R1(Aligned1) should equal R2(Aligned2)

clear NewFile

m = max(max(R1),max(R2));

N1 = File(1).NT(1);
N1.Number = '';
N1.Base = '';

for f = 1:2,
  NewFile(f).Filename      = File(f).Filename;
  NewFile(f).Edge          = sparse(zeros(m,m));
  NewFile(f).Crossing      = sparse(zeros(m,m));
  NewFile(f).BasePhosphate = sparse(zeros(m,m));
  NewFile(f).Covalent      = sparse(zeros(m,m));
  for n = 1:m,
    NewFile(f).NT(n) = N1;
  end
end


for r = 1:length(R1),
  NewFile(1).NT(R1(r)) = File(1).NT(r);
end

NewFile(1).Edge(R1,R1) = File(1).Edge;
NewFile(1).Crossing(R1,R1) = File(1).Crossing;
NewFile(1).BasePhosphate(R1,R1) = File(1).BasePhosphate;
NewFile(1).Covalent(R1,R1) = File(1).Covalent;

for r = 1:length(R2),
  NewFile(2).NT(R2(r)) = File(2).NT(r);
end

NewFile(2).Edge(R2,R2) = File(2).Edge;
NewFile(2).Crossing(R2,R2) = File(2).Crossing;
NewFile(2).BasePhosphate(R2,R2) = File(2).BasePhosphate;
NewFile(2).Covalent(R2,R2) = File(2).Covalent;

figure(1)
clf
zCircularDiagram(NewFile(1),0.5,[1 1 1 1 1 1 1 0]);
saveas(gcf,[NewFile(1).Filename '_alignment.pdf'],'pdf');

figure(2)
clf
zCircularDiagram(NewFile(2),0.5,[1 1 1 1 1 1 1 0]);
saveas(gcf,[NewFile(2).Filename '_alignment.pdf'],'pdf');

figure(3)
clf
Comp = NewFile(1);
Comp.Filename = [Comp.Filename ' conserved'];
Comp.Edge  = Comp.Edge .* (fix(abs(NewFile(1).Edge)) == fix(abs(NewFile(2).Edge)));
Comp.BasePhosphate  = Comp.BasePhosphate .* (fix(abs(NewFile(1).BasePhosphate/100)) == fix(abs(NewFile(2).BasePhosphate/100)));
zCircularDiagram(Comp,0.5,[1 1 1 1 1 1 1 0]);
saveas(gcf,[NewFile(1).Filename '-' NewFile(2).Filename '_alignment_conserved.pdf'],'pdf');
