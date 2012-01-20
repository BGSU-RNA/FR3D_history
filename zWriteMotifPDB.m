
function [void] = zWriteMotifPDB(M)

  A = {'N9' 'C4' 'N3' 'N1' 'C6' 'N6' 'C8' 'C5' 'C2' 'N7'};
  C = {'N1' 'C2' 'O2' 'N3' 'C4' 'N4' 'C6' 'C5'};
  G = {'N9' 'C4' 'N3' 'N1' 'C6' 'O6' 'C8' 'C5' 'C2' 'N7' 'N2'};
  U = {'N1' 'C2' 'O2' 'N3' 'C4' 'O4' 'C6' 'C5'};
  S = {'C1*' 'C2*' 'O2*' 'C3*' 'O3*' 'C4*' 'O4*' 'C5*' 'O5*' 'P' 'O1P' 'O2P'};

Lim(1,:) = [10 8 11 8];       % number of base atoms, excluding hydrogen
Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

fid = fopen([M.Filename '-' M.Bases{1} '.pdb'],'w');       % open for writing
a = 1;                                         % atom number
for i = 1:length(M.AllNT),
  for j = 1:Lim(1,M.AllNT(i).Code),
    fprintf(fid, 'ATOM  %5d', a);
    switch M.AllNT(i).Code,
      case 1, fprintf(fid,'  %-3s', A{j});
      case 2, fprintf(fid,'  %-3s', C{j});
      case 3, fprintf(fid,'  %-3s', G{j});
      case 4, fprintf(fid,'  %-3s', U{j});
    end
    fprintf(fid, '  %1s', M.AllNT(i).Base);
    fprintf(fid, ' %1s',  M.AllNT(i).Chain);
    fprintf(fid, '%4s',   M.AllNT(i).Number);
    fprintf(fid, '%8.3f', M.AllNT(i).Loc(j,:));
    fprintf(fid, '%6.2f', 0);
    fprintf(fid, '%6.2f\n', 0);
    a = a + 1;
  end
  for j = 1:length(M.AllNT(i).Sugar(:,1)),
    fprintf(fid, 'ATOM  %5d', a);
    fprintf(fid,'  %-3s', S{j});
    fprintf(fid, '  %1s', M.AllNT(i).Base);
    fprintf(fid, ' %1s',  M.AllNT(i).Chain);
    fprintf(fid, '%4s',   M.AllNT(i).Number);
    fprintf(fid, '%8.3f', M.AllNT(i).Sugar(j,:));
    fprintf(fid, '%6.2f', 0);
    fprintf(fid, '%6.2f\n', 0);
    a = a + 1;
  end
end

fclose(fid);

fprintf('Wrote %s\n', [M.Filename '-' M.Bases{1} '.pdb']);