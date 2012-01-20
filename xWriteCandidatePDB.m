
function [void] = xWriteCandidatePDB(File,Search)

N = Search.Query.NumNT;                        % number of nucleotides in each

if isfield(Search.Query,'LocWeight'),
  LW = Search.Query.LocWeight;
else
  LW = ones(N,1);
end

if isfield(Search.Query,'AngleWeight'),
  AW = Search.Query.AngleWeight;
else
  AW = ones(N,1);
end


A = {'N9' 'C4' 'N3' 'N1' 'C6' 'N6' 'C8' 'C5' 'C2' 'N7'};
C = {'N1' 'C2' 'O2' 'N3' 'C4' 'N4' 'C6' 'C5'};
G = {'N9' 'C4' 'N3' 'N1' 'C6' 'O6' 'C8' 'C5' 'C2' 'N7' 'N2'};
U = {'N1' 'C2' 'O2' 'N3' 'C4' 'O4' 'C6' 'C5'};
S = {'C1*' 'C2*' 'O2*' 'C3*' 'O3*' 'C4*' 'O4*' 'C5*' 'O5*' 'P' 'O1P' 'O2P'};

Lim(1,:) = [10 8 11 8];       % number of base atoms, excluding hydrogen
Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

if length(Search.Candidates) > 0,

fid = fopen([Search.SaveName '-Cand.pdb'],'w');       % open for writing

M = length(Search.Candidates(:,1));            % number of candidates

f     = Search.Candidates(1,N+1);
Model = File(f).NT(Search.Candidates(1,1:N));  % first cand, taken as model


a = 1;                                         % atom number

for c = 1:M,                                   % loop through candidates
 f     = Search.Candidates(c,N+1);             % file number, this candidate
 Cand  = File(f).NT(Search.Candidates(c,1:N)); % current candidate
 [R,Sh] = xSuperimposeCandidates(Model,Cand,LW,AW);

 for i = 1:N,                                  % loop through nucleotides
  NT = Cand(i);                                % current nucleotide
  for j = 1:Lim(1,NT.Code),                    % loop through base atoms
    fprintf(fid, 'ATOM  %5d', a);
    switch NT.Code,
      case 1, fprintf(fid,'  %-3s', A{j});
      case 2, fprintf(fid,'  %-3s', C{j});
      case 3, fprintf(fid,'  %-3s', G{j});
      case 4, fprintf(fid,'  %-3s', U{j});
    end
    fprintf(fid, '  %1s', NT.Base);
    fprintf(fid, ' %1s',  NT.Chain);
    fprintf(fid, '%4s',   NT.Number);
    L = (NT.Fit(j,:) - Sh)*R' + (c-1)*[20 0 0];
    fprintf(fid, '%8.3f', L);                  % write atom location
    fprintf(fid, '%6.2f', 0);
    fprintf(fid, '%6.2f\n', 0);
    a = a + 1;
  end
  for j = 1:length(NT.Sugar(:,1)),             % loop through sugar atoms
    fprintf(fid, 'ATOM  %5d', a);
    fprintf(fid,'  %-3s', S{j});
    fprintf(fid, '  %1s', NT.Base);
    fprintf(fid, ' %1s',  NT.Chain);
    fprintf(fid, '%4s',   NT.Number);
    L = (NT.Sugar(j,:) - Sh)*R' + (c-1)*[20 0 0];
    fprintf(fid, '%8.3f', L);
    fprintf(fid, '%6.2f', 0);
    fprintf(fid, '%6.2f\n', 0);
    a = a + 1;
  end
 end
end

fclose(fid);

fprintf('Wrote %s\n', [Search.SaveName '-Cand.pdb']);

end