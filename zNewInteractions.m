% zNewInteractions checks all nearby pairs of bases for base-phosphate
% interactions, and stores them in a sparse matrix field BasePhosphate

function [File] = zNewInteractions(File,Verbose,Atom)

if nargin == 1,
  Verbose = 0;
end

PHA = [];
PHC = [];
PHG = [];
PHU = [];
count = [0 0 0 0];

zStandardBases
Sugar = {'C1*','C2*','O2*','C3*','O3*','C4*','O4*','C5*','O5*','P','O1P','O2P'};

Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

for f = 1:length(File),

%File(f).BaseOther = sparse(zeros(File(f).NumNT));

% -------- First screening of base pairs ------------------------------------ 

DistCutoff = 16;                                % max distance for interaction

S = (File(f).Distance < DistCutoff) .* (File(f).Distance > 0) .* ((File(f).Edge == 0) + (fix(File(f).Edge) == 30)) .* (File(f).BasePhosphate == 0);

%S = (File(f).Distance < DistCutoff) .* (File(f).Distance > 0);

%S = (File(f).Distance < DistCutoff) .* (File(f).Distance > 0) .* ((File(f).Edge == 0) + (fix(File(f).Edge) == 30));

[i,j] = find(S);                                % screen by C-C distance

% -------- Screen and analyze base pairs ------------------------------------ 
% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

pc = 1;                                         % index for pairs

PHA = zeros(1,6);
PHC = PHA;
PHG = PHA;
PHU = PHA;

for k = 1:length(i),                            % loop through possible pairs

for Atom = 1:12,

  N1 = File(f).NT(i(k));                        % nucleotide i information
  N2 = File(f).NT(j(k));                        % nucleotide j information

  L = Lim(2,N1.Code);                           % number of atoms in N1

  if Atom <= 12,
    at = N2.Sugar(Atom,:);
    ph = (N2.Sugar(Atom,:)-N1.Center) * N1.Rot;     % atom displacement  
  else
    at = N2.Fit(Atom-12,:);
    ph = (N2.Fit(Atom-12,:)-N1.Center) * N1.Rot;
  end

  c = N1.Code;

  % preliminary elliptical screening based on the location of the atom

  if (ph*diag([1 1 6])*ph' < 60) && (abs(ph(3)) < 2),
    switch c
      case 1,                         % Base A
              h   = [11 12 14 15];    % rows of the base hydrogens
              m   = [ 9  7  6  6];    % rows of the corresponding massive atoms
              e   = [ 1  4  2  3];    % code for location of the interaction
      case 2, 
              h   = [10 11 12 13];
              m   = [ 7  8  6  6];
              e   = [ 4  3  3  2];
      case 3, 
              h   = [12 13 15 16];
              m   = [ 4  7 11 11];
              e   = [ 2  4  1  2];
      case 4, 
              h   = [ 9 11 12];
              m   = [ 8  4  7];
              e   = [ 3  2  4];
    end

    dis = zDistance(N1.Fit(1:L,:), at); % distances between base and new atom
    [ii,jj] = find(dis < 5);        % only keep those close enough

    g = 0;

    for kk = 1:length(ii),
      if Verbose > 1,
%          fprintf('%6s base %s%5s %3s phosphate %s%5s %3s length %6.2f angle %6.2f interaction %s\n', File(f).Filename, N1.Base, N1.Number, AtomNames{h(ii(kk)),c}, N2.Base, N2.Number, Sugar{p(jj(kk))}, dis(ii(kk),jj(kk)), Angle, zEdgeText(File(f).Edge(i(k),j(k))));

          a = [ph Atom dis(ii(kk),jj(kk)) File(f).Distance(i(k),j(k))];
          count(c) = count(c) + 1;
          switch c
            case 1,     PHA(count(c),:) = a;
            case 2,     PHC(count(c),:) = a;
            case 3,     PHG(count(c),:) = a;
            case 4,     PHU(count(c),:) = a;
          end

      end
    end

    if Verbose > 1,
%      fprintf('\n');
    end
  end
  
end   % loop over pairs

end   % loop over Atoms

end   % loop over files



if Verbose > 0,

for v = 1:4,
  figure(v)
  clf
  hold on
  switch v,
    case 1,     c = PHA(:,4);
                scatter3(PHA(:,1), PHA(:,2), PHA(:,3),4,c,'filled')
    case 2,     c = PHC(:,4);
                scatter3(PHC(:,1), PHC(:,2), PHC(:,3),4,c,'filled')
    case 3,     c = PHG(:,4);
                scatter3(PHG(:,1), PHG(:,2), PHG(:,3),4,c,'filled')
    case 4,     c = PHU(:,4);
                scatter3(PHU(:,1), PHU(:,2), PHU(:,3),4,c,'filled')
  end

  L = {'A','C','G','U'};

%  map = colormap;
%  map(1,:) = [0 0 0];
%  colormap(map);

  caxis([1 12]);

  zPlotStandardBase(v,1,1);                % plot base at the origin
  rotate3d on
  grid off
  axis equal
  view(2)
  title('Nearby sugar atoms, blue=C1,C2,O2,C3,O3,C4,O4,C5,O5,P,O1P,O2P=red');
  saveas(gcf,['Phosphate Interactions\BaseSugarInteractions_KnownReactionsRemoved' L{v} '.fig'],'fig')
end

end

return

for v = 1:4,
  figure(v+4)
  clf
  switch v,
    case 1,     hist(PHA(:,4),30);
    case 2,     hist(PHC(:,4),30);
    case 3,     hist(PHG(:,4),30);
    case 4,     hist(PHU(:,4),30);
  end
end

for v = 1:4,
  figure(v+8)
  clf
  switch v,
    case 1,     hist(PHA(:,5),30);
    case 2,     hist(PHC(:,5),30);
    case 3,     hist(PHG(:,5),30);
    case 4,     hist(PHU(:,5),30);
  end
end

for v = 1:4,
  figure(v+12)
  clf
  switch v,
    case 1,     plot(PHA(:,4),PHA(:,5),'.');
    case 2,     plot(PHC(:,4),PHC(:,5),'.');
    case 3,     plot(PHG(:,4),PHG(:,5),'.');
    case 4,     plot(PHU(:,4),PHU(:,5),'.');
  end
end


clf
hist(nonzeros(File(1).BaseOther),30)
pause
hist(nonzeros(File(2).BaseOther),30)
pause
hist(nonzeros(File(3).BaseOther),30)
pause
hist(nonzeros(File(4).BaseOther),30)
