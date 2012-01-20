% zPhosphateInteractions checks all nearby pairs of bases for base-phosphate
% interactions, and stores them in a sparse matrix field BasePhosphate

function [File] = zPhosphateInteractions(File,Verbose)

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

File(f).BasePhosphate = sparse(zeros(File(f).NumNT));

% -------- First screening of base pairs ------------------------------------ 

DistCutoff = 16;                                % max distance for interaction
[i,j] = find((File(f).Distance < DistCutoff).*(File(f).Distance > 0)); 
                                                % screen by C-C distance

% -------- Screen and analyze base pairs ------------------------------------ 
% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

pc = 1;                                         % index for pairs

for k = 1:length(i),                            % loop through possible pairs

  N1 = File(f).NT(i(k));                        % nucleotide i information
  N2 = File(f).NT(j(k));                        % nucleotide j information

  ph = (N2.Sugar(10,:)-N1.Center) * N1.Rot;     % phosphorus displacement  

  c = N1.Code;

  % preliminary elliptical screening based on the location of the phosphorus

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

    p   = [9 11 12 13];               % rows of the phosphate oxygens
    dis = zDistance(N1.Fit(h,:), N2.Sugar(p,:)); % distances between hyd & O's

    [ii,jj] = find(dis < 3.2);        % only keep those within 3.2 Angstroms

    g = 0;

    for kk = 1:length(ii),
      Angle=zAngle(N1.Fit(m(ii(kk)),:),N1.Fit(h(ii(kk)),:),N2.Sugar(p(jj(kk)),:));
      if Angle > 120,

        if ((Angle < 150) || (dis(ii(kk),jj(kk)) > 2.6)) %near case
          if (g == 0)
            g = e(ii(kk)) + 100;
          end
        elseif (g < 10) && (g > 0),
          g = g + e(ii(kk))/10;
        else
          g = e(ii(kk));
        end

        File(f).BasePhosphate(i(k),j(k)) =  g;

        if Verbose > 1,
%[kk ii(kk) jj(kk) f i(k) j(k) p(jj(kk))]
%size(jj)
%size(ii)
%h(ii(kk))
%p(jj(kk))
%dis(ii(kk),jj(kk))
%          fprintf('%6s base %s%5s %3s phosphate %s%5s %3s length %6.2f angle %6.2f interaction %s\n', File(f).Filename, N1.Base, N1.Number, AtomNames{h(ii(kk)),c}, N2.Base, N2.Number, Sugar{p(jj(kk))}, dis(ii(kk),jj(kk)), Angle, zEdgeText(File(f).Edge(i(k),j(k))));

          ox = (N2.Sugar(p(jj(kk)),:)-N1.Center) * N1.Rot; % oxygen displ

          a = [ox dis(ii(kk),jj(kk)) Angle File(f).Distance(i(k),j(k)) ph];
          count(c) = count(c) + 1;
          switch c
            case 1,     PHA(count(c),:) = a;
            case 2,     PHC(count(c),:) = a;
            case 3,     PHG(count(c),:) = a;
            case 4,     PHU(count(c),:) = a;
          end

        end

      end
    end

    if Verbose > 1,
      fprintf('\n');
    end
  end
  
end   % loop over pairs
end   % loop over files



if Verbose > 0,

for v = 1:4,
  figure(v)
  clf
  switch v,
    case 1,     c = PHA(:,4);
                scatter3(PHA(:,1), PHA(:,2), PHA(:,3),4,c,'filled')
                hold on
                scatter3(PHA(:,7), PHA(:,8), PHA(:,9),4,0*c,'filled')
    case 2,     c = PHC(:,4);
                scatter3(PHC(:,1), PHC(:,2), PHC(:,3),4,c,'filled')
                hold on
                scatter3(PHC(:,7), PHC(:,8), PHC(:,9),4,0*c,'filled')
    case 3,     c = PHG(:,4);
                scatter3(PHG(:,1), PHG(:,2), PHG(:,3),4,c,'filled')
                hold on
                scatter3(PHG(:,7), PHG(:,8), PHG(:,9),4,0*c,'filled')
    case 4,     c = PHU(:,4);
                scatter3(PHU(:,1), PHU(:,2), PHU(:,3),4,c,'filled')
                hold on
                scatter3(PHU(:,7), PHU(:,8), PHU(:,9),4,0*c,'filled')
  end

  L = {'A','C','G','U'};

  map = colormap;
  map(1,:) = [0 0 0];
  colormap(map);

  caxis([1 4]);

  zPlotStandardBase(v,1,1);                % plot base at the origin
  rotate3d on
  grid off
  axis equal
  view(2)
  saveas(gcf,['Phosphate Interactions\CodedPhosphateInteractions_' L{v} '.fig'],'fig')
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
hist(nonzeros(File(1).BasePhosphate),30)
pause
hist(nonzeros(File(2).BasePhosphate),30)
pause
hist(nonzeros(File(3).BasePhosphate),30)
pause
hist(nonzeros(File(4).BasePhosphate),30)
