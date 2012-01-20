% zPhosphateInteractions checks all nearby pairs of bases for base-phosphate
% interactions, and stores them in a sparse matrix field BasePhosphate

function [File,PH] = zPhosphateInteractions(File,Verbose)

if nargin == 1,
  Verbose = 0;
end

D = [];                   % where to save data if Verbose
count = [0 0 0 0];

zStandardBases
Sugar = {'C1*','C2*','O2*','C3*','O3*','C4*','O4*','C5*','O5*','P','O1P','O2P','O3 of next'};

Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

t = cputime;

for f = 1:length(File),

  if isempty(File(f).Distance),
    c = cat(1,File(f).NT(1:File(f).NumNT).Center); % nucleotide centers
    File(f).Distance = zMutualDistance(c,16); % compute distances < 16 Angs
  end

  File(f).BasePhosphate = sparse(zeros(File(f).NumNT));

  % -------- First screening of base pairs ----------------------------------- 

  DistCutoff = 16;                              % max distance for interaction
  [i,j] = find((File(f).Distance < DistCutoff).*(File(f).Distance > 0)); 
                                                % screen by C-C distance

  i = [i; (1:length(File(f).NT))'];             % allow self interactions
  j = [j; (1:length(File(f).NT))'];             % allow self interactions

  % -------- Screen and analyze pairs ----------------------------------------

  pc = 1;                                       % counter for valid pairs
  p   = [9 11 12 13];                           % rows of the phosphate oxygens
  pn  = {'O5*','O1P','O2P','O3*'};              % names of phosphate oxygens

  for k = 1:length(i),                          % loop through possible pairs

    N1 = File(f).NT(i(k));                      % nucleotide i information
    N2 = File(f).NT(j(k));                      % nucleotide j information

    switch N1.Code
      case 1,                         % Base A
              h   = [11 12 14 15];    % rows of the base hydrogens
              hn  = {'H2','H8','1H6','2H6'}; % names of the base hydrogens
              m   = [ 9  7  6  6];    % rows of the corresponding massive atoms
              e   = [ 1  4  2  3];    % code for location of the interaction
      case 2,                         % Base C
              h   = [10 11 12 13];
              hn  = {'H6','H5','1H4','2H4'}; % names of the base hydrogens
              m   = [ 7  8  6  6];
              e   = [ 9  8  6  5];
      case 3,                         % Base G
              h   = [12 13 15 16];
              hn  = {'H1','H8','1H2','2H2'}; % names of the base hydrogens
              m   = [ 4  7 11 11];
              e   = [13 14 10 11];
      case 4,                         % Base U
              h   = [ 9 11 12];
              hn  = {'H5','H3','H6'}; % names of the base hydrogens
              m   = [ 8  4  7];
              e   = [16 15 17];
    end

    dis = zDistance(N1.Fit(m,:), N2.Sugar(p,:)); % distances between mass & O's
    [mm,pp] = find(dis < 4.5);        % massive-oxygen pairs close enough

    g = [];                           % internal classification number

    for n = 1:length(mm),             % loop through potential matches
      Angle=zAngle(N1.Fit(m(mm(n)),:),N1.Fit(h(mm(n)),:),N2.Sugar(p(pp(n)),:));
                                      % base massive - hydrogen - oxygen angle
      if Angle > 100,                 % good enough to be "near" base-phosph

        if ((Angle < 150) || (dis(mm(n),pp(n)) > 3.6)) % something imperfect
          if isempty(g)               % no classification yet
            g = e(mm(n)) + 100;       % > 100 means "near"
          end
        else
          g = [g e(mm(n))];           % assign a first, non-near class.
        end

        % store information for later display

        if Verbose > 1,
          ox = (N2.Sugar(p(pp(n)),:)-N1.Fit(1,:)) * N1.Rot; % oxygen displ
          ph2= (N2.Sugar(10,:)-N1.Fit(1,:)) * N1.Rot; % phosphorus displacement

          a = [f i(k) j(k) N1.Code g(end) mm(n) pp(n) Angle dis(mm(n),pp(n)) ox ph2 File(f).Distance(i(k),j(k))];

          D = [D; a];                  % append data to data matrix

          if Verbose > 3,
            fprintf('%6s base %s%5s %3s %3d phosphate %s%5s %13s length %6.2f angle %6.2f interaction %s\n', File(f).Filename, N1.Base, N1.Number, AtomNames{h(mm(n)),N1.Code}, g(end), N2.Base, N2.Number, Sugar{p(pp(n))}, dis(mm(n),pp(n)), Angle, zEdgeText(File(f).Edge(i(k),j(k))));
          end
        end

      end

      if length(g) > 0,
        if (min(g) < 100) && (max(g) > 100),
          g = g(find(g < 100));
        end
        if length(g) > 1,
          g = sort(g);
          if (g(1) == 6) && (g(end) == 8),
            g = 7;
          elseif (g(1) == 11) && (g(end) == 13),
            g = 12;
          else
            fprintf('Another case to consider\n');
            g
            [f i(k) j(k)]
          end
        end

        File(f).BasePhosphate(i(k),j(k)) =   g(1);  % record classification
      end

    end   % loop over massive atom and oxygen pairs
  end     % loop over nucleotide pairs
end       % loop over files

if Verbose > 1,

fprintf('Classifying base-phosphate interactions took %8.2f minutes\n', (cputime-t)/60);

for v = 1:4,
  figure(v)
  clf
  switch v,
    case 1,     c = PHA(:,4);
                scatter3(PHA(:,7), PHA(:,8), PHA(:,9),4,0*c,'filled')
                hold on
                scatter3(PHA(:,1), PHA(:,2), PHA(:,3),4,c,'filled')
                text(8,0,'2BP');
                text(5,8,'6BP');
                text(-3,7.2,'7BP');
                text(-3.5,-3,'0BP');
    case 2,     c = PHC(:,4);
                scatter3(PHC(:,7), PHC(:,8), PHC(:,9),4,0*c,'filled')
                hold on
                scatter3(PHC(:,1), PHC(:,2), PHC(:,3),4,c,'filled')
                text(5,6.1,'6BP');
                text(-2.3,8,'7BP');
                text(-4.5,6.3,'8BP');
                text(-5.8,4.7,'9BP');
                text(-4,-2.8,'0BP');
    case 3,     c = PHG(:,4);
                scatter3(PHG(:,7), PHG(:,8), PHG(:,9),4,0*c,'filled')
                hold on
                scatter3(PHG(:,1), PHG(:,2), PHG(:,3),4,c,'filled')
                text(8,-2.4,'1BP');
                text(8.5,3,'3BP');
                text(7.3,4.8,'4BP');
                text(6,6,'5BP');
                text(-4.5,-2.9,'0BP');
    case 4,     c = PHU(:,4);
                scatter3(PHU(:,7), PHU(:,8), PHU(:,9),4,0*c,'filled')
                hold on
                scatter3(PHU(:,1), PHU(:,2), PHU(:,3),4,c,'filled')
                text(6.2,3,'5BP');
                text(-3.9,6,'8BP');
                text(-2.9,-2,'0BP');
  end

  L = {'A','C','G','U'};

  map = colormap;
  map(1,:) = [0 0 0];
  colormap(map);

  caxis([1 4]);

  zPlotStandardBase(v,1,0);                % plot base at the origin
  title(['Phosphate interactions with ' L{v}]);
  rotate3d on
  grid off
  axis([-6 9 -3.5 8.5]);
%  axis equal
  view(2)
  saveas(gcf,['Phosphate Interactions\PhosphateInteractions_' L{v} '.fig'],'fig')

  saveas(gcf,['Phosphate interactions' filesep 'Phosphate with ' L{v} '.png'],'png');

end

figure(5)
clf
hist(sdist,30);
max(sdist)

figure(6)
clf

subplot(2,2,1)
a = mod(PHA(:,10),100);
b = sort(unique(a));
c(b) = 1:length(b);
a = c(a);
scatter(PHA(:,4),PHA(:,5),14*ones(size(a)),a,'Filled');
hold on
plot([1 2.6 2.6], [150 150 180], 'k')
plot([1 3.2 3.2], [120 120 180], 'k')
title('Oxygen location parameters for A');
axis([1 4 100 180]);

subplot(2,2,2)
a = mod(PHC(:,10),100);
b = sort(unique(a));
c(b) = 1:length(b);
a = c(a);
scatter(PHC(:,4),PHC(:,5),14*ones(size(a)),a,'Filled');
hold on
plot([1 2.6 2.6], [150 150 180], 'k')
plot([1 3.2 3.2], [120 120 180], 'k')
title('Oxygen location parameters for C');
axis([1 4 100 180]);

subplot(2,2,3)
a = mod(PHG(:,10),100);
b = sort(unique(a));
c(b) = 1:length(b);
a = c(a);
scatter(PHG(:,4),PHG(:,5),14*ones(size(a)),a,'Filled');
hold on
plot([1 2.6 2.6], [150 150 180], 'k')
plot([1 3.2 3.2], [120 120 180], 'k')
title('Oxygen location parameters for G');
axis([1 4 100 180]);

subplot(2,2,4)
a = mod(PHU(:,10),100);
b = sort(unique(a));
c(b) = 1:length(b);
a = c(a);
scatter(PHU(:,4),PHU(:,5),14*ones(size(a)),a,'Filled');
hold on
plot([1 2.6 2.6], [150 150 180], 'k')
plot([1 3.2 3.2], [120 120 180], 'k')
title('Oxygen location parameters for U');
axis([1 4 100 180]);

end




if Verbose > 2,

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

end

return

File = zAddNTData({'1s72','1j5e','2avy','2aw4','2j01'});
zPhosphateInteractions(File,3);
