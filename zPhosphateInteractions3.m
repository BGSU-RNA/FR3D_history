
%function [void] = zPhosphateInteractions(File)

PHA = [];
PHC = [];
PHG = [];
PHU = [];
count = [1 1 1 1];

Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

for f = 1:length(File),

% -------- First screening of base pairs ------------------------------------ 

DistCutoff = 15;                                % max distance for interaction
[i,j] = find((File(f).Distance < DistCutoff).*(File(f).Distance > 0)); 
                                                % screen by C-C distance
k = find(abs(i-j) > 1);
i = i(k);                                       % exclude sequential pairs
j = j(k);

k = find(j > 1);
i = i(k);
j = j(k);

% -------- Screen and analyze base pairs ------------------------------------ 
% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

pc = 1;                                         % index for pairs

for k = 1:length(i),                            % loop through possible pairs

  N1 = File(f).NT(i(k));                           % nucleotide i information
  N2 = File(f).NT(j(k));                           % nucleotide j information
  N3 = File(f).NT(max(1,j(k)-1));

  ph = (N2.Sugar(10,:)-N1.Center) * N1.Rot; % phosphate displacement  
  P = 'O1P'; S = N2.Sugar(11,:);
  o1 = (S - N1.Center) * N1.Rot;
  P = 'O2P'; S = N2.Sugar(12,:);
  o2 = (S - N1.Center) * N1.Rot;
  P = 'O5'; S = N2.Sugar( 9,:);
  o5 = (S - N1.Center) * N1.Rot;
  P = 'O3'; S = N3.Sugar( 5,:);
  o3 = (S - N1.Center) * N1.Rot;

  P = 'All';

  c = N1.Code;

  if (ph*diag([1 1 6])*ph' < 60) && (abs(ph(3)) < 2) && norm(o3) < 10,
    col = min(zDistance(N1.Fit(1:Lim(2,N1.Code),:), S)); 
    a = [[ph 0]; [o1 1]; [o2 2]; [o5 3]; [o3 4]];
    switch c
      case 1,     PHA(count(c):count(c)+4,:) = a;
      case 2,     PHC(count(c):count(c)+4,:) = a;
      case 3,     PHG(count(c):count(c)+4,:) = a;
      case 4,     PHU(count(c):count(c)+4,:) = a;
    end
    count(c) = count(c) + 5;
  end

end   % loop over pairs
end   % loop over files

for v = 1:4,
  figure(v)
  clf
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

  map = colormap;
  map(1,:) = [0 0 1];
  colormap(map);

  caxis([0 4.2]);

  L = {'A','C','G','U'};

  zPlotStandardBase(v,1,1);                % plot base at the origin
  rotate3d on
  axis equal
  view(2)
  saveas(gcf,['PhosphateInteractions_' P '_' L{v} '.pdf'],'pdf')
  saveas(gcf,['PhosphateInteractions_' P '_' L{v} '.fig'],'fig')

end
