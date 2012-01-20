
L = {'A','C','G','U'};

% display interacting oxygens and phosphorus atoms together with bases

figure(1)
clf
for v = 1:4,
  subplot(2,2,v);

  r = find(D(:,4) == v);               % select the current nucleotide code
  DD = D(r,:);                         % use only these lines of data

  s = 4*(DD(:,5) < 100) + 1*(DD(:,5) > 100);% true BP are large, near are small

  scatter3(DD(:,10),DD(:,11),DD(:,12), s,   DD(:,9), 'filled');
  hold on
%  scatter3(DD(:,13),DD(:,14),DD(:,15), s, 0*DD(:,9), 'filled');

  switch v,
    case 1,     text(8,0,'2BP');
                text(5,8,'6BP');
                text(-3,7.2,'7BP');
                text(-3.5,-3,'0BP');
    case 2,     text(5,6.1,'6BP');
                text(-2.3,8,'7BP');
                text(-4.5,6.3,'8BP');
                text(-5.8,4.7,'9BP');
                text(-4,-2.8,'0BP');
    case 3,     text(8,-2.4,'1BP');
                text(8.5,3,'3BP');
                text(7.3,4.8,'4BP');
                text(6,6,'5BP');
                text(-4.5,-2.9,'0BP');
    case 4,     text(6.2,3,'5BP');
                text(-3.9,6,'9BP');
                text(-2.9,-2,'0BP');
  end

  map = colormap;
  map(1,:) = [0 0 0];
  colormap(map);

  caxis([1 5]);

  zPlotStandardBase(v,1,0);                % plot base at the origin
  title(['Phosphate interactions with ' L{v}]);
  rotate3d on
  grid off
  axis([-6 9 -3.5 8.5]);
%  axis equal
  view(2)
%  saveas(gcf,['Phosphate Interactions\PhosphateInteractions_' L{v} '.fig'],'fig')

%  saveas(gcf,['Phosphate interactions' filesep 'Phosphate with ' L{v} '.png'],'png');

end

figure(2)
clf
hist(D(:,16),30);
max(D(:,16))

% display bond length and angle by nucleotide

figure(3)
clf
for v = 1:4,
  subplot(2,2,v)

  r = find(D(:,4) == v);               % select the current nucleotide code
  DD = D(r,:);                         % use only these lines of data

% s = 4*(DD(:,5) < 100) + 1*(DD(:,5) > 100);% true BP are large, near are small
  s = 4*(DD(:,17) == 1) + 1*(DD(:,17) == 0);% best oxygen is large

  scatter(DD(:,9),DD(:,8), s,   DD(:,6), 'filled');
  hold on
  plot([1 3.6 3.6], [150 150 180], 'k')
  plot([1 4.0 4.0], [120 120 180], 'k')
  title(['Oxygen location parameters for ', L{v}]);
  axis([2 5 100 180]);
end

% display bond length and angle by type of massive base

MT = {'Carbon','Self','Ring Nitrogen','Amino Nitrogen'};

ICode{1} = [1 8];                   % ring carbon, not C6/C8
ICode{2} = [4 9 14 17];             % carbon C6/C8, mostly self interactions 
ICode{3} = [2 3 5 6 10 11];         % amino nitrogen
ICode{4} = [13 15];                 % ring nitrogen

figure(4)
clf
for v = 1:4,
  subplot(2,2,v)

  r = [];
  for k = 1:length(ICode{v}),
    r = [r (find(mod(D(:,5),100) == ICode{v}(k)))'];  % append matches to this massive atom
  end

  DD = D(r,:);                         % use only these lines of data

%  s = 4*(DD(:,5) < 100) + 1*(DD(:,5) > 100);% true BP are large, near are small
  s = 4*(DD(:,17) == 1) + 1*(DD(:,17) == 0);% best oxygen is large

  scatter(DD(:,9),DD(:,8), s,   DD(:,7), 'filled');
  hold on
  plot([1 3.6 3.6], [150 150 180], 'k')
  plot([1 4.0 4.0], [120 120 180], 'k')
  title(['Oxygen location parameters for ', MT{v}]);
  axis([2 5 100 180]);
end

