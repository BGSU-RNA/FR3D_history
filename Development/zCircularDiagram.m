% zCircularDiagram(File,Thickness) plots the pairwise interactions in File
% using colored chords around a unit circle.  Thickness controls the
% thickness of the lines, for different graphic output formats.

% zCircularDiagram('1s72') will load the file and display the diagram.

% Helpful suggestions for how to save the figure as a graphic file:
%  clf
%  zCircularDiagram(File(f),1);
%  saveas(gcf,[mypath FN '_circular_diagram.png'],'png');
%  [X,map] = imread([mypath FN '_circular_diagram.png']);
%  Y = X(30:830,210:1030,:);
%  imwrite(Y,[mypath FN '_circular_diagram.png']);

%  clf
%  zCircularDiagram(File(f),0.1);
%  saveas(gcf,[mypath FN '_circular_diagram.pdf'],'pdf');

%  zCircularDiagram(File,1,[1 1 0 0 0 1 1]);

%  View tells what to display
%  1 = nested cWW
%  2 = nested non-cWW
%  3 = non-nested cWW
%  4 = non-nested non-cWW
%  5 = stack
%  6 = BPh
%  7 = explanatory text
%  8 = leave gaps in diagram for breaks in nucleotide chain

function [void] = zCircularDiagram(File,Thickness,View)

if nargin < 2,
  Thickness = 1;
end

if nargin < 3,
  View = [1 1 1 1 1 1 1 1];
end

while length(View) < 8,
  View = [View 1];
end

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

E  = fix(abs(File.Edge));
B  = E .* (E > 0) .* (E < 24);                 % pairs and stacks
%R  = File.Range;
C  = File.Crossing;

Color = (B==1).*(C==0) + 2*(B>1).*(B<14).*(C==0) + 3*(B==1).*(C>0) + 4*(B > 1).*(B < 14) .*(C>0) + 5*(B > 20) .* (B < 25);
                                        % disjoint possibilities

BP = abs(File.BasePhosphate);           % 
BP = (BP > 0) .* (BP < 100);            % exclude near BP and self interactions

[i,j,c] = find(triu(Color));

[ii,jj,cc] = find(6*BP);                % handle this separately, don't add

%[ii,jj,cc] = find(6*BP.*(R > 1));     % range 0 BPh only

k = find(ii ~= jj);                     % eliminate self interactions

i = [i; ii(k)];
j = [j; jj(k)];
c = [c; cc(k)];

if length(i) > 0,

N = length(File.NT);

% ------------------------------------- Determine where to plot each NT

A = zeros(1,N);           % angle around circle for nucleotides
A(1) = 1;

if View(8) == 1,
  spaces = [1 4 8 12];
else
  spaces = [1 0 0 12];
end

nc = [1];                                % new chain starting
for t = 1:(N-1),
  A(t+1) = A(t) + spaces(1);
  if File.Covalent(t,t+1) == 0,
    A(t+1) = A(t+1) + spaces(2);
  end
  if (File.NT(t).Chain ~= File.NT(t+1).Chain),
    A(t+1) = A(t+1) + spaces(3);
    nc = [nc t+1];
  end
end
A(end+1) = A(end) + spaces(4);
mA = max(A);
d = find(diff(A) > 1);                  % locations of jumps in A
d = [0 d];                              % prepend for the first nucleotide

A = A * (-2*pi) / mA;                   % convert A to radians

% ---------------------------------- Draw the outside circle

for m = 1:(length(d)-1),
  theta = [A(d(m)+1):-0.01:A(d(m+1)) A(d(m+1))];
  plot(cos(theta),sin(theta),'k');
  hold on
end

for n = 1:length(nc),
  theta = A(nc(n));

  if cos(theta) > 0,
    angle = 180*theta/pi;
    ha    = 'left';
    va    = 'top';
  else
    angle = 180*(theta - pi)/pi;
    ha    = 'right';
    va    = 'bottom';
  end

  text(1.11*cos(theta), 1.11*sin(theta), ['Chain ' File.NT(nc(n)).Chain], 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', va, 'FontSize', 6);

end

if N > 1000,
  step  = 50;
  sstep = 10;
elseif N > 500,
  step = 20;
  sstep = 5;
elseif N > 300,
  step = 10;
  sstep = 2;
elseif N > 100,
  step = 5;
  sstep = 1;
else
  step = 1;
  sstep = 1;
end

kk = 1;

flag = 0;
for k = 1:N,
  kkk = str2num(File.NT(k).Number);
  if ~isempty(kkk),                               % it's really a number
    kk = kkk;
    flag = 0;                                     % OK to use next w/ letter
  else                                            % it has a letter in it
    kkk = str2num(File.NT(k).Number(1:(end-1)));  % omit last character
    if ~isempty(kkk),
      kk = kkk;                                   % use this for display
      flag = 1;                                   % but only once
    end
  end
  theta = A(k);
  if cos(theta) > 0,
    angle = 180*theta/pi;
    ha    = 'left';
  else
    angle = 180*(theta - pi)/pi;
    ha    = 'right';
  end

  if (mod(kk,sstep) == 0) && (mod(kk,step) ~= 0) && (Thickness < 1) && (flag == 0),
    plot([cos(theta) 1.02*cos(theta)], [sin(theta), 1.02*sin(theta)],'k');
    text(1.04*cos(theta), 1.04*sin(theta), File.NT(k).Number,'FontSize',4, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
  end
  if mod(kk,step) == 0 && (flag == 0),
    plot([cos(theta) 1.04*cos(theta)], [sin(theta), 1.04*sin(theta)],'k');
    text(1.11*cos(theta), 1.11*sin(theta), File.NT(k).Number, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
  end

end

% ---------------------------------------- Draw the interactions in color

[y,k] = sort(-abs(c));               % sort by decreasing color, for overlap
%[y,k] = sort(abs(c));               % sort by increasing color, for overlap
i = i(k);
j = j(k);
c = c(k);

color = 'bcrgymk';

shift = 2*pi*[0 0 0 0 -0.2 0.2]/mA;

for k = 1:length(i),
  if View(c(k)) > 0,
    sh = shift(c(k));
    zUnitCircleArc([cos(A(i(k))+sh) cos(A(j(k))+sh)], [sin(A(i(k))+sh) sin(A(j(k))+sh)],c(k),Thickness);

%  if mod(i(k),5) == 2,
      thetai = A(i(k));
      thetaj = A(j(k));

% make the radius depend on how many have already been written here; no overlap

    r = 1.005 + 0.01*c(k);
%    text(r*cos(thetai), r*sin(thetai), num2str(File.Crossing(i(k),j(k))),'FontSize',1, 'Rotation', 0, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle','Color',color(c(k)));
%    text(r*cos(thetaj), r*sin(thetaj), num2str(File.Crossing(i(k),j(k))),'FontSize',1, 'Rotation', 0, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle','Color',color(c(k)));

%[i(k) j(k) thetai thetaj]

%  end

    if c(k) == 6,                 % label base-phosphate interactions
      thetai = A(i(k));
      thetaj = A(j(k));

      if cos(thetai) > 0,
        angle = 180*thetai/pi;
        ha    = 'left';
      else
        angle = 180*(thetai - pi)/pi;
        ha    = 'right';
      end

      if File.BasePhosphate(i(k),j(k)) > 0,
        text(1.03*cos(thetai), 1.03*sin(thetai), 'B','FontSize',1, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle','Color','m');
%        text(1.03*cos(thetaj), 1.03*sin(thetaj), 'Ph','FontSize',1, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
      else
        text(1.03*cos(thetai), 1.03*sin(thetai), 'Ph','FontSize',1, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
%        text(1.03*cos(thetai), 1.03*sin(thetai), 'B','FontSize',1, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
      end

    end
  end
end
axis equal
axis([-1.2 1.2 -2 1.2]);
axis off

if View(7) > 0,

% Color = (B==1).*(C==0) + 2*(B>1).*(B<14).*(C==0) + 3*(B==1).*(C>0) + 4*(B > 1).*(B < 14) .*(C>0) + 5*(B > 20) .* (B < 25);

  text(-1.2,1.2,File.Filename,'HorizontalAlignment','Left');

  if View(1) > 0,
    cww = length(find(c == 1));
    text(-1.3,-1.4,['Dark blue chords indicate the ' num2str(cww) ' nested Watson-Crick basepairs']);
  end

  if View(3) > 0,
    nonnestcww = length(find(c == 3));
    text(-1.3,-1.5,['Red chords indicate the ' num2str(nonnestcww) ' non-nested Watson-Crick basepairs']);
  end

  if View(2) > 0,
    noncww = length(find(c == 2));
    text(-1.3,-1.6,['Cyan chords indicate the ' num2str(noncww) ' nested non-Watson-Crick basepairs']);
  end

  if View(4) > 0,
    nonnestnoncww = length(find(c == 4));
    text(-1.3,-1.7,['Green chords indicate the ' num2str(nonnestnoncww) ' non-nested non-Watson-Crick basepairs']);
  end

  if View(5) > 0,
    stack = length(find(c == 5));
    text(-1.3,-1.8,['Yellow chords indicate the ' num2str(stack) ' stacking interactions']);
  end

  if View(6) > 0,
    bph = length(find(c == 6));
    text(-1.3,-1.9,['Magenta chords indicate the ' num2str(bph) ' base-phosphate interactions']);
  end
end

end
