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
%  9 = display nucleotide numbers in the diagram
% 10 = display nucleotide letters in the diagram
% 11 = display B for BPh interactions in the diagram

function [Tally] = zCircularDiagram(File,Thickness,View)

Tally = zeros(1,6);

if nargin < 2,
  Thickness = 1;
end

if nargin < 3,
  View = [1 1 1 1 1 1 1 1 1 1 1];
end

while length(View) < 9,
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

Color = (B==1).*(C==0) + 2*(B>1).*(B<13).*(C==0) + 3*(B==1).*(C>0) + 4*(B > 1).*(B < 13) .*(C>0) + 5*(B > 20) .* (B < 25);
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
  [A,mA] = zNumberCircularDiagram(File,View,Thickness,1);

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

      if File.BasePhosphate(i(k),j(k)) > 0 && View(11) > 0,
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


text(-1.2,1.2,File.Filename,'HorizontalAlignment','Left');


if View(7) > 0,

% Color = (B==1).*(C==0) + 2*(B>1).*(B<13).*(C==0) + 3*(B==1).*(C>0) + 4*(B > 1).*(B < 13) .*(C>0) + 5*(B > 20) .* (B < 25);


  cww = length(find(c == 1));
  Tally(1,1) = cww;
  noncww = length(find(c == 2));
  Tally(1,2) = noncww;
  nonnestcww = length(find(c == 3));
  Tally(1,3) = nonnestcww;
  nonnestnoncww = length(find(c == 4));
  Tally(1,4) = nonnestnoncww;
  stack = length(find(c == 5));
  Tally(1,5) = stack;
  bph = length(find(c == 6));
  Tally(1,6) = bph;

  if View(1) > 0,
    text(-1.3,-1.4,['Dark blue chords indicate the ' num2str(cww) ' nested Watson-Crick basepairs']);
  end

  if View(2) > 0,
    text(-1.3,-1.6,['Cyan chords indicate the ' num2str(noncww) ' nested non-Watson-Crick basepairs']);
  end

  if View(3) > 0,
    text(-1.3,-1.5,['Red chords indicate the ' num2str(nonnestcww) ' non-nested Watson-Crick basepairs']);
  end

  if View(4) > 0,
    text(-1.3,-1.7,['Green chords indicate the ' num2str(nonnestnoncww) ' non-nested non-Watson-Crick basepairs']);
  end

  if View(5) > 0,
    text(-1.3,-1.8,['Yellow chords indicate the ' num2str(stack) ' stacking interactions']);
  end

  if View(6) > 0,
    text(-1.3,-1.9,['Magenta chords indicate the ' num2str(bph) ' base-phosphate interactions']);
  end
end

end
