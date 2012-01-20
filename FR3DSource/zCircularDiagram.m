% zCircularDiagram(Edge,Color) plots the basepairs in Edge using Color as
% chords of a circle

function [void] = zCircularDiagram(File,Edge,Color,Thickness)

if nargin < 2,
  Edge = File.Edge;
end

if nargin < 3,
  Color = sparse(zeros(size(Edge)));
  Color(1,1) = 1;
end

if nargin < 4,
  Thickness = 1;
end

[i,j] = find(Edge);

if length(i) > 0,

for k = 1:length(i),
  c(k) = Color(i(k),j(k));
end

[y,k] = sort(-abs(c));           % sort by decreasing color
i = i(k);
j = j(k);

[s,t] = size(Edge);

tp = -6.2;
theta = 0:-0.01:tp;
plot(cos(theta), sin(theta),'k');
hold on

if s > 1000,
  step  = 50;
  sstep = 10;
elseif s > 500,
  step = 20;
  sstep = 5;
elseif s > 100,
  step = 5;
  sstep = 1;
else
  step = 1;
  sstep = 1;
end

kk = 1;

for k = 1:s,
  kkk = str2num(File.NT(k).Number);
  if ~isempty(kkk),
    kk = kkk;
  end
  theta = k*tp/s;
  if cos(theta) > 0,
    angle = 180*theta/pi;
    ha    = 'left';
  else
    angle = 180*(theta - pi)/pi;
    ha    = 'right';
  end

  if (mod(kk,sstep) == 0) && (mod(kk,step) ~= 0) && (Thickness < 1),
    plot([cos(theta) 1.02*cos(theta)], [sin(theta), 1.02*sin(theta)],'k');
    text(1.04*cos(theta), 1.04*sin(theta), File.NT(k).Number,'FontSize',4, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
  end
  if mod(kk,step) == 0,
    plot([cos(theta) 1.04*cos(theta)], [sin(theta), 1.04*sin(theta)],'k');
    text(1.11*cos(theta), 1.11*sin(theta), File.NT(k).Number, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
  end

end


map = colormap;

colorletter = 'bcrgymbcrgymbcrgymbcrgym';

for k = 1:length(i),
%  c   = 1 + fix(63*(Color(i(k),j(k)))/max(max(Color)));
  if Color(i(k),j(k)) > 0,
    c = colorletter(Color(i(k),j(k)));

    plot([cos(i(k)*tp/s) cos(j(k)*tp/s)], [sin(i(k)*tp/s) sin(j(k)*tp/s)],c,'LineWidth',Thickness);
  end
end
axis square
axis([-1.2 1.2 -1.2 1.2]);
axis off

end
