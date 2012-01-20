% zClusterGraph(D,Lab,W) treats D as the distances between instances with labels given by Lab.  It re-orders the instances to group them into clusters, and puts nearby clusters next to one another.  Then it displays the distances graphically.  W(1) is the number of characters of Lab to use on the vertical axis, W(2) on the horizontal.  pp is an optional user-supplied ordering vector.

function [void] = zClusterGraph(D,Lab,W,pp)

% ----------------------------------------- Cluster analysis

d = diag(D);                                   % save diagonal for later

for i = 1:length(D(:,1)),
  D(i,i) = 0;                                  % set diagonal to zero
end

Y = squareform(full(D));                       % convert to a vector
Z = linkage(Y,'average');                      % compute cluster tree

% ------------------------------------------ Print table

DD = full(D);
p = zOrderGroups(Y,Z,DD);               % put instances in order

if nargin == 4,
  p = pp;                               % use user-supplied ordering
end

DDD = DD(p,p);                          % re-order according to p
d = d(p);                               % re-order diagonal too

for i = 1:length(DDD),
  DDD(i,i) = d(i);
end

[s,t] = size(DDD);

fprintf('                      ');
for j=1:t,
  w = min(length(Lab{p(j)}), W);
  fprintf('%6s', Lab{p(j)}(1:w));
end
fprintf('\n');
for i=1:s,
  fprintf('%24s ', Lab{p(i)});
  for j= 1:t,
    fprintf('%5.2f ', DDD(i,j));
  end
  fprintf('\n');
end
fprintf('\n');

% ------------------------------------------ Display graph of isodiscrepancies

DDDD = zeros(s+1,t+1);
DDDD(1:s,1:t) = DDD;
figure
pcolor(DDDD)
shading flat
axis ij
view(2)
colormap('default');
map = colormap;
map = map((end-8):-1:8,:);
%map = map((end-8):-1:end,:);
colormap(map);
caxis([0 16]);
colorbar('location','eastoutside');

for i = 1:t,
  SLab{i} = Lab{i}(1:W(1));
  SSLab{i} = Lab{i}(1:W(2));
end

FS = 12;
if length(Lab) > 20,
  FS = 8;
elseif length(Lab) > 50,
  FS = 4;
elseif length(Lab) > 100,
  FS = 3;
end

set(gca,'XTick',(1:s)+0.5)
set(gca,'XTickLabel',SSLab(p),'FontSize',FS)

set(gca,'YTick',(1:s)+0.5)
set(gca,'YTickLabel',SLab(p),'FontSize',FS)
