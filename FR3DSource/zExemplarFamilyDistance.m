
load PairExemplars

% Pair codes:  1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU

pcodes = [1 5 6 7 9 11 13 14 15 16];    % pair codes to work on

[s,t] = size(Exemplar);

for p = 1:length(pcodes),
  for q = p:length(pcodes),
    for a = 1:s,                        % loop through exemplars for pcodes(a)
      P1 = Exemplar(a,pcodes(p));
      c1 = P1.Class;
      if c1 == fix(c1),                 % not a subcategory of this family
        for b = 1:s,                    % loop through exemplars for pcodes(b)
          P2 = Exemplar(b,pcodes(q));
          c2 = P2.Class;
          if c2 == fix(c2),             % not a subcategory of this family
            t1 = zIsoDiscrepancy(P1.NT1,P1.NT2,P2.NT1,P2.NT2);
            t2 = zIsoDiscrepancy(P1.NT1,P1.NT2,P2.NT2,P2.NT1);
            dists{c1,c2} = [dists{c1,c2} min(t1,t2)];
          end
        end
      end
    end
  end
end

for c = 1:12,
  D(c,c) = mean(dists{c,c});
  for d = (c+1):12,
    D(c,d) = mean([dists{c,d} dists{d,c}]);
    D(d,c) = mean([dists{c,d} dists{d,c}]);
  end
end


Lab{1}  = 'cWW';
Lab{2}  = 'tWW';
Lab{3}  = 'cWH';
Lab{4}  = 'tWH';
Lab{5}  = 'cWS';
Lab{6}  = 'tWS';
Lab{7}  = 'cHH';
Lab{8}  = 'tHH';
Lab{9}  = 'cHS';
Lab{10} = 'tHS';
Lab{11} = 'tSS';
Lab{12} = 'cSS';



% ----------------------------------------- Cluster analysis and graph

Y = squareform(full(D));                       % convert to a vector

DD = full(D);

Z = linkage(Y,'average');                      % compute cluster tree
figure(fix(Category(ca))+13)
[H,T,p] = dendrogram(Z,0,'colorthreshold',threshold,'orientation','left','labels',Lab);

for j =1:length(Lab),
  fprintf('%s ', Lab{j}(1:2));
end
fprintf('\n');

h = gcf;
hh = gca;
if length(Lab) > 80,
  set(hh,'FontSize',3);
elseif length(Lab) > 40,
  set(hh,'FontSize',6);
else
  set(hh,'FontSize',8);
end
orient landscape
%print(h,'-depsc2',['Isostericity\ClusterIsoDisc' zEdgeText(Category)]);


% ------------------------------- Print table of isodiscrepancies

p = zOrderGroups(Y,Z,DD);

DDD = DD(p,p);                          % re-order according to p
[s,t] = size(DDD);

if length(Category) == 1,
  fprintf('Table %d\n',Category(ca));
else
  fprintf('Tables ');
  for c = 1:length(Category),
    fprintf('%d ', Category(c));
  end
  fprintf('\n');
end
fprintf('               ');
for j=1:t,
  fprintf('%6s ', Lab{p(j)}(1:6));
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

Title = [];
for i = 1:length(Category),
  Title = [Title ' ' zEdgeText(abs(Category(i)))];
end
Title = [Title 'family isodiscrepancy map'];
if Subcat > 0,
  Title = [Title ' with subcategories'];
end
title(Title);

for i = 1:t,
  SLab{i} = Lab{i}(1:6);
  SSLab{i} = Lab{i}(1:2);
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

saveas(gcf,['Isostericity' filesep Title '.pdf'],'pdf');
saveas(gcf,['Isostericity' filesep Title '.png'],'png');
