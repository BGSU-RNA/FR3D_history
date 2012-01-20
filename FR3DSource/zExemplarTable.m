% zExemplarTable(Cateogry) displays the best known representatives for interactions involving all pairs in interaction category(ies) Category
% Example:  zExemplarTable(1,3.5)
% Example:  zExemplarTable(1:12,3.5)

% threshold controls coloring of the dendrogram
% Subcat = 1, include subcategories, Subcat = 0, don't

function [Y,Z,DD] = zExemplarTable(Category,threshold,Subcat,Verbose)

close all

if nargin < 1,
  Category = 1;
end

if nargin < 2,
  threshold = 'default';
end

if nargin < 3,
  Subcat = 1;
end

if nargin < 4,
  Verbose = 2;
end

% load exemplars -------------------------------------

load('PairExemplars','Exemplar');

% loop through computer classifications, accumulating exemplars ------------

B(1) = Exemplar(1,1);
B(1).HydrogenClass = 0;
B(1).subplot  = 0;
B(1).original = 0;

Lab = [];
Cat = [];

% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

for c1 = 1:4,
 for c2 = 1:4,
  pc  = 4*(c2-1)+c1;                       % current paircode
  for r = 1:length(Exemplar(:,1)),         % loop through rows of Exemplar

    E = Exemplar(r,pc);                    % current exemplar

    if ~isempty(E.NT1),                    % non-empty entry of Exemplar
    if  any(abs(E.Class) == Category) || ...
       (any(fix(abs(E.Class)) == Category) && (Subcat == 1)),
    if (E.Count >= 0),
    
      [B,Lab,Cat] = AddExemplar(E,B,Lab,Cat,Subcat);

      % ------- in some symmetric families, produce the symmetric version
      if (any(pc == [5 7 9 13 14 15])) && any(fix(E.Class) == [1 2 7 8]),
        E.Class = -E.Class;
        [B,Lab,Cat] = AddExemplar(E,B,Lab,Cat,Subcat);
        E.Class = -E.Class;
      end

      % ------- in some symmetric families, store AA, CC, GG, UU pairs twice

      if (E.NT1.Code == E.NT2.Code) && any(fix(E.Class) == [1 7 8 11 12]),
        E.HydrogenClass = 0;
        E.subplot  = 0;
        E.original = 0;
        E.Class = -E.Class;
        [B,Lab,Cat] = AddExemplar(E,B,Lab,Cat,Subcat);
      end


    end
    end
    end
  end
 end
end

% -------------------------------------- Sort exemplars by category

B   = B(2:end);
Lab = Lab(2:end);
Cat = Cat(2:end);

for m = 1:length(B),
  categ(m) = abs(B(m).Class);
end

[y,i] = sort(categ);

B = B(i);
Lab = Lab(i);
Cat = Cat(i);

% specify parameters for viewing -------------------------------------------

ViewParam.Mode      = 1; 
ViewParam.Normal    = 1;
ViewParam.ColorAxis = [-12 30];
ViewParam.SortKeys  = [];
ViewParam.Nearby    = 0;
ViewParam.Sugar     = 1;
ViewParam.ConnectSugar = 0;
ViewParam.AtOrigin  = 1;
ViewParam.Hydrogen  = 1;
ViewParam.Sort      = 0;
ViewParam.LabelBases= 8;                    % font size

% -------------------------------------- Plot exemplars

if Verbose > 1,
  for ca = 1:length(Category),
    figure(fix(Category(ca)))
    clf
    plotted = zeros(16,1);                   % keep track of which are plotted
    for m = 1:length(B),
     if (B(m).original == 1) && ((abs(B(m).Class) == Category(ca)) ...
        || ((abs(fix(B(m).Class)) == Category(ca)) && (Subcat == 1))),

        pc2 = B(m).subplot;
        E   = B(m);
        % display the exemplar pair -----------------------------------------

       if abs(E.Class - fix(E.Class)) == 0,
         ViewParam.LineStyle = '-';
       elseif abs(E.Class - fix(E.Class)) > 0.29,
         ViewParam.LineStyle = '.';
       elseif abs(E.Class - fix(E.Class)) > 0.19,
         ViewParam.LineStyle = ':';
       elseif abs(E.Class - fix(E.Class)) > 0.09,
         ViewParam.LineStyle = '--';
       end

       subplot(4,4,pc2);

       if plotted(pc2) > 0,
         ViewParam.LabelBases = 0;                  % don't display now
         xlab{pc2} = [xlab{pc2} ', ' num2str(E.Count)];  % append subcat count
       else,
         ViewParam.LabelBases = 8;                  % font size
         xlab{pc2} = ['Count: ' num2str(E.Count)];
         Title{pc2} = [E.NT1.Base E.NT2.Base zEdgeText(E.Class,Subcat,E.NT1.Code,E.NT2.Code) ' ' strrep(E.Filename,'_','\_') ' '];
         Title{pc2} = [Title{pc2} E.NT1.Base E.NT1.Number '-' E.NT2.Base E.NT2.Number];
         CP = norm(E.NT1.Sugar(1,:) - E.NT2.Sugar(1,:));     % c1'-c1' dist
         Title{pc2} = [Title{pc2} ' ' num2str(CP)];
       end

       F.NT(1) = E.NT1;
       F.NT(2) = E.NT2;
       F.Filename = E.Filename;
       zDisplayNT(F,[1 2],ViewParam);
       zPlotHydrogenBonds(E.NT1,E.NT2,E.HydrogenClass,E.NT1.Rot,E.NT1.Fit(1,:));

       view(2)
       grid off
       axis equal
       %axis tight
       axis fill
       xlabel(xlab{pc2});
       title(Title{pc2});

%       a = axis;
%       text(0.3*a(1)+0.7*a(2),0.8*a(3)+0.2*a(4),['Count: ' num2str(E.Count)]);

       rotate3d on
       plotted(pc2) = 1;

     end
    end
  end

 h = gcf;
 orient landscape
% print(h,'-depsc2',['Isostericity\Exemplars' zEdgeText(Category(ca))]);

end

% -------------------------------------- Compare basepairs against each other

if exist('B'),
D = [];
G = zeros((length(B)^2-length(B))/2,7);
j = 1;

for a = 1:length(B),
  for b = (a+1):length(B),

    % calculate isodiscrepancy

    D(a,b) = zIsoDiscrepancy(B(a).NT1,B(a).NT2,B(b).NT1,B(b).NT2);

if B(a).Class == B(b).Class,
%  LookupTable{B(a).Class}(B(a).NT1.Code,B(a).NT2.Code,B(b).NT1.Code,B(b).NT2.Code) = D(a,b);
%  LookupTable{B(a).Class}(B(a).NT2.Code,B(a).NT1.Code,B(b).NT2.Code,B(b).NT1.Code) = D(a,b);
end
  

    D(b,a) = D(a,b);

    RD(a,b) = zIsoDiscrepancy(B(a).NT1,B(a).NT2,B(b).NT2,B(b).NT1);
    RD(b,a) = RD(a,b);               % reverse order, for family distance

    j = j + 1;

  end
end

% ----------------------------------------- Cluster analysis

Y = squareform(full(D));                       % convert to a vector
Z = linkage(Y,'average');                      % compute cluster tree

% ----------------------------------------- Display dendrogram

figure
[H,T,p] = dendrogram(Z,0,'colorthreshold',threshold,'orientation','left','labels',Lab);

h = gcf;
hh = gca;
if length(Lab) > 80,
  set(hh,'FontSize',3);
elseif length(Lab) > 40,
  set(hh,'FontSize',6);
else
  set(hh,'FontSize',8);
end
%orient landscape
%print(h,'-depsc2',['Isostericity\ClusterIsoDisc' zEdgeText(Category)]);

% ------------------------------------------ Print table of isodiscrepancies
% ------------------------------------------ Display graph of isodiscrepancies

zClusterGraph(D, Lab, [12 2], [], 0);

Title = [];
for i = 1:length(Category),
  Title = [Title ' ' zEdgeText(abs(Category(i)))];
end
Title = [Title 'family isodiscrepancy map'];
if Subcat > 0,
  Title = [Title ' with subcategories'];
end
title(Title);

saveas(gcf,['Isostericity' filesep Title '.pdf'],'pdf');
saveas(gcf,['Isostericity' filesep Title '.png'],'png');

% ----------------------------------------- Calculate avg dist btw families

if length(Category) > 1,

Fam{1}  = 'cWW';
Fam{2}  = 'tWW';
Fam{3}  = 'cWH';
Fam{4}  = 'tWH';
Fam{5}  = 'cWS';
Fam{6}  = 'tWS';
Fam{7}  = 'cHH';
Fam{8}  = 'tHH';
Fam{9}  = 'cHS';
Fam{10} = 'tHS';
Fam{11} = 'tSS';
Fam{12} = 'cSS';

MeanD = [];

minD = min(D,RD);

for p = 1:12,
  i = find(ismember(Cat,upper(Fam{p})));
  if ~isempty(i),
    for q = 1:12,
      j = find(ismember(Cat,upper(Fam{q})));
      if ~isempty(j),
%        MeanD(p,q) = mean(mean(minD(i,j)));
        MeanD2(p,q) = mean(mean(D(i,j)));
        MeanD3(p,q) = mean(mean(RD(i,j)));
        MeanD(p,q)  = min(MeanD2(p,q),MeanD3(p,q));
      end
    end
  end
end

zClusterGraph(MeanD,Fam,[3 3]);

end

return

% -------------------------------------------------------------------------
% ----------------------------------------- Display nearby exemplars

% ------------ display in dendrogram order

[s,t] = size(DDD);
if length(Category) == 1,
  fprintf('Table %d\n',Category(1));
else
  fprintf('Tables ');
  for c = 1:length(Category),
    fprintf('%d ', Category(c));
  end
  fprintf('\n');
end

for i=1:s,
  fprintf('%24s ', Lab{p(i)});
  [d,h] = sort(DDD(i,:));
  a = 2;
  d = [d 99999];

  while (d(a) < 10),
    fprintf('%24s %5.2f ', Lab{p(h(a))}, d(a));
    a = a + 1;
  end

  fprintf('\n');
end
fprintf('\n');

% ------------ display in zClassLimits order

[s,t] = size(DD);
if length(Category) == 1,
  fprintf('Table %d\n',Category(1));
else
  fprintf('Tables ');
  for c = 1:length(Category),
    fprintf('%d ', Category(c));
  end
  fprintf('\n');
end

for i=1:s,
  fprintf('%24s ', Lab{i});
  [d,h] = sort(DD(i,:));
  a = 2;
  d = [d 99999];

  while (d(a) < 10),
    fprintf('%24s %5.2f ', Lab{h(a)}, d(a));
    a = a + 1;
  end
  fprintf('\n');
end
fprintf('\n');



return

% ----------------------------------------- Components of isodiscrepancy

Lab{1} = 'Angle';
Lab{2} = 'Vector diff 1';
Lab{3} = 'Vector diff 2';
Lab{4} = 'C1* difference';

g = 4;
figure(fix(Category(ca))+26)
clf
for i=2:g,
  for j=1:(i-1),
    subplot(g-1,g-1,(i-2)*(g-1)+j);
    plot(G(:,i),G(:,j),'.');
    xlabel(Lab{i});
    ylabel(Lab{j});
  end
end

end

% ----------------------------------------- Display basepairs by IsoDisc

fprintf('Comparison of basepairs, lowest IsoDiscrepancy first\n')

G = sortrows(G,7);
for i = 1:length(G(:,1)),
  a = G(i,5);
  b = G(i,6);
  figure(fix(Category(ca))+39)
  clf
  F.NT(1) = B(a).NT1;
  F.NT(2) = B(a).NT2;
  F.Filename = B(a).Filename;
  subplot(1,2,1)
  ViewParam.LineStyle = '-';
  zDisplayNT(F,[1 2],ViewParam);
  hold on

  subplot(1,2,2)
  hold on
  ViewParam.LineStyle = '--';
  zDisplayNT(F,[2 1],ViewParam);

  F.NT(1) = B(b).NT1;
  F.NT(2) = B(b).NT2;
  F.Filename = B(b).Filename;
  subplot(1,2,1)
  ViewParam.LineStyle = '--';
  zDisplayNT(F,[1 2],ViewParam);
  hold on
  view(2)
  subplot(1,2,2)
  rotate3d on
  ViewParam.LineStyle = '-';
  zDisplayNT(F,[2 1],ViewParam);
  hold on
  xlab = ['Count: ' num2str(B(b).Count)];
  Title = [B(b).NT1.Base B(b).NT2.Base zEdgeText(B(b).Class,Subcat,B(b).NT1.Code,B(b).NT2.Code) ' ' num2str(B(b).Pair.Class,3) ' ' strrep(B(b).Filename,'_','\_') ' '];
  Title = [Title B(b).NT1.Base B(b).NT1.Number '-' B(b).NT2.Base B(b).NT2.Number];
  CP = norm(B(b).NT1.Sugar(1,:) - B(b).NT2.Sugar(1,:));     % c1'-c1' dist
  Title = [Title ' ' num2str(CP)];
  xlabel(xlab);
  title(Title);
  view(2)

  subplot(1,2,1)
  rotate3d on
  xlab = ['Count: ' num2str(B(a).Count)];
  Title = [B(a).NT1.Base B(a).NT2.Base zEdgeText(B(a).Class,Subcat,B(a).NT1.Code,B(a).NT2.Code) ' ' num2str(B(a).Pair.Class,3) ' ' strrep(B(a).Filename,'_','\_') ' '];
  Title = [Title B(a).NT1.Base B(a).NT1.Number '-' B(a).NT2.Base B(a).NT2.Number];
  CP = norm(B(a).NT1.Sugar(1,:) - B(a).NT2.Sugar(1,:));     % c1'-c1' dist
  Title = [Title ' ' num2str(CP)];
  xlabel(xlab);
  title(Title);

  fprintf('IsoDiscrepancy %6.3f\n', G(i,7));
  fprintf('Angle contribution          %6.3f\n', G(i,1));
  fprintf('C1* location contribution 1 %6.3f\n', G(i,2));
  fprintf('C1* location contribution 2 %6.3f\n', G(i,3));
  fprintf('C1* distance contribution   %6.3f\n', G(i,4));
  fprintf('Press a key\n');

  pause
  
end
  
return

zExemplarTable(1,3.5,0,0);
zExemplarTable(2,3.5,0,0);
zExemplarTable(3,3.5,0,0);
zExemplarTable(4,3.5,0,0);
zExemplarTable(5,3.5,0,0);
zExemplarTable(6,3.5,0,0);
zExemplarTable(7,3.5,0,0);
zExemplarTable(8,3.5,0,0);
zExemplarTable(9,3.5,0,0);
zExemplarTable(10,3.5,0,0);
zExemplarTable(11,3.5,0,0);
zExemplarTable(12,3.5,0,0);
zExemplarTable(13,3.5,0,0);
zExemplarTable([1:12],3.5,0,0);
zExemplarTable([1 5],3.5,0);
zExemplarTable(1:6,3.5,0);
zExemplarTable(7:12,3.5,0);


% --------------------------------------------------------------------------

function [B,Lab,Cat] = AddExemplar(E,B,Lab,Cat,Subcat)

m = length(B) + 1;

pc = 4*(E.NT2.Code-1)+E.NT1.Code;                  % paircode

E.HydrogenClass = E.Class;

if (E.Class < 0),
  T     = E.NT1;
  E.NT1 = E.NT2;
  E.NT2 = T;
  if E.NT1.Code ~= E.NT2.Code,
    E.Class = -E.Class;
  end
end

if (pc == 7),                                      % reverse GC pairs
  T     = E.NT1;
  E.NT1 = E.NT2;
  E.NT2 = T;
end

pc2 = 4*(E.NT1.Code-1)+E.NT2.Code;                 % paircode when reversed

E.subplot = pc2;
if ~isfield(E,'original'),
  E.original = 1;
end

ic = zIsostericSubgroups(E.NT1.Code,E.NT2.Code,abs(E.Class));

B(m) = E;                    % store this exemplar for isodisc calc
Lab{m} = [E.NT1.Base E.NT2.Base zEdgeText(E.Class,Subcat,E.NT1.Code,E.NT2.Code) ' ' ic sprintf('%5.1f',E.Pair.Class) ' ' E.Filename ' ' sprintf('%4d',E.Count)];
Cat{m} = upper(zEdgeText(E.Class,Subcat,E.NT1.Code,E.NT2.Code));
