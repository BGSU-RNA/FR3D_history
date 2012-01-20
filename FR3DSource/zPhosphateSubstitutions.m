% zPhosphateSubstitutions reads Jesse's alignments of 16S and base-phosphate interactions and substitution data from sequence alignments, then makes graphs for each base-phosphate interaction category of the pattern of substitutions

[n,t] = xlsread('16S_Ec_Tt_BPh_Aln_12_6_08_Seq_Data.xls');
t = t(3:end,:);

[nn,tt] = xlsread('23S_Ec_Tt_BPh_Aln_12_6_08_Seq_Data.xls');
tt = tt(3:end,:);

tt = [t; tt];
nn = [n; nn];

[a,b] = size(nn);

for i = 1:a,
  for j = 1:b,
    if isnan(nn(i,j)),
      nn(i,j) = 0;
    end
  end
end

A = [0 1];
C = [1 1];
G = [1 0];
U = [0 0];
Letter = {'A','C','G','U'};
color = [[1 0 0]; [0.9 0.9 0]; [0 0.9 0]; [0 0 1]];
types = {'0BPh','1BPh','2BPh','3BPh','4BPh','4bBPh','5BPh','6BPh','7BPh','8BPh','8bBPh','9BPh'};
figure(1)
clf
for t = 1:length(types),

  subplot(3,4,t);

  % ----------------------- aligned positions plotted as dots
  p = find(ismember(tt(:,10),types{t}));                % rows with OK BPh type
  i = zeros(length(p),1);
  for j = 1:length(p),
    if ~isempty(tt{p(j),19}),
      i(j) = 1;
    end
  end
  p = p(find(i));
  counts = sum(nn(p,22:25),2);                          % sequence counts
  L = nn(p,22)*A + nn(p,23)*C + nn(p,24)*G + nn(p,25)*U ;  % sequence counts
  L = L ./ (counts * [1 1]);                            % weighted average
  for i = 1:4,
    j = find(cat(1,tt{p,5}) == Letter{i});              % base in 3D structure
    scatter(L(j,1),L(j,2),10,color(i,:),'filled');      % color differently
    hold on
  end
  
  pp = length(p);

  % ----------------------- non-aligned positions plotted as plusses
  p = find(ismember(tt(:,10),types{t}));                % BPh type
  i = zeros(length(p),1);
  for j = 1:length(p),
    if isempty(tt{p(j),19}),
      i(j) = 1;
    end
  end
  p = p(find(i));
  counts = sum(nn(p,22:25),2);                          % sequence counts
  L = nn(p,22)*A + nn(p,23)*C + nn(p,24)*G + nn(p,25)*U ;  % sequence counts
  L = L ./ (counts * [1 1]);                            % weighted average
  for i = 1:4,
    j = find(cat(1,tt{p,5}) == Letter{i});              % base in 3D structure
    plot(L(j,1),L(j,2),'+','color',color(i,:));      % color differently
    hold on
  end

  axis([0 1 0 1]);
  axis off
  s = 0.07;
  text(-s,1,'A','HorizontalAlignment','center','Color',color(1,:));
  text(1+s,1,'C','HorizontalAlignment','center','Color',color(2,:));
  text(1+s,0,'G','HorizontalAlignment','center','Color',color(3,:));
  text(-s,0,'U','HorizontalAlignment','center','Color',color(4,:));
  plot([0 1 1 0 0], [0 0 1 1 0], 'k');
  title([num2str(pp) ' (' num2str(length(p)) ') instances of ' types{t}]);
end

saveas(gcf,['Phosphate Interactions\PhosphateSubstitutions.pdf'],'pdf');
