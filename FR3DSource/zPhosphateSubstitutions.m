
[n,t] = xlsread('Seq_Freqs_for_Conserved_Ec_BPhs_JS.xls');

t = t(3:end,:);

[nn,tt] = xlsread('Seq_Freqs_for_Conserved_Ec_BPhs_JS_23S.xls');

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
types = {'0BPh','1BPh','2BPh','3BPh','4BPh','5BPh','6BPh','7BPh','8BPh','9BPh'};
figure(1)
clf
for t = 1:length(types),
  p = find(ismember(tt(:,4),types{t}));  
  counts = sum(nn(p,4:7),2);
  L = nn(p,4)*A + nn(p,5)*C + nn(p,6)*G + nn(p,7)*U ;
  L = L ./ (counts * [1 1]);
  subplot(3,4,t);
  for i = 1:4,
    j = find(cat(1,tt{p,2}) == Letter{i});
    scatter(L(j,1),L(j,2),10,color(i,:),'filled');
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
  title([num2str(length(p)) ' instances of ' types{t}]);
end

saveas(gcf,['Phosphate Interactions\PhosphateSubstitutions.pdf'],'pdf');
