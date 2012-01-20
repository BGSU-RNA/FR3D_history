
clear S

pModelScores;                            % read model scores copied from JAR3D

[s,t] = size(S);                         % s by 2s matrix

for i = 1:s,
  S(i,:) = S(i,:) - max(S(i,:));         % subtract the maximum; make it zero
  SLab{i} = Names{2*i-1}(1:8);
end

for i = 1:s,
  for j = 1:s,
    ii = [2*i-1 2*i];                % indices of forward and reversed
    jj = [2*j-1 2*j];
    D(i,j) = min(min([max([0 0],(S(i,ii)-S(i,jj))) max([0 0],(S(j,jj)-S(j,ii)))]));
  end
end



ModelNames = Names(1:2:(t-1));

figure(2)
p = zOrderbySimilarity(D);
zGraphDistanceMatrix(abs(D(p,p)),ModelNames(p),8);
caxis([0 8])
colorbar('eastoutside');
title('Similarity of models according to how they score each others sequences');

figure(1)
q = [];
for i = 1:s,
  q = [q 2*p(i)-1 2*p(i)];
end
pcolor(S(p,q));
axis ij
shading flat
caxis([-8 0]);
colormap('default')
map = colormap;
map = flipud(map);
colormap(map);
colorbar('eastoutside');
set(gca,'YTick',(1:s)+0.5)
set(gca,'YTickLabel',SLab(p),'FontSize',10)
title('Scores of sequences against models and reversed models, relative to max score');
xlabel('Models and reversed models, paired, in the order on the y axis');
