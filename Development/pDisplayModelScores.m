% pDisplayModelScores displays the results of pMakeModelsFromLibrary
% S(i,2k-1) tells the score that model k gave to sequences from model i
% S(i,2k)   tells the score that model k gave to reversed sequences from i
% We would hope that S(i,2i-1) would be the largest score in row i.

[s,t] = size(S);                         % s by 2s matrix

% ---------------------------------------- Read model names

n = 1;

fid = fopen(['models' filesep loopType '_Models.txt']);
if fid > 0,
  L = 1;
  while L > -1,
    L = fgetl(fid);                      % read a line
    if L > -1,
      Names{n} = L;
      switch loopType
      case 'IL'
        Names{n+1} = ['R ' L];
        n = n + 2;
      case 'HL'
        n = n + 1;
      end
    end
  end
else
  fprintf('Could not open file of model names\n');
end

% ------------------------- Recast scores as relative to max across all models

RelSeq = [];
RelModel = [];

for i = 1:s,
  RelSeq(i,:) = S(i,:) - max(S(i,:));    % subtract the maximum; make it zero
  SLab{i} = Names{2*i-1}(1:8);           % short labels
end

for j = 1:t,
  RelModel(:,j) = S(:,j) - max(S(:,j));  % subtract the maximum; make it zero
end

for i = 1:s,
  for j = 1:s,
    ii = [2*i-1 2*i];                % indices of forward and reversed
    jj = [2*j-1 2*j];
    D(i,j) = min(min([max([0 0],(S(i,ii)-S(i,jj))) max([0 0],(S(j,jj)-S(j,ii)))]));
  end
end

ModelNames = Names(1:2:(t-1));
p = zOrderbySimilarity(D);

figure(10)
clf
q = [];
for i = 1:s,
  q = [q 2*p(i)-1 2*p(i)];          % how to permute models, keep together
  plot([2*i 2*i+2 2*i+2 2*i 2*i]-1,[i i i+1 i+1 i],'k');
  hold on
end

p = 1:s;
q = 1:t;

T = RelSeq(p,q);
T = RelSeq;
T(s+1,t+1) = 0;
pcolor(T);
hold on
axis ij
shading flat
axis([1 2*s+1 1 s+1]);
caxis([-8 0]);
colormap('default')
map = colormap;
map = map(8:56,:);
colormap(map);
colorbar('eastoutside');
set(gca,'YTick',(1:s)+0.5)
set(gca,'YTickLabel',SLab(p),'FontSize',10)
title('Scores of sequences against models and reversed models, relative to max score for those sequences');
xlabel('Models and reversed models, paired, in the order on the y axis');

figure(11)
clf
p = zOrderbySimilarity(D);
zGraphDistanceMatrix(abs(D(p,p)),ModelNames(p),8);
caxis([0 8])
colormap(map);
colorbar('eastoutside');
title('Similarity of models according to how they score each others sequences');

