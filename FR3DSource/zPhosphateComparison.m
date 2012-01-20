% zPhosphateComparison takes phosphate classifications from multiple files and calculates average phosphate interaction values

% File = zAddNTData('1s72_equiv');             % load data
% [File,D] = zPhosphateInteractions(File,2,0); % analyze interactions
% save PhosphateComparisonD D

load PhosphateComparisonD

i = find(D(:,17) == 1);                        % use best oxygen only
D = D(i,:);

D(:,23) = mod(D(:,5),100);                      % don't distinguish near BPh

D = sortrows(D,5);                              % sort by classification

i = zeros(size(D(:,1)));

% --------------------------------------------- restrict to chains O, A only

OK = '0A';

for m = 1:length(D(:,1)),
  if any(File(D(m,1)).NT(D(m,2)).Chain == OK) && any(File(D(m,1)).NT(D(m,3)).Chain == OK),
    i(m) = 1;
  end
end

D = D(find(i),:);

[s,t] = size(D);

n1n2 = unique(D(:,[21 22 23]), 'rows');         % unique BPh interactions

clear E

for m = 1:length(n1n2(:,1)),                   % loop through unique inters
  i = find((D(:,21) == n1n2(m,1)) .* (D(:,22) == n1n2(m,2)) .* (D(:,23) == n1n2(m,3)));

%  i = find((D(:,21) == n1n2(m,1)) .* (D(:,22) == n1n2(m,2)) );

%length(i)
%D(i,[2 3 5])
%pause  

  if length(i) > 1,
    E(m,1:t) = mean(D(i,:));                   % average parameters
  else
    E(m,1:t) = D(i,:);                         % unless there is just one
  end
  E(m,1)   = D(i(1),1);                        % file number of first
  E(m,2)   = D(i(1),2);
  E(m,3)   = D(i(1),3);
  E(m,5)   = D(i(1),5);

  E(m,t+1) = length(i);                        % store number of such inters
  E(m,t+2) = std(D(i,9));                      % variation in bond length
  E(m,t+3) = std(D(i,8));                      % variation in angle
end

% --------------------------------------- remove cases of only near BPh

i = find(E(:,5) < 100);                          % first happens to be true BP
E = E(i,:);                                      % use these only


figure(1)
clf
n = hist(E(:,t+1),0.5:1:(length(File)-0.5));
hist(E(:,t+1),0.5:1:(length(File)-0.5));
xlabel('Number of structures with each base-phosphate interaction (near or true)');
axis([0 length(File) 0 max(n)*1.1]);

figure(2)
clf
scatter(E(:,9),E(:,8),6,E(:,t+1),'filled');
caxis([1 1.1*max(D(:,1))]);
colorbar('eastoutside');
title('BPh parameters averaged over corresponding instances')
xlabel('Colored by number of structures having the same interaction');

figure(3)
clf
scatter(E(:,9),E(:,8),6,E(:,t+2),'filled');
colorbar('eastoutside');
title('BPh parameters averaged over corresponding instances')
xlabel('Colored by standard deviation of bond length');

figure(4)
clf
scatter(E(:,9),E(:,8),6,E(:,t+3),'filled');
colorbar('eastoutside');
title('BPh parameters averaged over corresponding instances')
xlabel('Colored by standard deviation of angle');

figure(5)
clf
plot(E(:,23)+0.5*rand(size(E(:,23)))-0.25,E(:,t+1),'.');
title('Which interactions are conserved between structures');
xlabel('Internal interaction code');
ylabel('Number of structures');

figure(6)
clf

n = hist(E(:,t+1),0.5:1:(length(File)-0.5));
hist(E(:,t+1),0.5:1:(length(File)-0.5));
xlabel('Number of structures with each base-phosphate interaction');
axis([0 length(File) 0 max(n)*1.1]);

break

% ------------------------------------- Display nucleotides in all structures
% ------------------------------------- when one has a true Base-Phosphate

% -------- Note: with 1s72_equiv, the large chain is 0 in some, A in others
% -------- The small chain is 9 or B.  So it is hard to find all the 
% -------- nucleotides in all structures.
% -------- Above, we include only chain 0 or A
% -------- Below, we look up the index using only nucleotide number

i = find(E(:,5) < 100);                          % first happens to be true BP
E = E(i,:);                                      % use these only

E = sortrows(E,[t+1 21 22]);                     % sort by number of inters

for r = 1:length(E(:,1)),                        % loop through instances
  Cand = [];

  for f = 1:length(File),                        % loop through files
    i = zIndexLookup(File(f),File(E(r,1)).NT(E(r,2)).Number,'0',0);
    i = [i zIndexLookup(File(f),File(E(r,1)).NT(E(r,2)).Number,'A',0)];

    j = zIndexLookup(File(f),File(E(r,1)).NT(E(r,3)).Number,'0',0);
    j = [j zIndexLookup(File(f),File(E(r,1)).NT(E(r,3)).Number,'A',0)];
    if ~isempty(i) && ~isempty(j),
      Cand = [Cand; [i j f]];
    end
  end

  xDisplayCandidates(File,Cand);
end


break

for f = 1:length(File),
  [b,i,j] = unique(cat(1,File(f).NT.Chain));
  fprintf('File %s has ', File(f).Filename);
  for m = 1:length(b),
    fprintf('%d in chain %s, ', sum(j==m), b(m));
  end
  fprintf('\n');
end
