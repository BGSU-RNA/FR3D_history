
function [void] = SPIN_comparison(D)

[N,M] = size(D);

% --------------------------- define the weight matrix W

c = 1*(N/4)^2;                % standard deviation is N/4
w = exp(-(((1-N):(N-1)).^2)/c);    % large near diagonal, 0 far away
w = [w w];

W = zeros(N,N);

for i = 1:N,
  a = w((N-i+1):(2*N-i));
  W(i,:) = a;        % normalize each row to sum to 1

  Lab{i} = '      ';

end

for i = 1:N,
  W(:,i) = W(:,i) / sum(W(:,i));
end

% ---------------------------- display distance matrices

figure(1)
clf

subplot(2,2,2)
zClusterGraph(D,Lab,[5 2],1:N,0);
axis ij
shading flat
title('Original distance matrix');

subplot(2,2,3)
zClusterGraph(D,Lab,[5 2],[],0);
axis ij
shading flat
title('Distance matrix from zClusterGraph');

p = SPIN_neighborhood(D,Lab,W);

subplot(2,2,4)
zClusterGraph(D,Lab,[5 2],p,0);
axis ij
shading flat
title('Distance matrix from SPIN\_neighborhood');
