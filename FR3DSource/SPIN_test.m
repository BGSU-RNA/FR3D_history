
N = 100;

pattern = 5;

switch pattern,
case 1,
  x = rand(N,1);
  y = 0.1*rand(N,1);
case 2,
  theta = rand(N,1)*6.28;
  x = cos(theta);
  y = sin(theta);
case 3,
  x = rand(N,1);
  y = (rand(N,1) > 0.5)*2;
case 4,
  theta = (2*pi/3)*floor(rand(N,1)*3) + (pi/20)*rand(N,1);
  r = rand(N,1);
  x = r .* cos(theta);
  y = r .* sin(theta);
case 5,
  theta = rand(N,1)*6.28;
  r = (rand(N,1) > 0.5) + 0.2*rand(N,1);
  x = r .* cos(theta);
  y = r .* sin(theta);
end

for i = 1:length(x),
  Lab{i} = sprintf('%7.4f', x(i));
end

D = zDistance([x y],[x y]);

% --------------------------- define the weight matrix W

c = 1*(N/4)^2;                % standard deviation is N/4
w = exp(-(((1-N):(N-1)).^2)/c);    % large near diagonal, 0 far away
w = [w w];

format long
w'
format short

W = zeros(N,N);

for i = 1:N,
  a = w((N-i+1):(2*N-i));
  W(i,:) = a;        % normalize each row to sum to 1
end

for i = 1:N,
  W(:,i) = W(:,i) / sum(W(:,i));
end

% ---------------------------- display distance matrices

figure(1)
clf

subplot(2,2,1)
plot(x,y,'.');
title('Data points');

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
