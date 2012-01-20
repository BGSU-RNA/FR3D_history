
N = 1000;

A = rand(N,3);
B = rand(2*N,3);

fprintf('Generated random matrix\n');

t = cputime;

D = zDistance(A,B);

fprintf('zDistance took %8.4f seconds\n', cputime-t);

t = cputime;

E = dist(A,B');

fprintf('dist took %8.4f seconds\n', cputime-t);

fprintf('Maximum difference is %18.18f \n', max(max(D-E)));

