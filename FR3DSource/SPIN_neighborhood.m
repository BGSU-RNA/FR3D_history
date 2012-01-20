
% SPIN_neighborhood_approximate(D) implements the approximate algorithm from the article referenced below, to find a permutation which, when applied to the rows and columns, concentrates the small entries of the symmetric nonnegative matrix D near the diagonal.  In so doing, it finds a natural ordering of points.

% Sorting points into neighborhoods (SPIN): data analysis and visualization by ordering distance matrices

% D. Tsafrir 1,*, I. Tsafrir 1, L. Ein-Dor 1, O. Zuk 1, D.A. Notterman 2 and E. Domany 1

% http://bioinformatics.oxfordjournals.org/cgi/content/full/21/10/2301

function [p,W] = SPIN_neighborhood(D,Lab,W)

[s,t] = size(D);

% --------------------------- refine permutations

p = 1:s;                      % initial permutation
M = D(p,p)*W;                 % score current D against weight matrix

previous_score = 0;           % previous score
current_score = 1;            % current score, to be updated

while current_score ~= previous_score,

  [i,c] = assignmentoptimal(M);

%i'
%  c
%M

  [z,q] = sort(i);            % permutation q tells how to sort i
  p = p(q);                   % compose previous permutation with this

  M = D(p,p)*W;               % score current D against weight matrix
  previous_score = current_score;
  current_score = trace(M);   % sum down the diagonal

%  subplot(2,2,4);
  zClusterGraph(D,Lab,5,p,0);
  shading flat
  axis ij
  drawnow

  current_score

end



