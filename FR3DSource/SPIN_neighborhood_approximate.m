
% SPIN_neighborhood_approximate(D) implements the approximate algorithm from the article referenced below, to find a permutation which, when applied to the rows and columns, concentrates the small entries of the symmetric nonnegative matrix D near the diagonal.  In so doing, it finds a natural ordering of points.

% Sorting points into neighborhoods (SPIN): data analysis and visualization by ordering distance matrices

% D. Tsafrir 1,*, I. Tsafrir 1, L. Ein-Dor 1, O. Zuk 1, D.A. Notterman 2 and E. Domany 1

% http://bioinformatics.oxfordjournals.org/cgi/content/full/21/10/2301

function [p] = SPIN_neighborhood_approximate(D,Lab)

[s,t] = size(D);

% --------------------------- define the weight matrix W

c = 1*(s/4)^2;                % standard deviation is s/4
w = exp(-(((1-s):(s-1)).^2)/c);    % large near diagonal, 0 far away
w = [w w];

format long
w'
format short

W = zeros(s,s);

for i = 1:s,
  a = w((s-i+1):(2*s-i));
  W(i,:) = a;        % normalize each row to sum to 1
end

for i = 1:s,
  W(:,i) = W(:,i) / sum(W(:,i));
end

% --------------------------- refine permutations

p = 1:s;                      % initial permutation
M = D(p,p)*W;                 % score current D against weight matrix

previous_score = 0;           % previous score
current_score = 1;            % current score, to be updated

while current_score ~= previous_score,

  [y,i] = min(M,[],2);           % i is position of minimum in each row
i'
%W
M
D(p,p)
  [z,q] = sort(i+0.1*rand(s,1)); % permutation q tells how to sort i
  p = p(q);                   % compose previous permutation with this

%D(p,p)

  M = D(p,p)*W;               % score current D against weight matrix
  current_score = trace(M);   % sum down the diagonal

  subplot(2,2,3);
  zClusterGraph(D,Lab,[5 2],p,0);
  shading flat
  axis ij

current_score
  pause

end



