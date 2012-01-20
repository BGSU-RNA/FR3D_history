% xSortByCentrality re-orders the candidates according to their centrality

function [Search] = xSortByCentrality(File,Search)

Search = xMutualDiscrepancy(File,Search);

s = length(find(Search.DiscComputed));            % number computed

% ----------------------------------- Sort by centrality

[z,j] = sort(sum(Search.Disc));
z = z / (s-1);                                % average discrepancy among these

% ----------------------------------- List and display results

fprintf('Candidates sorted by centrality within these candidates:\n');
fprintf('\n');

S.Query        = Search.Query;                      % set up new "Search" data
S.Candidates   = Search.Candidates(j,:);            % re-order candidates
S.Discrepancy  = Search.Discrepancy(j);
S.Disc         = Search.Disc(j,j);
S.DiscComputed = Search.Disc(1,j);
S.AvgDisc      = z;

xListCandidates(File,S,Inf);                 % show on screen
xDisplayCandidates(File,S,1);                % display, level 1

