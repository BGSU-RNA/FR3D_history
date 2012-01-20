% zUpdateDistanceToExemplars computes the distance between each pair in
% File and the exemplars for each category of interaction.  This indicates
% how close the pair is to the top three categories.

function [File] = zUpdateDistanceToExemplars(File)

if exist('PairExemplars.mat','file') > 0,

  load('PairExemplars','Exemplar');

  for p = 1:length(File.Pair),
    [c,d,h] = zDistanceToExemplars(Exemplar,File.Pair(p));
    File.Pair(p).Classes   = c(1:3);
    File.Pair(p).Distances = d(1:3);
    File.Pair(p).ExemIndex = h(1:3);
  end
else
  for p = 1:length(File.Pair),                  % fill in fictitious distances
    File.Pair(p).Classes   = 99 * ones(1,3);
    File.Pair(p).Distances = 99999999 * ones(1,3);
    File.Pair(p).ExemIndex = [1 2 3];
  end
end
