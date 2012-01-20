% zUpdateClassificationAndExemplars should be used after zClassLimits has
% been modified.  It redoes the classifications, finds new exemplars, and
% then updates the distance from each pair to the nearest exemplars.

for f=1:length(File),
  File(f) = zClassifyPairs(File(f));
  File(f).Modified = 0;
  zSaveNTData(File(f));
end

zFindExemplars;

for f=1:length(File),
  File(f) = zUpdateDistanceToExemplars(File(f));
  zSaveNTData(File(f));
end

