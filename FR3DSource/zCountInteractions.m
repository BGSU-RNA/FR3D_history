% zCountInteractions(File) tabulates the number of each type of interaction in File

CountFilename = 'CountsFromNonRedundantList.xls';

MaxCat = 25;

CI = [];

for f = 1:length(File),
  E = File(f).Edge;
  E = E .* (E > 0) .* (E < MaxCat);           % only certain interactions
  [i,j,k] = nonzeros(E);
  Code1 = cat(1,File(f).NT(j).Code);
  Code2 = cat(1,File(f).NT(k).Code);
  
  CI = [CI; [i Code1 Code2]];
end

CI = [floor(CI(:,1)*10) CI(:,2:3)];           % multiply classes by 10

Count = zeros(4,4,MaxCat);

for i = 1:length(CI(:,1)),
  Count(CI(i,2),CI(i,3),CI(i,1)) = Count(CI(i,2),CI(i,3),CI(i,1)) + 1;
end

zDisplayPairCounts(Count,CountFilename);