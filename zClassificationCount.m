% zClassificationCount loops through the molecules in File and counts how
% many instances of each pair type and each category occur, displaying the
% output in a table

Classes = [];
PairCodes = [];

for f=1:length(File),
  for i = 1:length(File(f).Pair),
    Classes = [Classes; File(f).Pair(i).Class];
    PairCodes = [PairCodes; File(f).Pair(i).Paircode];
  end
end

[Table, Chi2, P, Labels] = crosstab(Classes,PairCodes);

for i=1:length(Labels(:,2)),
 if ~isempty(Labels{i,2}),
  switch str2num(Labels{i,2}),
    case 1, Labels{i,2} = 'AA';
    case 5, Labels{i,2} = 'AC';
    case 6, Labels{i,2} = 'CC';
    case 7, Labels{i,2} = 'GC';
    case 9, Labels{i,2} = 'AG';
    case 11, Labels{i,2} = 'GG';
    case 13, Labels{i,2} = 'AU';
    case 14, Labels{i,2} = 'CU';
    case 15, Labels{i,2} = 'GU';
    case 16, Labels{i,2} = 'UU';
  end
 end
end

zShowTable(Labels(:,1), Labels(:,2), Table);

fprintf('\nBasepairing interactions only\n');

L = [];

for k=1:length(Labels(:,1)),
  if str2num(Labels{k,1}) < 15,
    L = [L k];
  end
end

zShowTable(Labels(L,1), Labels(:,2), Table(L,:));

fprintf('\nStacking interactions only\n');

L = [];

for k=1:length(Labels(:,1)),
  if (abs(str2num(Labels{k,1})) >= 21) & (abs(str2num(Labels{k,1})) <= 23),
    L = [L k];
  end
end

zShowTable(Labels(L,1), Labels(:,2), Table(L,:));

