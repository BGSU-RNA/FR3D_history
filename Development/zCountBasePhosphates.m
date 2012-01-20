% zCountBasePhosphates(File) tabulates the base-phosphate interactions in File
% It formats the counts by geometric family and, within each family, by basepair.
% It writes an Excel sheet with a standard format.
% CountFilename could be 'Basepair_counts_NonRedundant_2008_02_21_list.xls'
% If omitted, the filenames in File are used to construct the filename.

% File = zAddNTData('NonRedundant_2008_02_21_list');
% F    = zResolutionFilter(File,2.5);
% zCountBasePhosphates(F,'BPh_counts');

function [Count,CCount] = zCountBasePhosphates(File,CountFilename)

if nargin == 1,
  CountFilename = 'BasePhosphate_counts';
  for f = 1:length(File),
    CountFilename = [CountFilename '_' File(f).Filename];
  end
  CountFilename = [CountFilename '.xls'];
end

MaxCat = 25;                                  % maximum category number

CI  = [];                                     % count interactions

for f = 1:length(File),
  E = File(f).BasePhosphate;
  for i = 1:length(E(:,1)),
    E(i,i) = 0;                               % eliminate self interactions
  end
  E = E .* (E > 0) .* (E < MaxCat);           % only certain interactions
  [i,j,k] = find(E);
  Code1 = cat(1,File(f).NT(i).Code);
  Code2 = cat(1,File(f).NT(j).Code);
  
  CI  = [CI; [k Code1 Code2]];
end

Count = zeros(4,4,MaxCat);

for i = 1:length(CI(:,1)),
  Count(CI(i,2),CI(i,3),CI(i,1)) = Count(CI(i,2),CI(i,3),CI(i,1)) + 1;
end

BPCat = [2 6 7 0 6 7 8 9 0 1 3 4 5 0 5 9 0 7 4];  % updated 8-19-2008

% ---------------------------------------- start producing output

row = 1;
clear T
Letters = 'ACGU';

% ---------------------------------------- counts by interacting base

NC = [];

for i = 1:4,                               % base making the base-phosphate
  for j = 0:9,                             % interaction being made
    NC(i,j+1) = sum(sum(Count(i,:,find(BPCat == j))));
  end
end

T{row,1} = 'Base nucleotide';
T{row,12} = 'Total 1BP to 9BP';
T{row,13} = 'Total';

for j = 0:9,
  T{row,j+2} = [num2str(j) 'BP'];        % column labels
end

for i = 1:4,
  T{row+i,1} = Letters(i);
  for j = 0:9,
    T{row+i,2+j} = NC(i,j+1);
  end
  T{row+i,12} = sum(NC(i,2:10));
  T{row+i,13} = sum(NC(i,:));
end

% ---------------------------------------- counts by phosphate identity

NC = [];

row = row + 6;

for i = 1:4,                               % base making the base-phosphate
  for j = 0:9,                             % interaction being made
    NC(i,j+1) = sum(sum(Count(:,i,find(BPCat == j))));
  end
end

T{row,1} = 'Phosphate nucleotide';
T{row,12} = 'Total 1BP to 9BP';
T{row,13} = 'Total';

for j = 0:9,
  T{row,j+2} = [num2str(j) 'BP'];        % column labels
end

for i = 1:4,
  T{row+i,1} = Letters(i);
  for j = 0:9,
    T{row+i,2+j} = NC(i,j+1);
  end
  T{row+i,12} = sum(NC(i,2:10));
  T{row+i,13} = sum(NC(i,:));
end

% ---------------------------------------- counts by pair of bases

row = row + 6;
NC = [];

for i = 1:4,                               % base making the base-phosphate
 for k = 1:4,
  for j = 0:9,                             % interaction being made
    NC(i,k,j+1) = sum(Count(i,k,find(BPCat == j)));
  end
 end
end

T{row,1} = 'Pair of bases';
T{row,12} = 'Total 1BP to 9BP';
T{row,13} = 'Total';

for j = 0:9,
  T{row,j+2} = [num2str(j) 'BP'];        % column labels
end

for i = 1:4,
 for k = 1:4,
  p = 4*(i-1) + k;
  T{row+p,1} = Letters([i k]);
  for j = 0:9,
    T{row+p,2+j} = NC(i,k,j+1);
  end
  T{row+p,12} = sum(NC(i,k,2:10));
  T{row+p,13} = sum(NC(i,k,:));
 end
end

% --------------------------------------- write to spreadsheet

if ~isempty(CountFilename),
  xlswrite(CountFilename,T);
end


