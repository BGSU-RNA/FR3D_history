
% download new .pdb and .pdb1 files - manual at this time

% download a new report, save it as PDB_File_Information current-date.xls
% and modify the following line:

Date = '2008-07-01';

PDBInfoName = ['PDB_File_Information ' Date '.xls'];

% -------------------------------- Create PDBInfo.mat

[n,t,r] = xlsread(PDBInfoName);
t = t(2:end,:);                  % remove the header row
t{1,9} = '';                     % make space to save base sequence
n(1,3) = 0;                      % make space to save number of nucleotides

save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7

% -------------------------------- Remove duplicate, misleading entries

i = find(ismember(t(:,1),'2QOX'));
t{i,8} = 'Escherichia coli';

i = find(ismember(t(:,1),'2QOW'));
t{i(1),8} = 'Escherichia coli';

Keep = zeros(length(t(:,1)));
i = 1;
while i < length(t(:,1))-1,
  if strcmp(lower(t{i,1}),lower(t{i+1,1})),        % duplicate entry
    first = i;
    KeptOne = 0;
    while strcmp(lower(t{i,1}),lower(t{i+1,1})),   % file names are the same
      if ~isempty(strfind(lower(t{i,2}),lower(t{i,8}))),
        Keep(i) = 1;
        KeptOne = 1;
        fprintf('%s must be %s because of %s\n',t{i,1},t{i,8},t{i,2});
      else
        fprintf('%s is not %s because of %s\n',t{i,1},t{i,8},t{i,2});
      end
      i = i + 1;
      if ~isempty(strfind(lower(t{i,2}),lower(t{i,8}))),
        Keep(i) = 1;
        KeptOne = 1;
        fprintf('%s must be %s because of %s\n',t{i,1},t{i,8},t{i,2});
      else
        fprintf('%s is not %s because of %s\n',t{i,1},t{i,8},t{i,2});
      end
    end
    if KeptOne == 0,
      Keep(first) = 1;
    end
  else
    Keep(i) = 1;
  end
  i = i + 1;
end

j = find(Keep);
t = t(j,:);
n = n(j,:);

save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7

% ------------------------------------- Read each PDB file now, save data

for i = 1:length(t(:,1)),
  File = zAddNTData(t{i,1},4,[],1);
  if ~isempty(File.NT),
    n(i,2) = length(File.NT);
    E  = abs(triu(File.Edge));
    n(i,3) = full(sum(sum((E > 0) .* (E < 16)))); % number of pairs
    File = zOrderChains(File);            % largest chain first!
    t{i,9} = cat(2,File.NT.Base);
  end

  save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7
end

% save PDBInfo.mat n t -V6      % for compatibility with earlier versions

% ------------------------------------- Find a non-redundant list of PDBs

zFileRedundancy_1

% ------------------------------------- Find new exemplars using NR list

% fid = fopen([PDBFiles filesep 'Nonredundant_' Date '_list.pdb'],'r');

% zFindExemplars