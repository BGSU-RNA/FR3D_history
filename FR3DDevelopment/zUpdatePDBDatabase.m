
% download new .pdb and .pdb1 files - manually at this time

Verbose = 1;

% download a new report, save it as PDB_File_Information current-date.xls
% and modify the following line:

Date = '2010-05-19';

PDBInfoName = ['PDB_File_Information ' Date '.xls'];

% Columns should correspond to this list:
% 1 A	Structure ID
% 2 B	Descriptor / Structure title
% 3 C	Experimental Technique / Exp. method
% 4 D	Release Date
% 5 E	Authors
% 6 F	Keywords
% 7 G	Resolution (Å)
% 8 H	Source

% We should add an additional column for NDB ID

% Later in this program, other columns are set in t and in n:

% 11    Sequence of longest chain

% -------------------------------- Create PDBInfo.mat

[n,t] = xlsread(PDBInfoName);  % read Excel report on PDB files

t = t(2:end,:);                  % remove the header row
t{1,11} = '';                     % make space to save base sequence
n(1,5) = 0;                      % make space to save numeric values

for i = 1:length(t(:,7)),
  a = str2num(t{i,7});
  if isempty(a),
    n(i,1) = NaN;
  else
    n(i,1) = a;
  end
end

% save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7

% -------------------------------- Fix a few entries

if 0 > 1,
  i = find(ismember(t(:,1),'2QOX'));
  t{i(1),8} = 'Escherichia coli';                 % change source organism

  i = find(ismember(t(:,1),'2QOW'));
  t{i(1),8} = 'Escherichia coli';

  % Lorena Nasalean says that 2B90 is T.th, not E. coli.
 
  i = find(ismember(t(:,1),'2B9O'));
  t{i(1),8} = 'Thermus thermophilus HB8';         % change source organism

  i = find(ismember(t(:,1),'1W2B'));
  t{i(1),8} = 'HALOARCULA MARISMORTUI';
end


if 0 > 1,
% --------------------------------- Remove duplicate entries - before 2010
Keep = zeros(length(t(:,1)));                   % which entries to keep

i = 1;
while i < length(t(:,1))-1,                     % go through entries
  if strcmp(lower(t{i,1}),lower(t{i+1,1})),     % i, i+1 are duplicate entries
    first = i;
    KeptOne = 0;
    while strcmp(lower(t{i,1}),lower(t{i+1,1})),   % file names are the same
      if ~isempty(strfind(lower(t{i,2}),lower(t{i,8}))), % source matches desc
        Keep(i) = 1;
        KeptOne = 1;
        fprintf('%s must be %s because of %s\n',t{i,1},t{i,8},t{i,2});
      end
      i = i + 1;
    end
    
    if KeptOne == 0,
      Keep(first) = 1;
      fprintf('Can''t be sure which one of these annotations to keep:\n');
      for j = first:i,
        fprintf('%4s | %25s | %s\n', t{j,1}, t{j,8}, t{j,2});
        if j > first,
          t{first,8} = [t{first,8} ' | ' t{j,8}];  % append source organisms
        end
      end
    end

    fprintf('\n');
  else
    Keep(i) = 1;
  end
  i = i + 1;
end

j = find(Keep);
t = t(j,:);
n = n(j,:);

end

% --------------------------------- Remove duplicate entries - May 2010
Keep = zeros(length(t(:,1)),1);                   % which entries to keep

Keep(1) = 1;
first = 1;
i = 2;

while i < length(t(:,1)),
  while (i < length(t(:,1))) && strcmp(lower(t{first,1}),lower(t{i,1})),
                                            % first, i are duplicate entries
    if ~strcmp(t{first,2},t{i,2}),
      fprintf('Different descriptors for %s\n', t{first,1});
      fprintf('%s\n', t{first,2});
      fprintf('%s\n', t{i,2});
    end

    i = i + 1;

  end

%  fprintf('Kept %s\n', t{i,1});

  Keep(i) = 1;
  first = i;
  i = i + 1;

end

j = find(Keep);
t = t(j,:);
n = n(j,:);

save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7

% ------------------------------------- Read each PDB file now, save data

% load PDBInfo

current = 1;

zUpdatePDBDatabaseLoop

save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7

% save PDBInfo.mat n t -V6      % for compatibility with earlier versions

% ------------------------------------- Find a non-redundant list of PDBs

zFileRedundancy_2

% ------------------------------------- Find new exemplars using NR list

% fid = fopen([PDBFiles filesep 'Nonredundant_' Date '_list.pdb'],'r');

% zFindExemplars