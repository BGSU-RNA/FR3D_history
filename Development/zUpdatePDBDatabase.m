
% download new .pdb and .pdb1 files - manually at this time

% download a new report, save it as PDB_File_Information current-date.xls
% and modify the following line:

Verbose = 1;

Date = '2009-07-24';

PDBInfoName = ['PDB_File_Information ' Date '.xls'];

% -------------------------------- Create PDBInfo.mat

[n,t,r] = xlsread(PDBInfoName);  % read Excel report on PDB files
t = t(2:end,:);                  % remove the header row
t{1,9} = '';                     % make space to save base sequence
n(1,3) = 0;                      % make space to save number of nucleotides

% save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7

% -------------------------------- Remove duplicate and/or misleading entries

i = find(ismember(t(:,1),'2QOX'));
t{i(1),8} = 'Escherichia coli';                 % change source organism

i = find(ismember(t(:,1),'2QOW'));
t{i(1),8} = 'Escherichia coli';

% Lorena Nasalean says that 2B90 is T.th, not E. coli.

i = find(ismember(t(:,1),'2B9O'));
t{i(1),8} = 'Thermus thermophilus HB8';         % change source organism

i = find(ismember(t(:,1),'1W2B'));
t{i(1),8} = 'HALOARCULA MARISMORTUI';

Keep = zeros(length(t(:,1)));                   % which entries to keep

i = 1;
while i < length(t(:,1))-1,
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

save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7

% ------------------------------------- Read each PDB file now, save data

% load PDBInfo

coun = 0;

for i = 1:length(t(:,1)),
  File = zAddNTData(t{i,1},0,[],1);          % load file
  if ~isempty(File.NT),                      % if it has nucleotides,
    n(i,2) = length(File.NT);                % store the number of NT

    E  = abs(triu(File.Edge));
    n(i,3) = full(sum(sum((E > 0) .* (E < 16)))); % number of pairs

    [F,LC] = zMarkRedundantChains(File,Verbose);  % find redundant chains

    j = find(cat(2,File.NT.Chain) == LC{1});    % indices of longest chain
    t{i,11} = cat(2,File.NT(j).Base);        % bases in largest chain

    if Verbose > 0,
      fprintf('All      %s\n',cat(2,File.NT.Base));
      fprintf('Longest  %s\n',t{i,11});
    end

  end
end

save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7

% save PDBInfo.mat n t -V6      % for compatibility with earlier versions

% ------------------------------------- Find a non-redundant list of PDBs

zFileRedundancy_2

% ------------------------------------- Find new exemplars using NR list

% fid = fopen([PDBFiles filesep 'Nonredundant_' Date '_list.pdb'],'r');

% zFindExemplars