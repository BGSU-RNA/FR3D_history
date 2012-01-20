
for i = current:length(t(:,1)),
  current = i;

  try
    File = zAddNTData(t{i,1},0,[],1);          % load file
  catch
    delete(['PrecomputedData' filesep t{i,1} '.mat']);
    File = zAddNTData(t{i,1},0,[],1);          % load file
  end

  File = zMarkRedundantNucleotides(File,1);       % needed once in May 2010

%  File = zGetPDBInfo(File,n,t);          % get resolution and other info

  if ~isempty(File.NT),                      % if it has nucleotides,
    n(i,2) = length(File.NT);                % store the number of NT

    E  = abs(triu(File.Edge));
    n(i,3) = full(sum(sum((E > 0) .* (E < 13)))); % number of pairs

    LC = File.LongestChain;

    t{i,11} = cat(2,File.NT(LC(1):LC(2)).Base);    % bases in longest chain

    n(i,4) = length(t{i,11});         % number of nucleotides in longest chain
    n(i,5) = LC(1);                   % starting index of longest chain
    n(i,6) = LC(2);                   % end index of longest chain

    if Verbose > 0,
      fprintf('All      %s\n',cat(2,File.NT.Base));
      fprintf('Longest  %s\n',t{i,11});
    end
  end

  zSaveNTData(File);

end
