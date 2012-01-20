% xNeighborhood(File,Indices,v,MaxDiff,MaxInsert) returns indices "near"
% the given indices, in a way determined by the code v

function [NewIndices] = xNeighborhood(File,Indices,v,Display)

N = length(Indices);                                 % number of nucleotides

if nargin < 4,
  
else
  MaxDiff = Display.MaxDiff;
  MaxInsert = Display.MaxInsert;
end

if ~isfield(File,'Distance'),
  c = cat(1,File.NT(1:File.NumNT).Center); % nucleotide centers
  File.Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
end

if isempty(File.Distance),
  c = cat(1,File.NT(1:File.NumNT).Center); % nucleotide centers
  File.Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
end

NewIndices = Indices;                                % always start here

if v > 0,                                          % add basepairs
  for n = 1:N,
    e = abs(File.Edge(Indices(n),:));              % interactions with n
    j = find( (e > 0) .* (e < 14) );               % basepairing with n
    NewIndices = [NewIndices j];
  end
  for n = (N+1):length(NewIndices),                % add what these pair with
    e = abs(File.Edge(NewIndices(n),:));           % interactions with n
    j = find( (e > 0) .* (e < 14) );               % basepairing with n
    NewIndices = [NewIndices j];
  end
end

if v > 1,                                            % add intervening ones
  for n = 1:N,
    NewIndices = [NewIndices (Indices(n)-1):(Indices(n)+1)];
  end
end

if v > 2,
  d = [1 1 8 10 12];
  a = zeros(1,File.NumNT);
  for j=1:length(Indices),
    a = a + (File.Distance(Indices(j),:) < d(v)) .* ...
            (File.Distance(Indices(j),:) > 0);
  end
  NewIndices = [NewIndices find(a)];
end

NewIndices = unique(NewIndices);
NewIndices = sort(NewIndices);
