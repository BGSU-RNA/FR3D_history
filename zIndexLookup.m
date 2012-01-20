% zIndexLookup(File,Num,Chain) finds the base index for base having nucleotide 
% number Num and, if specified, Chain
% Num can be a cell array of N nucleotide numbers
% Any entry of that cell array can use the notation '1830:1835' for a range
% ind is a 1xN vector of the first match to each, the easy answer
% allindices is a 1xN cell array, each row being a 1xN vector
% allchains is a Cx1 cell array as well.

function [ind,allindices,allchains] = zIndexLookup(File,Num,Chain)

if nargin < 3,
  for k = 1:length(Num),
    Chain{k} = '';
  end
end

ind = [];
  
% if File is a text string (filename), load the file and display

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

Numbers = cat(1,{File.NT(:).Number});

for k = 1:length(Num)                      % loop through nucleotide numbers
  if isempty(strfind(Num{k},':')),         % a single number, not a range
    L   = LookUpOne(File,Numbers,Num{k},Chain{k});
    ind = [ind L];                         
    for j = 1:length(L),
      allindices{k,j} = L(j);
      allchains{k,j}  = File.NT(L(j)).Chain;
    end
  else                                     % process a range of numbers
    n = Num{k};                            % kth specified number or range
    i = strfind(Num{k},':');               % find the colon
    Num1 = n(1:(i-1));
    p = LookUpOne(File,Numbers,Num1,Chain{k});
    Num2 = n((i+1):length(n));
    q = LookUpOne(File,Numbers,Num2,Chain{k});
    c = 1;
    for j = 1:length(p),
      for jj = 1:length(q),
        if File.NT(p(j)).Chain == File.NT(q(jj)).Chain,
          ind = [ind p(j):q(jj)];
          allindices{k,c} = [p(j):q(jj)];
          allchains{k,c} = File.NT(p(j)).Chain;
          c = c + 1;
        end
      end
    end
  end
end


function [ind] = LookUpOne(File,Numbers,N,Chain)

    ind = [];
    p = find(ismember(Numbers,N));
    if length(p) == 0,
      fprintf('Could not find nucleotide %s in %s\n',N,File.Filename);
    elseif length(p) == 1 & length(Chain) == 0, % one match, no chain specified
      ind = [ind p];
    elseif length(p) > 1 & length(Chain) == 0,% two matches, no chain specified
      ind = [ind p];
      fprintf('Multiple matches found for %s in %s, consider specifying a chain\n', N, File.Filename);
      for a = 1:length(ind),
        fprintf('Nucleotide %5s Chain %5s Index %5d\n', File.NT(ind(a)).Number, File.NT(ind(a)).Chain, ind(a));
      end
    elseif length(Chain) > 0,                    % chain specified
      c = 0;
      for j = 1:length(p),
        if strcmp(File.NT(p(j)).Chain,Chain),
          ind = [ind p(j)];
          c = c + 1;
        end
      end
      if c == 0,
        fprintf('Could not find nucleotide %s in chain %s in %s\n',N,Chain,File.Filename);
      elseif c > 1,
        fprintf('Multiple matches found for %s in chain %s in %s\n', N,Chain,File.Filename);
      end
    end
