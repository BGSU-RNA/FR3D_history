% zIndexLookup(File,Num,Chain) finds the base index for base having nucleotide 
% number Num and, if specified, Chain

function [ind] = zIndexLookup(File,Num,Chain)

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
    ind = [ind LookUpOne(File,Numbers,Num{k},Chain{k})];    
  else                                     % process a range of numbers
    n = Num{k};
    i = strfind(Num{k},':');
    Num1 = n(1:(i-1));
    p = LookUpOne(File,Numbers,Num1,Chain{k});
    Num2 = n((i+1):length(n));
    q = LookUpOne(File,Numbers,Num2,Chain{k});
    if (length(p) > 0) & (length(q) > 0),
      ind = [ind p:q];
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
      fprintf('Multiple matches found for %s in %s, consider specifying a chain\n', N,File.Filename);
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
