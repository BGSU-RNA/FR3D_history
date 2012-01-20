% zIndexLookup(File,Num,Chain) finds the base index for base having nucleotide 
% number Num and, if specified, Chain
% Num can be a cell array of N nucleotide numbers
% Any entry of that cell array can use the notation '1830:1835' for a range
% Nucleotide numbers can be followed by (A) to indicate chain A
% ind is a 1xN vector of the first match to each, the easy answer
% allindices is a 1xN cell array, each row being a 1xN vector
% allchains is a Cx1 cell array as well.

function [ind,allindices,allchains] = zIndexLookup(File,Num,Chain)

if nargin < 3,
  for k = 1:length(Num),
    Chain{k} = '';
  end
end

if strcmp(class(Num),'char'),
  Num = {Num};
end

if strcmp(class(Chain),'char'),
  Chain = {Chain};
end

% split multiple specifications into separate cells

t = 1;
for k = 1:length(Num),
  N = regexprep(Num{k},';| ',',');    % replace other delimeters with commas
  while strfind(N,',,'),              % replace double commas with single
    N = regexprep(N,',,',',');
  end
  a = [1 1+strfind(N,',')];               % locations of commas
  for i = 1:(length(a)-1),
    Numb{t} = N(a(i):(a(i+1)-2));
    Chai{t} = Chain{k};
    t = t + 1;
  end
  Numb{t} = N(a(end):end);
  Chai{t} = Chain{k};
  t = t + 1;
end

% check for chain indicated in parentheses

for k = 1:length(Numb),
  if ~isempty(strfind(Numb{k},'(')),
    a = strfind(Numb{k},'(');
    b = strfind(Numb{k},')');
    Chai{k} = Numb{k}(a(1)+1:b(1)-1);     % extract chain
    if length(a) == 1,                    % one chain specified      
      Numb{k} = [Numb{k}(1:a(1)-1) Numb{k}(b(1)+1:end)];
    elseif length(a) == 2,                % range with two chains specified
      Numb{k} = [Numb{k}(1:a(1)-1) Numb{k}(b(1)+1:a(2)-1)];
    end
  end
end

% check for chain indicated by underscore

for k = 1:length(Numb),
  if ~isempty(strfind(Numb{k},'_')),
    a = strfind(Numb{k},'_');
    Chai{k} = Numb{k}(a(end)+1:end);        % extract chain
    Numb{k} = Numb{k}(1:a(end)-1);          % remove chain reference
  end
end

% check for indicated base - what if this happens like A65:G67?

for k = 1:length(Numb),
  Numb{k} = upper(Numb{k});
  aa = strfind(Numb{k},'A');
  ac = strfind(Numb{k},'C');
  ag = strfind(Numb{k},'G');
  au = strfind(Numb{k},'U');
  Base{k} = '';
  if ~isempty(aa),
    Numb{k} = strrep(Numb{k},'A','');
    Base{k} = 'A';
  end
  if ~isempty(ac),
    Numb{k} = strrep(Numb{k},'C','');
    Base{k} = 'C';
  end
  if ~isempty(ag),
    Numb{k} = strrep(Numb{k},'G','');
    Base{k} = 'G';
  end
  if ~isempty(au),
    Numb{k} = strrep(Numb{k},'U','');
    Base{k} = 'U';
  end
end

Numb
Chai
Base


ind = [];
  
% if File is a text string (filename), load the file

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
