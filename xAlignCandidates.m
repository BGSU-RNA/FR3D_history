% xAlignCandidates(File,Search,Direction) shows a multiple sequence
% alignment of the candidates in Search.  Bases which correspond in the
% geometric search are aligned with one another.  If a maximum distance
% has been specified, or the maximum gap is small, bases between the
% bases in the candidate are also displayed.

% Direction can be +1 or -1; it tells the order in which to put the 1st cand

function [void] = xAlignCandidates(File,Search,Direction)

Model       = Search.Query;
Candidates  = Search.Candidates;
N           = Model.NumNT;

[L,t] = size(Candidates);

[y,p] = sort(Direction*double(Candidates(1,1:N)));    
                                    % put in increasing or decreasing order
Cand = double(Candidates(:,p));               % re-order nucleotides
F    = Candidates(:,N+1);           % file numbers

if isfield(Model,'MaxDiffMat'),
  MaxDiff = diag(Model.MaxDiffMat(p,p),1);
else
  MaxDiff = Inf*ones(1,N-1);
end

if Direction > 0,
  MaxDiff = MaxDiff;
else
  MaxDiff = fliplr(MaxDiff);
end

% ---------------------------- Calculate maximum gaps between cand. nucleotides

maxinsert = zeros(1,N-1);
for c = 1:L,
  maxinsert = max(maxinsert,abs(diff(double(Cand(c,1:N))))-1);
end

% ---------------------------- Print header line

k = N;
fprintf('               ');
if Model.Geometric > 0,
  fprintf('           ');
end
for j=1:N,
  fprintf('        ');
end
fprintf('    ');
for n = 1:(N-1),
  fprintf('%d',mod(n,10));
  if (MaxDiff(n) < Inf) | (maxinsert(n) < 5),   % if only few insertions
    for i=1:maxinsert(n),
      fprintf(' ');
      k = k + 1;
    end
  else
    fprintf('    ');
  end
end
fprintf('%d\n', mod(N,10));

% ----------------------------- Print alignment

CodeList = zeros(c,k);

for c = 1:L,                                      % loop through candidates
  f = F(c);                                       % file number
  fprintf('%15s', File(f).Filename);
  if Model.Geometric > 0,
    if isfield(Search,'Discrepancy'),
      fprintf('%11.4f',Search.Discrepancy(c));
    elseif isfield(Search,'AvgDisc'),
      fprintf('%11.4f',Search.AvgDisc(c));
    end
  end
  for j=1:N,                            % print candidate
    fprintf('%3s',File(f).NT(Cand(c,j)).Base);    
    fprintf('%5s',File(f).NT(Cand(c,j)).Number);    
  end
  fprintf('    ');

  k = 1;                                        % which letter

  for n = 1:(N-1),                      % print alignment
    fprintf('%s', File(F(c)).NT(Cand(c,n)).Base);
    j = File(F(c)).NT(Cand(c,n)).Code;
    CodeList(c,k) = j;
    k = k + 1;
    if (MaxDiff(n) < Inf) | (maxinsert(n) < 5),   % if only few insertions
      if Cand(c,n+1) - Cand(c,n) > 1,             % increasing order
        for i = (Cand(c,n)+1):(Cand(c,n+1)-1),
          fprintf('%c', File(F(c)).NT(i).Base);   % show insertions
          j = File(F(c)).NT(i).Code;
          CodeList(c,k) = j;
          k = k + 1;
        end
      elseif Cand(c,n+1) - Cand(c,n) < -1,        % decreasing order
        for i = (Cand(c,n)-1):-1:(Cand(c,n+1)+1),
          fprintf('%c', File(F(c)).NT(i).Base);   % show insertions
          j = File(F(c)).NT(i).Code;
          CodeList(c,k) = j;
          k = k + 1;
        end
      end
      for i=1:(1 + maxinsert(n) - abs(Cand(c,n+1)-Cand(c,n))),
        fprintf('-');
        CodeList(c,k) = 5;                         % - is code 5
        k = k + 1;
      end
    else
      fprintf('....');
    end
  end
  fprintf('%s', File(F(c)).NT(Cand(c,N)).Base);
  j = File(F(c)).NT(Cand(c,N)).Code;
  CodeList(c,k) = j;
  fprintf('\n');

  drawnow
end

[m,n] = size(CodeList);
Symbols = {'A', 'C', 'G', 'U', '-'};

fprintf('Letter frequencies\n');
for i = 1:n,
  fprintf('Column %d percentages ',i);
  for j = 1:5,
    fprintf('%s %6.2f  ', Symbols{j}, 100*length(find(CodeList(:,i) == j))/m);
  end
  fprintf('\n');  
end

for i = 1:n,
  for j = 1:n,
    if i < j, 
      fprintf('Pairs made by columns %d and %d\n', i,j);
      Tallies = zeros(5,5);
      for k = 1:m,
        a = CodeList(k,i);
        b = CodeList(k,j);
        Tallies(a,b) = Tallies(a,b) + 1;
      end
      for a = 1:5,
        fprintf('%s ', Symbols{a});
        for b = 1:5,
          fprintf('%8.2f', 100*Tallies(a,b)/m);
        end
        fprintf('\n');
      end
    end
  end
end
