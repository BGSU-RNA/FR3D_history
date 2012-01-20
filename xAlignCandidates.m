% Direction can be +1 or -1; it tells the order in which to put the 1st cand

function [void] = xAlignCandidates(File,Model,Candidates,Discrepancy,Direction)

L = length(Discrepancy);

[y,p] = sort(Direction*double(Candidates(1,1:Model.NumNT)));    % put in increasing or decreasing order

if ~isfield(Model,'MaxDiff'),
  Model.MaxDiff = Inf*ones(1,Model.NumNT-1);
end

if Direction > 0,
  MaxDiff = Model.MaxDiff;
else
  MaxDiff = fliplr(Model.MaxDiff);
end

Cand = double(Candidates(:,p));               % re-order nucleotides
F    = Candidates(:,Model.NumNT+1);           % file numbers

maxinsert = zeros(1,Model.NumNT-1);

for c = 1:L,
  maxinsert = max(maxinsert,abs(diff(Cand(c,1:Model.NumNT)))-1);
end

for c = 1:L,
  f = F(c);
  fprintf('%15s', File(f).Filename);
  if Model.Geometric > 0,
    fprintf('%11.4f',Discrepancy(c));
  end
  for j=1:Model.NumNT,
    fprintf('%3s',File(f).NT(Cand(c,j)).Base);    
    fprintf('%4s',File(f).NT(Cand(c,j)).Number);    
  end
  fprintf('    ');
  for n = 1:(Model.NumNT-1),
    fprintf('%s', File(F(c)).NT(Cand(c,n)).Base);
    if MaxDiff(n) < Inf,
      if Cand(c,n+1) - Cand(c,n) > 1,
        for i = (Cand(c,n)+1):(Cand(c,n+1)-1),
          fprintf('%s', File(F(c)).NT(i).Base);
        end
      elseif Cand(c,n+1) - Cand(c,n) < -1,
        for i = (Cand(c,n)-1):-1:(Cand(c,n+1)+1),
          fprintf('%s', File(F(c)).NT(i).Base);
        end
      end
      for i=1:(1 + maxinsert(n) - abs(Cand(c,n+1)-Cand(c,n))),
        fprintf('-');
      end
    else
      fprintf('....');
    end
  end
  fprintf('%s', File(F(c)).NT(Cand(c,Model.NumNT)).Base);
  fprintf('\n');
  drawnow
end

