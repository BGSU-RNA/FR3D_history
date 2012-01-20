
function [Screen] = xPairwiseScreen(File,Codes,Model,p,q);

% --------- Screen according to interaction between nucleotides

D = File.Distance;
%D = tril(File.Distance) + tril(File.Distance)';  % only saved lower tri part

if isfield(Model,'ReqInter'),                 % if screening by interaction

  if length(Model.ReqInter{p,q}) > 0 & all(Model.ReqInter{p,q} ~= 0),
                                           % 0 means no screening for p,q
    E = zeros(size(D));
    for i=1:length(Model.ReqInter{p,q}),
      E = E + (fix(File.Inter) == Model.ReqInter{p,q}(i));
    end
    D = D .* (E > 0);
  end
end
    
[i,j] = find(D);         % nucleotide pairs with OK distances and interactions
d = nonzeros(D);         % distances below the large cutoff

% --------- Screen according to the sequential mask

if Model.Sequential > 0,
  k = find(abs(i-j) <= Model.DiffMat(p,q)); %retain pairs close enough together
  i = i(k);
  j = j(k);
  d = d(k);
end

% --------- Screen according to the nucleotide mask

if (Model.Mask(p) ~= 'N') | (Model.Mask(q) ~= 'N'),
  if (Model.Mask(p) == 'N') & (Model.Mask(q) ~= 'N'),
    k = find(Model.OKCodes{q}(Codes(j)));
  elseif (Model.Mask(p) ~= 'N') & (Model.Mask(q) == 'N'),
    k = find(Model.OKCodes{p}(Codes(i)));
  else
    k = find(Model.OKCodes{p}(Codes(i)) .* Model.OKCodes{q}(Codes(j)));
  end

  i = i(k);
  j = j(k);
  d = d(k);
end

% --------- Screen according to pairwise distance in model

if (Model.Geometric > 0),

  % --------- Compute square of distance difference from model

  d = (d - Model.Distance(p,q)).^2;   % squared difference in distances

  d = d + 0.00000001 * (d == 0);       % avoid rejecting model; make d nonzero

  % --------- Impose upper limit on distance differences; 2-nucleotide cutoff


  if Model.NumNT > 2,
    Wp = Model.LocWeight(p);
    Wq = Model.LocWeight(q);
    MaxD = (Wp + Wq) * (Model.NumNT * Model.DiscCutoff)^2 / (Wp * Wq);
  else
    Wp = 1;
    Wq = 1;
    MaxD = (Model.NumNT * Model.DiscCutoff)^2;
  end

  if isfield(Model,'Flex'),
    MaxD = max(MaxD,Model.Flex(p,q)^2);  % allow larger distance if desired
  end

  k = find(d <= MaxD);            % keep ones with small difference from model

  i = i(k);
  j = j(k);
  d = d(k) * Wp * Wq;

end

% --------- Construct sparse matrix of retained distance difference squares

% [length(i) length(j) length(d) max(k)];

Screen = sparse(i,j,d,File.NumNT,File.NumNT,length(i));



