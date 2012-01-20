% xPairwiseScreen returns a sparse matrix with non-zero entries corresponding to pairs of bases which satisfy all given constraints
% Codes are the base codes from File
% 

function [Screen] = xPairwiseScreen(File,Codes,Query,p,q,PC);

if Query.Geometric == 0,
  D = File.Distance .* (File.Distance < Query.Diameter);
                                        % cap distance for non-geometric search
else
  D = File.Distance;
end

% --------- Screen according to interaction between nucleotides

if isfield(Query,'EdgeNums'),           % if screening by edges, incl stacking
  if length(Query.EdgeNums{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.EdgeNums{p,q}),
      E = E + (fix(File.Edge) == Query.EdgeNums{p,q}(i));
    end
    D = D .* (E > 0);                         % include only those that match
  end
end
    
if isfield(Query,'ExcludeEdges'),                 % if excluding by edges
  if length(Query.ExcludeEdges{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.ExcludeEdges{p,q}),
      E = E + (fix(File.Edge) == Query.ExcludeEdges{p,q}(i));
    end
    D = D .* (E == 0);
  end
end
    
if isfield(Query,'OKPairs'),                 % if screening by paircode
  if length(Query.OKPairs{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.OKPairs{p,q}),
      E = E + sparse(PC == Query.OKPairs{p,q}(i));
    end
    D = D .* (E > 0);
  end
end
    
if isfield(Query,'ExPairs'),                 % if screening by paircode
  if length(Query.ExPairs{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.ExPairs{p,q}),
      E = E + sparse(PC == Query.ExPairs{p,q}(i));
    end
    D = D .* (E == 0);
  end
end
    
if isfield(Query,'BasePhos'),                 % if screening by base-phosphate
  if length(Query.BasePhos{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.BasePhos{p,q}),
      E = E + (fix(File.BasePhosphate) == Query.BasePhos{p,q}(i));
    end
    D = D .* (E > 0);
  end
  if length(Query.BasePhos{q,p} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.BasePhos{q,p}),
      E = E + (fix(File.BasePhosphate') == Query.BasePhos{q,p}(i));
    end
    D = D .* (E > 0);
  end
end
    
if isfield(Query,'ExcludeBasePhos'),                 % if excluding by edges
  if length(Query.ExcludeBasePhos{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.ExcludeBasePhos{p,q}),
      E = E + (fix(File.BasePhosphate) == Query.ExcludeBasePhos{p,q}(i));
    end
    D = D .* (E == 0);
  end
  if length(Query.ExcludeBasePhos{q,p} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.ExcludeBasePhos{q,p}),
      E = E + (fix(File.BasePhosphate') == Query.ExcludeBasePhos{q,p}(i));
    end
    D = D .* (E == 0);
  end
end

if isfield(Query,'Flank'),
  if ~isempty(Query.Flank{p,q}),
    E = triu(fix(abs(File.Edge))==1) .* (File.Range == 0); % nested cWW's
    H = sparse(zeros(File.NumNT,File.NumNT));
    [i,j] = find(E);                     % indices of NT's making nested cWW's
    a = [i; j];                           % all indices of nested cWW pairs
    a = sort(a);
    for k = 1:(length(a)-1),
      if a(k+1) - a(k) > 1,
        H(a(k),a(k+1)) = 1;              % these two flank something
      end
    end
    H = H + H';
    D = D .* H;                          % only keep pairs between nested cWW

  end
end

if isfield(Query,'Range'),               % screen by range restriction
  if ~isempty(Query.Range{p,q}),
    R = Query.Range{p,q};
    D = D .* (File.Range >= R(1)) .* (File.Range <= R(2));
  end
end

if isfield(Query,'Coplanar'),            % screen by coplanarity
 if ~isempty(Query.Coplanar{p,q}),
  y = 0.4;
  [i,j] = find(D);                       % all remaining pairs near each other
  CP = sparse(zeros(File.NumNT,File.NumNT));

  for k = 1:length(i);
    N1 = File.NT(i(k));
    N2 = File.NT(j(k));
    mind  = min(min(zDistance(N1.Fit,N2.Fit)));
    if mind < 5,                         % only keep those within 5 Angstroms
      v     = N1.Center - N2.Center;
      v     = v / norm(v);                  % normalize

      if (abs(v * N1.Rot(:,3)) < y) && (abs(v * N2.Rot(:,3)) < y), % angle > 45 degrees

        Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

        d = zDistance(N2.Fit(1:Lim(2,N2.Code),:), N1.Center); 
                                           % distances to base 1 center
        [y,m] = min(d);                          % identify the closest atom
        m = m(1);                                % in case of a tie, use the first
        Gap1 = N1.Rot(:,3)'*(N2.Fit(m,:)-N1.Center)';% height above plane of 1

        if abs(Gap1) < 2,

          d = zDistance(N1.Fit(1:Lim(2,N1.Code),:), N2.Center); 
                                           % distances to base 1 center
          [y,m] = min(d);                          % identify the closest atom
          m = m(1);                                % in case of a tie, use the first
          Gap2 = N2.Rot(:,3)'*(N1.Fit(m,:)-N2.Center)';% height above plane of 1

          if abs(Gap2) < 2,

            CP(i(k),j(k)) = 1;
            CP(j(k),i(k)) = 1;
          end
        end
      end
    end
  end

  if Query.Coplanar{p,q} > 0,
    D = D .* CP;
  elseif Query.Coplanar{p,q} == 0,
    D = D .* (CP == 0);
  end
 end
end
    
[i,j] = find(D);         % nucleotide pairs with OK distances and interactions
d = nonzeros(D);         % distances below the large cutoff

% --------- Screen according to maximum difference in nucleotide numbers

if isfield(Query,'MaxDiffMat'),
  if Query.MaxDiffMat(p,q) < Inf,
    k = find(abs(i-j) <= Query.MaxDiffMat(p,q)); %retain pairs close enough together
    i = i(k);
    j = j(k);
    d = d(k);
  end
end

% --------- Screen according to minimum difference in nucleotide numbers

if isfield(Query,'MinDiffMat'),
  if Query.MinDiffMat(p,q) > 1,
    k = find(abs(i-j) >= Query.MinDiffMat(p,q)); %retain pairs far enough apart
    i = i(k);
    j = j(k);
    d = d(k);
  end
end

% --------- Screen according to sign of nucleotide number difference

if isfield(Query,'DifferenceSignMat'),
  if Query.DifferenceSignMat(p,q) < 0,
    k = find(j > i);                          % retain pairs with num2>num1
    i = i(k);
    j = j(k);
    d = d(k);
  elseif Query.DifferenceSignMat(p,q) > 0,
    k = find(j < i);                          % retain pairs with num2<num1
    i = i(k);
    j = j(k);
    d = d(k);
  end
end

% --------- Screen according to the nucleotide mask

if (min(Query.OKCodes{p}) == 0) | (min(Query.OKCodes{q}) == 0),
  if (min(Query.OKCodes{p}) == 1) & (min(Query.OKCodes{q}) == 0),
    k = find(Query.OKCodes{q}(Codes(j)));
  elseif (min(Query.OKCodes{p}) == 0) & (min(Query.OKCodes{q}) == 1),
    k = find(Query.OKCodes{p}(Codes(i)));
  else
    k = find(Query.OKCodes{p}(Codes(i)) .* Query.OKCodes{q}(Codes(j)));
  end

  i = i(k); 
  j = j(k);
  d = d(k);
end

% --------- Screen according to configuration (syn or anti)

if length(Query.Config{p}) > 0,
  switch Query.Config{p}
    case 'syn'
      k = find(cat(1,File.NT(i).Syn) == 1);
    case 'anti'
      k = find(cat(1,File.NT(i).Syn) == 0);
    otherwise
      k = 1:length(i);
  end

  i = i(k);
  j = j(k);
  d = d(k);
end

if length(Query.Config{q}) > 0,
  switch Query.Config{q}
    case 'syn'
      k = find(cat(1,File.NT(j).Syn) == 1);
    case 'anti'
      k = find(cat(1,File.NT(j).Syn) == 0);
    otherwise
      k = 1:length(j);
  end

  i = i(k);
  j = j(k);
  d = d(k);
end

% --------- Screen according to pairwise distance in model

if (Query.Geometric > 0),

  % --------- Compute square of distance difference from model

  d = (d - Query.Distance(p,q)).^2;   % squared difference in distances

  d = d + 0.00000001 * (d == 0);       % avoid rejecting model; make d nonzero

  % --------- Impose upper limit on distance differences; 2-nucleotide cutoff


  if Query.NumNT > 2,
    Wp = Query.LocWeight(p);
    Wq = Query.LocWeight(q);
    MaxD = (Wp + Wq) * (Query.NumNT * Query.DiscCutoff)^2 / (Wp * Wq);
  else
    Wp = 1;
    Wq = 1;
    MaxD = (Query.NumNT * Query.DiscCutoff)^2;
  end

  if isfield(Query,'Flex'),
    MaxD = max(MaxD,Query.Flex(p,q)^2);  % allow larger distance if desired
  end

  k = find(d <= MaxD);            % keep ones with small difference from model

  i = i(k);
  j = j(k);
  d = d(k) * Wp * Wq;

end

% --------- Construct sparse matrix of retained distance difference squares

% [length(i) length(j) length(d) max(k)];

Screen = sparse(i,j,d,File.NumNT,File.NumNT,length(i));



