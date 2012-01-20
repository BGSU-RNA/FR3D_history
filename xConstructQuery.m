% xConstructModel(Model,File,f) fills in details of Model from File

% defaults:
% Model.Geometric = 1 if Model.Filename and Model.NTList are defined
% Model.Geometric = 0 if not; in that case, Model.Mask, Model.ReqInter, or
%      Model.MaxDiff need to be defined, or there is no focus for the search
% Model.ExcludeOverlap = 1 if there are 7 or more nucleotides
% Model.Mask is "NNNNN..."
% Model.ReqInter is empty
% Model.Sequential is 1 if Model.MaxDiff is non-empty
% Model.MaxDiff is [Inf Inf ...] by default
% Model.LocWeight is [1 1 1 1...]
% Model.AngleWeight is [1 1 1 1...]
% Model.DiscCutoff is 0.4
% Model.RelCutoff is Model.DiscCutoff, if not overridden
% Model.ChainList is not needed unless there is ambiguity in nucleotide
%                 numbers, and if there is, xConstructModel will tell you

function [Model] = xConstructModel(Model,File,f)

if ~isfield(Model,'Geometric'),
  if isfield(Model,'Filename') & isfield(Model,'NTList'),
    Model.Geometric = 1;
  else
    Model.Geometric = 0;
  end
end

if Model.Geometric == 1,
  Model.NumNT     = length(Model.NTList);
end  

if Model.Geometric == 0,
  if isfield(Model,'ReqInter'),
    [s,t] = size(Model.ReqInter);
    Model.NumNT = max(s,t);               % number of nucleotides
  elseif isfield(Model,'Mask'),
    Model.NumNT = length(Model.Mask);
  elseif isfield(Model,'MaxDiff'),
    Model.NumNT = length(Model.MaxDiff)+1;
  else
    fprintf('Not enough information to form a query motif.\n');
    fprintf('No search was conducted.\n');
  end
end

if isfield(Model,'NumNT'),

% --------- Set default values

if ~isfield(Model,'ExcludeOverlap'),
  if Model.NumNT >= 7,
    Model.ExcludeOverlap = 1;
  else
    Model.ExcludeOverlap = 0;
  end
end

if ~isfield(Model,'MaxDiff'),
  Model.Sequential = 0;
end

if ~isfield(Model,'Sequential') & isfield(Model,'MaxDiff'),
  Model.Sequential = 1;
end

if ~isfield(Model,'Mask'),
  Model.Mask = char('N' * ones(1,Model.NumNT));
end

if ~isfield(Model,'Name'),
  Model.Name = num2str(Model.Number);
end

% ---------------------------------------------------------------------

if Model.Geometric == 1,                    % model comes from a file
    
  if ~isfield(Model,'LocWeight'),
    Model.LocWeight = ones(1,Model.NumNT);
  end

  if ~isfield(Model,'AngleWeight'),
    Model.AngleWeight = ones(1,Model.NumNT);
  end

  if ~isfield(Model,'DiscCutoff'),
    Model.DiscCutoff = 0.4;
  end

  if ~isfield(Model,'RelCutoff'),
    Model.RelCutoff = Model.DiscCutoff;
  end

  % --------- basic parameters for the model

  if isfield(Model,'ChainList'),
    Model.Indices  = zIndexLookup(File,Model.NTList,Model.ChainList);
  else
    Model.Indices  = zIndexLookup(File,Model.NTList);
  end

  for i=1:Model.NumNT,
    Model.NT(i) = File.NT(Model.Indices(i));
  end

  Model.Inter = File.Inter(Model.Indices,Model.Indices);

  Model.Centers = cat(1,Model.NT.Center);
  Model.Distance = zMutualDistance(Model.Centers,Inf);
end

% --------- Make interaction matrix symmetric

if isfield(Model,'ReqInter'),
  Model.ReqInter{Model.NumNT,Model.NumNT} = 0;
  for i=1:Model.NumNT,
    for j=(i+1):Model.NumNT,
      if ~isempty(Model.ReqInter{j,i}),
        Model.ReqInter{i,j} = Model.ReqInter{j,i};
      elseif ~isempty(Model.ReqInter{i,j}),
        Model.ReqInter{j,i} = Model.ReqInter{i,j};
      else
        Model.ReqInter{i,j} = 0;
        Model.ReqInter{j,i} = 0;
      end
      if length(Model.ReqInter{i,j}) > 0,
        Model.ReqInter{i,j} = Model.ReqInter{i,j}(find(Model.ReqInter{i,j}));
        Model.ReqInter{j,i} = Model.ReqInter{i,j};
      end
    end
  end
end

% --------- precompute parameters for geometric screening and ranking

if (Model.Geometric > 0) & (Model.NumNT >= 2),

  Model.DistCutoff = max(max(Model.Distance)) + sqrt(2)*Model.NumNT*Model.DiscCutoff;
                                     % distances this large needed in File

  Model.LocWeight = Model.NumNT * Model.LocWeight / sum(Model.LocWeight);
                                        % weights sum to Model.NumNT

  Model.SSCutoff  = (Model.NumNT^2)*(Model.RelCutoff^2)*cumsum(Model.LocWeight);
                                    % cutoffs for first 1, 2, 3, ... nucleotides

  Model.WeightedCenter = Model.LocWeight * Model.Centers / Model.NumNT;

  Model.CenteredCenters = Model.Centers-ones(Model.NumNT,1)*Model.WeightedCenter;
  Model.WeightedCenteredCenters = diag(Model.LocWeight)* Model.CenteredCenters;

  Model.LDiscCutoff = (Model.NumNT*Model.RelCutoff)^2;

  if isfield(Model,'Flex'),
    Model.DistanceScreen = 0;            % cannot use sums of squares with flex
    Model.Flex(Model.NumNT,Model.NumNT) = 0; % make big enough
    Model.Flex = Model.Flex + Model.Flex';   % make symmetric
  end
else
  Model.SSCutoff = Inf * ones(1,Model.NumNT);
end

% --------- Special calculations for two-nucleotide motifs

if Model.NumNT == 2 & Model.Geometric > 0,
  Model.DistCutoff = Model.Distance(1,2) + 2 * Model.DiscCutoff;
  Model.LocWeight  = ones(1,Model.NumNT);
  Model.R  = Model.NT(2).Rot' * Model.NT(1).Rot;
  Model.T1 = (Model.NT(2).Center - Model.NT(1).Center)*Model.NT(1).Rot;
  Model.T2 = (Model.NT(1).Center - Model.NT(2).Center)*Model.NT(2).Rot;
end


% --------- implement sequential screening regardless of specified order

if Model.Sequential > 0,
  D = zeros(Model.NumNT,Model.NumNT);

  for i=1:(Model.NumNT-1),
    D(i,i+1) = Model.MaxDiff(i);
  end

  for d=3:Model.NumNT,                     % diagonal to consider
    for i=1:(Model.NumNT-d+1),
      j = i + d - 1;
      D(i,j) = D(i,j-1) + D(i+1,j);
    end
  end

  Model.DiffMat = D + D';

  if isfield(Model,'Filename'),
    [y,i] = sort(Model.Indices);
  else
    i = 1:Model.NumNT;
  end

  k(i) = 1:Model.NumNT;

  Model.DiffMat = Model.DiffMat(k,k);
end

% --------- Interpret model mask

% For nucleotide i, the vector Model.OKCodes{i} has 0's and 1's 
% which tell which nucleotide codes (A=1, C=2,...) meet the mask

for i=1:Model.NumNT,
  switch Model.Mask(i)
   case 'A', OK = [1 0 0 0];
   case 'C', OK = [0 1 0 0];
   case 'G', OK = [0 0 1 0];
   case 'U', OK = [0 0 0 1];
   case 'M', OK = [1 1 0 0];
   case 'R', OK = [1 0 1 0];
   case 'W', OK = [1 0 0 1];
   case 'S', OK = [0 1 1 0];
   case 'Y', OK = [0 1 0 1];
   case 'K', OK = [0 0 1 1];
   case 'V', OK = [1 1 1 0];
   case 'H', OK = [1 1 0 1];
   case 'D', OK = [1 0 1 1];
   case 'B', OK = [0 1 1 1];
   case 'N', OK = [1 1 1 1];
  end
  Model.OKCodes{i} = OK;
end

end
