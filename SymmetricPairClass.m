
  To get started, select MATLAB Help or Demos from the Help menu.

>> type pdist

function Y = pdist(X,dist,varargin)
%PDIST Pairwise distance between observations.
%   Y = PDIST(X) returns a vector Y containing the Euclidean distances
%   between each pair of observations in the N-by-P data matrix X.  Rows of
%   X correspond to observations, columns correspond to variables.  Y is a
%   1-by-(N*(N-1)/2) row vector, corresponding to the N*(N-1)/2 pairs of
%   observations in X.
%
%   Y = PDIST(X, DISTANCE) computes Y using DISTANCE.  Choices are:
%
%       'euclidean'   - Euclidean distance
%       'seuclidean'  - Standardized Euclidean distance, each coordinate
%                       in the sum of squares is inverse weighted by the
%                       sample variance of that coordinate
%       'cityblock'   - City Block distance
%       'mahalanobis' - Mahalanobis distance
%       'minkowski'   - Minkowski distance with exponent 2
%       'cosine'      - One minus the cosine of the included angle
%                       between observations (treated as vectors)
%       'correlation' - One minus the sample linear correlation between
%                       observations (treated as sequences of values).
%       'spearman'    - One minus the sample Spearman's rank correlation
%                       between observations (treated as sequences of values).
%       'hamming'     - Hamming distance, percentage of coordinates
%                       that differ
%       'jaccard'     - One minus the Jaccard coefficient, the
%                       percentage of nonzero coordinates that differ
%       'chebychev'   - Chebychev distance (maximum coordinate difference)
%       function      - A distance function specified using @, for
%                       example @DISTFUN
%
%   A distance function must be of the form
%
%         function D = DISTFUN(XI, XJ, P1, P2, ...),
%
%   taking as arguments a 1-by-P vector XI containing a single row of X, an
%   M-by-P matrix XJ containing multiple rows of X, and zero or more
%   additional problem-dependent arguments P1, P2, ..., and returning an
%   M-by-1 vector of distances D, whose Kth element is the distance between
%   the observations XI and XJ(K,:).
%
%   Y = PDIST(X, DISTFUN, P1, P2, ...) passes the arguments P1, P2, ...
%   directly to the function DISTFUN.
%
%   Y = PDIST(X, 'minkowski', P) computes Minkowski distance using the
%   positive scalar exponent P.
%
%   The output Y is arranged in the order of ((1,2),(1,3),..., (1,N),
%   (2,3),...(2,N),.....(N-1,N)), i.e. the lower left triangle of the full
%   N-by-N distance matrix in column order.  To get the distance between
%   the Ith and Jth observations (I < J), either use the formula
%   Y((I-1)*(N-I/2)+J-I), or use the helper function Z = SQUAREFORM(Y),
%   which returns an N-by-N square symmetric matrix, with the (I,J) entry
%   equal to distance between observation I and observation J.
%
%   Example:
%
%      X = randn(100, 5);                 % some random points
%      Y = pdist(X, 'euclidean');         % unweighted distance
%      Wgts = [.1 .3 .3 .2 .1];           % coordinate weights
%      Ywgt = pdist(X, @weucldist, Wgts); % weighted distance
%
%      function d = weucldist(XI, XJ, W) % weighted euclidean distance
%      d = sqrt((repmat(XI,size(XJ,1),1)-XJ).^2 * W');
%
%   See also SQUAREFORM, LINKAGE, SILHOUETTE.

%   An example of distance for data with missing elements:
%
%      X = randn(100, 5);     % some random points
%      X(unidrnd(prod(size(X)),1,20)) = NaN; % scatter in some NaNs
%      D = pdist(X, @naneucdist);
%
%      function d = naneucdist(XI, XJ) % euclidean distance, ignoring NaNs
%      [m,p] = size(XJ);
%      sqdx = (repmat(XI,m,1) - XJ) .^ 2;
%      pstar = sum(~isnan(sqdx),2); % correction for missing coords
%      pstar(pstar == 0) = NaN;
%      d = sqrt(nansum(sqdx,2) .* p ./ pstar);
%
%
%   For a large number of observations, it is sometimes faster to compute
%   the distances by looping over coordinates of the data (though the code
%   is more complicated):
%
%      function d = nanhamdist(XI, XJ) % hamming distance, ignoring NaNs
%      [m,p] = size(XJ);
%      nesum = zeros(m,1);
%      pstar = zeros(m,1);
%      for q = 1:p
%          notnan = ~(isnan((XI(q)) | isnan(XJ(:,q)));
%          nesum = nesum + (XI(q) ~= XJ(:,q)) & notnan;
%          pstar = pstar + notnan;
%      end
%      nesum(any() | nans((i+1):n)) = NaN;
%      Y(k:(k+n-i-1)) = nesum ./ pstar;

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.15.4.10 $ $Date: 2004/08/20 20:06:05 $

if nargin < 2
    dist = 'euc';
else
    if ischar(dist)
        methods = {'euclidean'; 'seuclidean'; 'cityblock'; 'chebychev'; ...
                   'mahalanobis'; 'minkowski'; 'cosine'; 'correlation'; ...
                   'spearman'; 'hamming'; 'jaccard'};
        i = strmatch(lower(dist), methods);
        if length(i) > 1
            error('stats:pdist:BadDistance',...
                  'Ambiguous ''DISTANCE'' argument:  %s.', dist);
        elseif isempty(i)
            % Assume an unrecognized string is a user-supplied distance
            % function name, change it to a handle.
            distfun = str2func(dist);
            distargs = varargin;
            dist = 'usr';
        else
            dist = lower(methods{i}(1:3));
        end
    elseif isa(dist, 'function_handle') ||  isa(dist, 'inline')
        distfun = dist;
        distargs = varargin;
        dist = 'usr';
    else
        error('stats:pdist:BadDistance',...
              'The ''DISTANCE'' argument must be a string or a function.');
    end
end

% Integer/logical/char/anything data may be handled by a caller-defined
% distance function, otherwise it is converted to double.  Complex floating
% point data must also be handled by a caller-defined distance function.
if ~strcmp(dist,'usr')
    if ~isfloat(X)
        warning('stats:pdist:DataConversion', ...
                'Converting %s data to double.',class(X));
        X = double(X);
    elseif any(imag(X(:)))
        error('stats:pdist:InvalidData', ...
              'PDIST does not accept complex data for built-in distances.');
    end
end

[n,p] = size(X);

% Degenerate case, just return an empty of the proper size.
if n < 2
    if ~strcmp(dist,'usr')
        Y = zeros(1,0,class(X)); % X was single/double, or cast to double
    else
        Y = zeros(1,0);
    end
    return;
end

switch dist
case 'seu' % Standardized Euclidean weights by coordinate variance
    additionalArg = 1 ./ var(X)';
case 'mah' % Mahalanobis
    additionalArg = cov(X) \ eye(p); %inv(cov(X));
case 'min' % Minkowski distance needs a third argument
    if nargin < 3  % use default value for exponent
        additionalArg = 2;
    elseif varargin{1} > 0
        additionalArg = varargin{1}; % get exponent from input args
    else
        error('stats:pdist:InvalidExponent',...
              'The exponent for the Minkowski metric must be positive.');
    end
case 'cos' % Cosine
    Xnorm = sqrt(sum(X.^2, 2));
    if min(Xnorm) <= eps(max(Xnorm))
        error('stats:pdist:InappropriateDistance',...
              ['Some points have small relative magnitudes, making them ', ...
               'effectively zero.\nEither remove those points, or choose a ', ...
               'distance other than cosine.'], []);
    end
    X = X ./ Xnorm(:,ones(1,p));
    additionalArg = [];
case 'cor' % Correlation
    X = X - repmat(mean(X,2),1,p);
    Xnorm = sqrt(sum(X.^2, 2));
    if min(Xnorm) <= eps(max(Xnorm))
        error('stats:pdist:InappropriateDistance',...
              ['Some points have small relative standard deviations, making ', ...
               'them effectively constant.\nEither remove those points, or ', ...
               'choose a distance other than correlation.'], []);
    end
    X = X ./ Xnorm(:,ones(1,p));
    additionalArg = [];
case 'spe'
    X = tiedrank(X')'; % treat rows as a series
    X = X - (p+1)/2; % subtract off the (constant) mean
    Xnorm = sqrt(sum(X.^2, 2));
    if min(Xnorm) <= eps(max(Xnorm))
        error('stats:pdist:InappropriateDistance',...
              ['Some points have too many ties, making them effectively ', ...
               'constant.\nEither remove those points, or choose a ', ...
               'distance other than rank correlation.'], []);
    end
    X = X ./ Xnorm(:,ones(1,p));
    additionalArg = [];
otherwise
    additionalArg = [];
end

% Call a mex file to compute distances for the standard distance measures
% and real double or single data.
if ~strcmp(dist,'usr') && isfloat(X) % ~usr => ~complex
    Y = pdistmex(X',dist,additionalArg);

% This M equivalent assumes real single or double.  It is not currently
% ever called, but it may be useful as a template for customization.
elseif ~strcmp(dist,'usr')
    if strmatch(dist, {'ham' 'jac' 'che'})
        nans = any(isnan(X),2);
    end
    outClass = class(X);
    Y = zeros(1,n*(n-1)./2, outClass);
    k = 1;
    for i = 1:n-1
        switch dist
        case 'euc'    % Euclidean
            dsq = zeros(n-i,1,outClass);
            for q = 1:p
                dsq = dsq + (X(i,q) - X((i+1):n,q)).^2;
            end
            Y(k:(k+n-i-1)) = sqrt(dsq);

        case 'seu'    % Standardized Euclidean
            wgts = additionalArg;
            dsq = zeros(n-i,1,outClass);
            for q = 1:p
                dsq = dsq + wgts(q) .* (X(i,q) - X((i+1):n,q)).^2;
            end
            Y(k:(k+n-i-1)) = sqrt(dsq);

        case 'cit'    % City Block
            d = zeros(n-i,1,outClass);
            for q = 1:p
                d = d + abs(X(i,q) - X((i+1):n,q));
            end
            Y(k:(k+n-i-1)) = d;

        case 'mah'    % Mahalanobis
            invcov = additionalArg;
            del = repmat(X(i,:),n-i,1) - X((i+1):n,:);
            dsq = sum((del*invcov).*del,2);
            Y(k:(k+n-i-1)) = sqrt(dsq);

        case 'min'    % Minkowski
            expon = additionalArg;
            dpow = zeros(n-i,1,outClass);
            for q = 1:p
                dpow = dpow + abs(X(i,q) - X((i+1):n,q)).^expon;
            end
            Y(k:(k+n-i-1)) = dpow .^ (1./expon);

        case {'cos' 'cor' 'spe'}   % Cosine, Correlation, Rank Correlation
            % This assumes that data have been appropriately preprocessed
            d = zeros(n-i,1,outClass);
            for q = 1:p
                d = d + (X(i,q).*X((i+1):n,q));
            end
            Y(k:(k+n-i-1)) = 1 - d;

        case 'ham'    % Hamming
            nesum = zeros(n-i,1,outClass);
            for q = 1:p
                nesum = nesum + (X(i,q) ~= X((i+1):n,q));
            end
            nesum(nans(i) | nans((i+1):n)) = NaN;
            Y(k:(k+n-i-1)) = nesum ./ p;

        case 'jac'    % Jaccard
            nzsum = zeros(n-i,1,outClass);
            nesum = zeros(n-i,1,outClass);
            for q = 1:p
                nz = (X(i,q) ~= 0 | X((i+1):n,q) ~= 0);
                ne = (X(i,q) ~= X((i+1):n,q));
                nzsum = nzsum + nz;
                nesum = nesum + (nz & ne);
            end
            nesum(nans(i) | nans((i+1):n)) = NaN;
            Y(k:(k+n-i-1)) = nesum ./ nzsum;

        case 'che'    % Chebychev
            dmax = zeros(n-i,1,outClass);
            for q = 1:p
                dmax = max(dmax, abs(X(i,q) - X((i+1):n,q)));
            end
            dmax(nans(i) | nans((i+1):n)) = NaN;
            Y(k:(k+n-i-1)) = dmax;

        end
        k = k + (n-i);
    end

% Compute distances for a caller-defined distance function.
else % if strcmp(dist,'usr')
    warning('stats:pdist:APIChanged', ...
            'The input arguments for caller-defined distance functions has\nchanged beginning in R14.  See the help for details.');
    try
        Y = feval(distfun,X(1,:),X(2,:),distargs{:})';
    catch
        [errMsg,errID] = lasterr;
        if strcmp('MATLAB:UndefinedFunction', errID) ...
                && ~isempty(strfind(errMsg, func2str(distfun)))
            error('stats:pdist:DistanceFunctionNotFound',...
                  'The distance function ''%s'' was not found.', func2str(distfun));
        end
        % Otherwise, let the catch block below generate the error message
        Y = [];
    end

    % Make the return have whichever numeric type the distance function
    % returns, or logical.
    if islogical(Y)
        Y = false(1,n*(n-1)./2);
    else % isnumeric
        Y = zeros(1,n*(n-1)./2, class(Y));
    end

    k = 1;
    for i = 1:n-1
        try
            Y(k:(k+n-i-1)) = feval(distfun,X(i,:),X((i+1):n,:),distargs{:})';
        catch
            [errMsg,errID] = lasterr;
            if isa(distfun, 'inline')
                error('stats:pdist:DistanceFunctionError',...
                      ['The inline distance function generated the following ', ...
                       'error:\n%s'], lasterr);
            else
                error('stats:pdist:DistanceFunctionError',...
                      ['The distance function ''%s'' generated the following ', ...
                       'error:\n%s'], func2str(distfun),lasterr);
            end
        end
        k = k + (n-i);
    end
end

>> A = rand(100,3)

A =

    0.9501    0.5828    0.4398
    0.2311    0.4235    0.3400
    0.6068    0.5155    0.3142
    0.4860    0.3340    0.3651
    0.8913    0.4329    0.3932
    0.7621    0.2259    0.5915
    0.4565    0.5798    0.1197
    0.0185    0.7604    0.0381
    0.8214    0.5298    0.4586
    0.4447    0.6405    0.8699
    0.6154    0.2091    0.9342
    0.7919    0.3798    0.2644
    0.9218    0.7833    0.1603
    0.7382    0.6808    0.8729
    0.1763    0.4611    0.2379
    0.4057    0.5678    0.6458
    0.9355    0.7942    0.9669
    0.9169    0.0592    0.6649
    0.4103    0.6029    0.8704
    0.8936    0.0503    0.0099
    0.0579    0.4154    0.1370
    0.3529    0.3050    0.8188
    0.8132    0.8744    0.4302
    0.0099    0.0150    0.8903
    0.1389    0.7680    0.7349
    0.2028    0.9708    0.6873
    0.1987    0.9901    0.3461
    0.6038    0.7889    0.1660
    0.2722    0.4387    0.1556
    0.1988    0.4983    0.1911
    0.0153    0.2140    0.4225
    0.7468    0.6435    0.8560
    0.4451    0.3200    0.4902
    0.9318    0.9601    0.8159
    0.4660    0.7266    0.4608
    0.4186    0.4120    0.4574
    0.8462    0.7446    0.4507
    0.5252    0.2679    0.4122
    0.2026    0.4399    0.9016
    0.6721    0.9334    0.0056
    0.8381    0.6833    0.2974
    0.0196    0.2126    0.0492
    0.6813    0.8392    0.6932
    0.3795    0.6288    0.6501
    0.8318    0.1338    0.9830
    0.5028    0.2071    0.5527
    0.7095    0.6072    0.4001
    0.4289    0.6299    0.1988
    0.3046    0.3705    0.6252
    0.1897    0.5751    0.7334
    0.1934    0.4514    0.3759
    0.6822    0.0439    0.0099
    0.3028    0.0272    0.4199
    0.5417    0.3127    0.7537
    0.1509    0.0129    0.7939
    0.6979    0.3840    0.9200
    0.3784    0.6831    0.8447
    0.8600    0.0928    0.3678
    0.8537    0.0353    0.6208
    0.5936    0.6124    0.7313
    0.4966    0.6085    0.1939
    0.8998    0.0158    0.9048
    0.8216    0.0164    0.5692
    0.6449    0.1901    0.6318
    0.8180    0.5869    0.2344
    0.6602    0.0576    0.5488
    0.3420    0.3676    0.9316
    0.2897    0.6315    0.3352
    0.3412    0.7176    0.6555
    0.5341    0.6927    0.3919
    0.7271    0.0841    0.6273
    0.3093    0.4544    0.6991
    0.8385    0.4418    0.3972
    0.5681    0.3533    0.4136
    0.3704    0.1536    0.6552
    0.7027    0.6756    0.8376
    0.5466    0.6992    0.3716
    0.4449    0.7275    0.4253
    0.6946    0.4784    0.5947
    0.6213    0.5548    0.5657
    0.7948    0.1210    0.7165
    0.9568    0.4508    0.5113
    0.5226    0.7159    0.7764
    0.8801    0.8928    0.4893
    0.1730    0.2731    0.1859
    0.9797    0.2548    0.7006
    0.2714    0.8656    0.9827
    0.2523    0.2324    0.8066
    0.8757    0.8049    0.7036
    0.7373    0.9084    0.4850
    0.1365    0.2319    0.1146
    0.0118    0.2393    0.6649
    0.8939    0.0498    0.3654
    0.1991    0.0784    0.1400
    0.2987    0.6408    0.5668
    0.6614    0.1909    0.8230
    0.2844    0.8439    0.6739
    0.4692    0.1739    0.9994
    0.0648    0.1708    0.9616
    0.9883    0.9943    0.0589

>> AA = zDistance(A);
>> BB = zDistance(A,A);
>> max(max(abs(AA-BB)))

ans =

     0

>> size(AA)

ans =

   100   100

>> size(BB)

ans =

   100   100

>> CC = pdist(A);
>> max(max(abs(AA-CC)))
??? Error using ==> minus
Matrix dimensions must agree.

>> size(CC)

ans =

           1        4950

>> help pdist
 PDIST Pairwise distance between observations.
    Y = PDIST(X) returns a vector Y containing the Euclidean distances
    between each pair of observations in the N-by-P data matrix X.  Rows of
    X correspond to observations, columns correspond to variables.  Y is a
    1-by-(N*(N-1)/2) row vector, corresponding to the N*(N-1)/2 pairs of
    observations in X.
 
    Y = PDIST(X, DISTANCE) computes Y using DISTANCE.  Choices are:
 
        'euclidean'   - Euclidean distance
        'seuclidean'  - Standardized Euclidean distance, each coordinate
                        in the sum of squares is inverse weighted by the
                        sample variance of that coordinate
        'cityblock'   - City Block distance
        'mahalanobis' - Mahalanobis distance
        'minkowski'   - Minkowski distance with exponent 2
        'cosine'      - One minus the cosine of the included angle
                        between observations (treated as vectors)
        'correlation' - One minus the sample linear correlation between
                        observations (treated as sequences of values).
        'spearman'    - One minus the sample Spearman's rank correlation
                        between observations (treated as sequences of values).
        'hamming'     - Hamming distance, percentage of coordinates
                        that differ
        'jaccard'     - One minus the Jaccard coefficient, the
                        percentage of nonzero coordinates that differ
        'chebychev'   - Chebychev distance (maximum coordinate difference)
        function      - A distance function specified using @, for
                        example @DISTFUN
 
    A distance function must be of the form
 
          function D = DISTFUN(XI, XJ, P1, P2, ...),
 
    taking as arguments a 1-by-P vector XI containing a single row of X, an
    M-by-P matrix XJ containing multiple rows of X, and zero or more
    additional problem-dependent arguments P1, P2, ..., and returning an
    M-by-1 vector of distances D, whose Kth element is the distance between
    the observations XI and XJ(K,:).
 
    Y = PDIST(X, DISTFUN, P1, P2, ...) passes the arguments P1, P2, ...
    directly to the function DISTFUN.
 
    Y = PDIST(X, 'minkowski', P) computes Minkowski distance using the
    positive scalar exponent P.
 
    The output Y is arranged in the order of ((1,2),(1,3),..., (1,N),
    (2,3),...(2,N),.....(N-1,N)), i.e. the lower left triangle of the full
    N-by-N distance matrix in column order.  To get the distance between
    the Ith and Jth observations (I < J), either use the formula
    Y((I-1)*(N-I/2)+J-I), or use the helper function Z = SQUAREFORM(Y),
    which returns an N-by-N square symmetric matrix, with the (I,J) entry
    equal to distance between observation I and observation J.
 
    Example:
 
       X = randn(100, 5);                 % some random points
       Y = pdist(X, 'euclidean');         % unweighted distance
       Wgts = [.1 .3 .3 .2 .1];           % coordinate weights
       Ywgt = pdist(X, @weucldist, Wgts); % weighted distance
 
       function d = weucldist(XI, XJ, W) % weighted euclidean distance
       d = sqrt((repmat(XI,size(XJ,1),1)-XJ).^2 * W');
 
    See also squareform, linkage, silhouette.

    Overloaded functions or methods (ones with the same name in other directories)
       help phytree/pdist.m

    Reference page in Help browser
       doc pdist



>> CC = squareform(pdist(A));
>> max(max(abs(AA-CC)))

ans =

  2.9802e-008

>> max(diag(AA))

ans =

  2.9802e-008

>> diag(AA) = zeros(1,100)
??? Subscript indices must either be real positive integers or logicals.

>> help diag
 DIAG Diagonal matrices and diagonals of a matrix.
    DIAG(V,K) when V is a vector with N components is a square matrix
    of order N+ABS(K) with the elements of V on the K-th diagonal. K = 0
    is the main diagonal, K > 0 is above the main diagonal and K < 0
    is below the main diagonal. 
 
    DIAG(V) is the same as DIAG(V,0) and puts V on the main diagonal.
 
    DIAG(X,K) when X is a matrix is a column vector formed from
    the elements of the K-th diagonal of X.
 
    DIAG(X) is the main diagonal of X. DIAG(DIAG(X)) is a diagonal matrix.
 
    Example
       m = 5;
       diag(-m:m) + diag(ones(2*m,1),1) + diag(ones(2*m,1),-1)
    produces a tridiagonal matrix of order 2*m+1.
 
    See also spdiags, triu, tril.


    Reference page in Help browser
       doc diag



>> AA = zDistance(A);
>> max(max(abs(AA-CC)))

ans =

   92.8462

>> AA = zDistance(A);
>> max(max(abs(AA-CC)))

ans =

  7.9312e-015

>> AA = zDistance(A);
>> max(max(abs(AA-CC)))

ans =

  7.9312e-015

>> AA = zDistance(A);
??? Error using ==> plus
Matrix dimensions must agree.

Error in ==> zDistance at 17
  D = sqrt(X + X' - 2*G);             % |u-v| = |u|^2 + |v|^2 - 2 u . v

>> AA = zDistance(A);
>> max(max(abs(AA-CC)))

ans =

  4.2050e-015

>> AA = zDistance(A);

ans =

     0

>> type pdist

function Y = pdist(X,dist,varargin)
%PDIST Pairwise distance between observations.
%   Y = PDIST(X) returns a vector Y containing the Euclidean distances
%   between each pair of observations in the N-by-P data matrix X.  Rows of
%   X correspond to observations, columns correspond to variables.  Y is a
%   1-by-(N*(N-1)/2) row vector, corresponding to the N*(N-1)/2 pairs of
%   observations in X.
%
%   Y = PDIST(X, DISTANCE) computes Y using DISTANCE.  Choices are:
%
%       'euclidean'   - Euclidean distance
%       'seuclidean'  - Standardized Euclidean distance, each coordinate
%                       in the sum of squares is inverse weighted by the
%                       sample variance of that coordinate
%       'cityblock'   - City Block distance
%       'mahalanobis' - Mahalanobis distance
%       'minkowski'   - Minkowski distance with exponent 2
%       'cosine'      - One minus the cosine of the included angle
%                       between observations (treated as vectors)
%       'correlation' - One minus the sample linear correlation between
%                       observations (treated as sequences of values).
%       'spearman'    - One minus the sample Spearman's rank correlation
%                       between observations (treated as sequences of values).
%       'hamming'     - Hamming distance, percentage of coordinates
%                       that differ
%       'jaccard'     - One minus the Jaccard coefficient, the
%                       percentage of nonzero coordinates that differ
%       'chebychev'   - Chebychev distance (maximum coordinate difference)
%       function      - A distance function specified using @, for
%                       example @DISTFUN
%
%   A distance function must be of the form
%
%         function D = DISTFUN(XI, XJ, P1, P2, ...),
%
%   taking as arguments a 1-by-P vector XI containing a single row of X, an
%   M-by-P matrix XJ containing multiple rows of X, and zero or more
%   additional problem-dependent arguments P1, P2, ..., and returning an
%   M-by-1 vector of distances D, whose Kth element is the distance between
%   the observations XI and XJ(K,:).
%
%   Y = PDIST(X, DISTFUN, P1, P2, ...) passes the arguments P1, P2, ...
%   directly to the function DISTFUN.
%
%   Y = PDIST(X, 'minkowski', P) computes Minkowski distance using the
%   positive scalar exponent P.
%
%   The output Y is arranged in the order of ((1,2),(1,3),..., (1,N),
%   (2,3),...(2,N),.....(N-1,N)), i.e. the lower left triangle of the full
%   N-by-N distance matrix in column order.  To get the distance between
%   the Ith and Jth observations (I < J), either use the formula
%   Y((I-1)*(N-I/2)+J-I), or use the helper function Z = SQUAREFORM(Y),
%   which returns an N-by-N square symmetric matrix, with the (I,J) entry
%   equal to distance between observation I and observation J.
%
%   Example:
%
%      X = randn(100, 5);                 % some random points
%      Y = pdist(X, 'euclidean');         % unweighted distance
%      Wgts = [.1 .3 .3 .2 .1];           % coordinate weights
%      Ywgt = pdist(X, @weucldist, Wgts); % weighted distance
%
%      function d = weucldist(XI, XJ, W) % weighted euclidean distance
%      d = sqrt((repmat(XI,size(XJ,1),1)-XJ).^2 * W');
%
%   See also SQUAREFORM, LINKAGE, SILHOUETTE.

%   An example of distance for data with missing elements:
%
%      X = randn(100, 5);     % some random points
%      X(unidrnd(prod(size(X)),1,20)) = NaN; % scatter in some NaNs
%      D = pdist(X, @naneucdist);
%
%      function d = naneucdist(XI, XJ) % euclidean distance, ignoring NaNs
%      [m,p] = size(XJ);
%      sqdx = (repmat(XI,m,1) - XJ) .^ 2;
%      pstar = sum(~isnan(sqdx),2); % correction for missing coords
%      pstar(pstar == 0) = NaN;
%      d = sqrt(nansum(sqdx,2) .* p ./ pstar);
%
%
%   For a large number of observations, it is sometimes faster to compute
%   the distances by looping over coordinates of the data (though the code
%   is more complicated):
%
%      function d = nanhamdist(XI, XJ) % hamming distance, ignoring NaNs
%      [m,p] = size(XJ);
%      nesum = zeros(m,1);
%      pstar = zeros(m,1);
%      for q = 1:p
%          notnan = ~(isnan((XI(q)) | isnan(XJ(:,q)));
%          nesum = nesum + (XI(q) ~= XJ(:,q)) & notnan;
%          pstar = pstar + notnan;
%      end
%      nesum(any() | nans((i+1):n)) = NaN;
%      Y(k:(k+n-i-1)) = nesum ./ pstar;

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.15.4.10 $ $Date: 2004/08/20 20:06:05 $

if nargin < 2
    dist = 'euc';
else
    if ischar(dist)
        methods = {'euclidean'; 'seuclidean'; 'cityblock'; 'chebychev'; ...
                   'mahalanobis'; 'minkowski'; 'cosine'; 'correlation'; ...
                   'spearman'; 'hamming'; 'jaccard'};
        i = strmatch(lower(dist), methods);
        if length(i) > 1
            error('stats:pdist:BadDistance',...
                  'Ambiguous ''DISTANCE'' argument:  %s.', dist);
        elseif isempty(i)
            % Assume an unrecognized string is a user-supplied distance
            % function name, change it to a handle.
            distfun = str2func(dist);
            distargs = varargin;
            dist = 'usr';
        else
            dist = lower(methods{i}(1:3));
        end
    elseif isa(dist, 'function_handle') ||  isa(dist, 'inline')
        distfun = dist;
        distargs = varargin;
        dist = 'usr';
    else
        error('stats:pdist:BadDistance',...
              'The ''DISTANCE'' argument must be a string or a function.');
    end
end

% Integer/logical/char/anything data may be handled by a caller-defined
% distance function, otherwise it is converted to double.  Complex floating
% point data must also be handled by a caller-defined distance function.
if ~strcmp(dist,'usr')
    if ~isfloat(X)
        warning('stats:pdist:DataConversion', ...
                'Converting %s data to double.',class(X));
        X = double(X);
    elseif any(imag(X(:)))
        error('stats:pdist:InvalidData', ...
              'PDIST does not accept complex data for built-in distances.');
    end
end

[n,p] = size(X);

% Degenerate case, just return an empty of the proper size.
if n < 2
    if ~strcmp(dist,'usr')
        Y = zeros(1,0,class(X)); % X was single/double, or cast to double
    else
        Y = zeros(1,0);
    end
    return;
end

switch dist
case 'seu' % Standardized Euclidean weights by coordinate variance
    additionalArg = 1 ./ var(X)';
case 'mah' % Mahalanobis
    additionalArg = cov(X) \ eye(p); %inv(cov(X));
case 'min' % Minkowski distance needs a third argument
    if nargin < 3  % use default value for exponent
        additionalArg = 2;
    elseif varargin{1} > 0
        additionalArg = varargin{1}; % get exponent from input args
    else
        error('stats:pdist:InvalidExponent',...
              'The exponent for the Minkowski metric must be positive.');
    end
case 'cos' % Cosine
    Xnorm = sqrt(sum(X.^2, 2));
    if min(Xnorm) <= eps(max(Xnorm))
        error('stats:pdist:InappropriateDistance',...
              ['Some points have small relative magnitudes, making them ', ...
               'effectively zero.\nEither remove those points, or choose a ', ...
               'distance other than cosine.'], []);
    end
    X = X ./ Xnorm(:,ones(1,p));
    additionalArg = [];
case 'cor' % Correlation
    X = X - repmat(mean(X,2),1,p);
    Xnorm = sqrt(sum(X.^2, 2));
    if min(Xnorm) <= eps(max(Xnorm))
        error('stats:pdist:InappropriateDistance',...
              ['Some points have small relative standard deviations, making ', ...
               'them effectively constant.\nEither remove those points, or ', ...
               'choose a distance other than correlation.'], []);
    end
    X = X ./ Xnorm(:,ones(1,p));
    additionalArg = [];
case 'spe'
    X = tiedrank(X')'; % treat rows as a series
    X = X - (p+1)/2; % subtract off the (constant) mean
    Xnorm = sqrt(sum(X.^2, 2));
    if min(Xnorm) <= eps(max(Xnorm))
        error('stats:pdist:InappropriateDistance',...
              ['Some points have too many ties, making them effectively ', ...
               'constant.\nEither remove those points, or choose a ', ...
               'distance other than rank correlation.'], []);
    end
    X = X ./ Xnorm(:,ones(1,p));
    additionalArg = [];
otherwise
    additionalArg = [];
end

% Call a mex file to compute distances for the standard distance measures
% and real double or single data.
if ~strcmp(dist,'usr') && isfloat(X) % ~usr => ~complex
    Y = pdistmex(X',dist,additionalArg);

% This M equivalent assumes real single or double.  It is not currently
% ever called, but it may be useful as a template for customization.
elseif ~strcmp(dist,'usr')
    if strmatch(dist, {'ham' 'jac' 'che'})
        nans = any(isnan(X),2);
    end
    outClass = class(X);
    Y = zeros(1,n*(n-1)./2, outClass);
    k = 1;
    for i = 1:n-1
        switch dist
        case 'euc'    % Euclidean
            dsq = zeros(n-i,1,outClass);
            for q = 1:p
                dsq = dsq + (X(i,q) - X((i+1):n,q)).^2;
            end
            Y(k:(k+n-i-1)) = sqrt(dsq);

        case 'seu'    % Standardized Euclidean
            wgts = additionalArg;
            dsq = zeros(n-i,1,outClass);
            for q = 1:p
                dsq = dsq + wgts(q) .* (X(i,q) - X((i+1):n,q)).^2;
            end
            Y(k:(k+n-i-1)) = sqrt(dsq);

        case 'cit'    % City Block
            d = zeros(n-i,1,outClass);
            for q = 1:p
                d = d + abs(X(i,q) - X((i+1):n,q));
            end
            Y(k:(k+n-i-1)) = d;

        case 'mah'    % Mahalanobis
            invcov = additionalArg;
            del = repmat(X(i,:),n-i,1) - X((i+1):n,:);
            dsq = sum((del*invcov).*del,2);
            Y(k:(k+n-i-1)) = sqrt(dsq);

        case 'min'    % Minkowski
            expon = additionalArg;
            dpow = zeros(n-i,1,outClass);
            for q = 1:p
                dpow = dpow + abs(X(i,q) - X((i+1):n,q)).^expon;
            end
            Y(k:(k+n-i-1)) = dpow .^ (1./expon);

        case {'cos' 'cor' 'spe'}   % Cosine, Correlation, Rank Correlation
            % This assumes that data have been appropriately preprocessed
            d = zeros(n-i,1,outClass);
            for q = 1:p
                d = d + (X(i,q).*X((i+1):n,q));
            end
            Y(k:(k+n-i-1)) = 1 - d;

        case 'ham'    % Hamming
            nesum = zeros(n-i,1,outClass);
            for q = 1:p
                nesum = nesum + (X(i,q) ~= X((i+1):n,q));
            end
            nesum(nans(i) | nans((i+1):n)) = NaN;
            Y(k:(k+n-i-1)) = nesum ./ p;

        case 'jac'    % Jaccard
            nzsum = zeros(n-i,1,outClass);
            nesum = zeros(n-i,1,outClass);
            for q = 1:p
                nz = (X(i,q) ~= 0 | X((i+1):n,q) ~= 0);
                ne = (X(i,q) ~= X((i+1):n,q));
                nzsum = nzsum + nz;
                nesum = nesum + (nz & ne);
            end
            nesum(nans(i) | nans((i+1):n)) = NaN;
            Y(k:(k+n-i-1)) = nesum ./ nzsum;

        case 'che'    % Chebychev
            dmax = zeros(n-i,1,outClass);
            for q = 1:p
                dmax = max(dmax, abs(X(i,q) - X((i+1):n,q)));
            end
            dmax(nans(i) | nans((i+1):n)) = NaN;
            Y(k:(k+n-i-1)) = dmax;

        end
        k = k + (n-i);
    end

% Compute distances for a caller-defined distance function.
else % if strcmp(dist,'usr')
    warning('stats:pdist:APIChanged', ...
            'The input arguments for caller-defined distance functions has\nchanged beginning in R14.  See the help for details.');
    try
        Y = feval(distfun,X(1,:),X(2,:),distargs{:})';
    catch
        [errMsg,errID] = lasterr;
        if strcmp('MATLAB:UndefinedFunction', errID) ...
                && ~isempty(strfind(errMsg, func2str(distfun)))
            error('stats:pdist:DistanceFunctionNotFound',...
                  'The distance function ''%s'' was not found.', func2str(distfun));
        end
        % Otherwise, let the catch block below generate the error message
        Y = [];
    end

    % Make the return have whichever numeric type the distance function
    % returns, or logical.
    if islogical(Y)
        Y = false(1,n*(n-1)./2);
    else % isnumeric
        Y = zeros(1,n*(n-1)./2, class(Y));
    end

    k = 1;
    for i = 1:n-1
        try
            Y(k:(k+n-i-1)) = feval(distfun,X(i,:),X((i+1):n,:),distargs{:})';
        catch
            [errMsg,errID] = lasterr;
            if isa(distfun, 'inline')
                error('stats:pdist:DistanceFunctionError',...
                      ['The inline distance function generated the following ', ...
                       'error:\n%s'], lasterr);
            else
                error('stats:pdist:DistanceFunctionError',...
                      ['The distance function ''%s'' generated the following ', ...
                       'error:\n%s'], func2str(distfun),lasterr);
            end
        end
        k = k + (n-i);
    end
end

>> AA = zDistance(A);
>> AA = zDistance(A);
>> max(max(abs(AA-CC)))

ans =

  4.2050e-015

>> max(diag(AA))

ans =

     0

>> max(diag(CC))

ans =

     0

>> Exem
Wrote PairExemplarPDB.pdb
>> zDisplayNT('PairExemplarPDB')
Read PairExemplarPDB.pdb for header information
Classifying interactions
Found   332 bases within 15 Angstroms from 3d structure
Found   266 pairs that are possibly interacting
Classification took 0.03 minutes, or 7797 classifications per minute
Saved PairExemplarPDB.mat
>> Exemplar

Exemplar = 

32x16 struct array with fields:
    Filename
    Class
    NT1
    NT2
    Pair
    Count
    Displ
    Rot
    Base1Index
    Base2Index
    Base1
    Base2

>> Exemplar(1,1)

ans = 

      Filename: '1j5e'
         Class: 1
           NT1: [1x1 struct]
           NT2: [1x1 struct]
          Pair: [1x1 struct]
         Count: 4
         Displ: [5.3239 9.7762 -0.2345]
           Rot: [3x3 double]
    Base1Index: 654
    Base2Index: 694
         Base1: '675'
         Base2: '715'

>> Exemplar(1,1)help fix
>> help fix
 FIX    Round towards zero.
    FIX(X) rounds the elements of X to the nearest integers
    towards zero.
 
    See also floor, round, ceil.


    Reference page in Help browser
       doc fix



>> zWriteExemplarPDB(1)
??? Error: File: zWriteExemplarPDB.m Line: 60 Column: 107
Unbalanced or misused parentheses or brackets.

>> zWriteExemplarPDB(1)
??? Undefined function or variable "c".

Error in ==> zWriteExemplarPDB at 62
        a = zWriteNucleotidePDB(fid,E.NT1,a,c,R,sh);

>> zWriteExemplarPDB(1)
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
Wrote PairExemplarPDB.pdb
>> zWriteExemplarPDB(1)
>> zWriteExemplarPDB(1)
>> zWriteExemplarPDB(1)
>> zWriteExemplarPDB(1)
>> zWriteExemplarPDB(1)
>> zDisplayNT('PairExemplar_tWW_UU_2AW4_1782B_2586B')
Read PairExemplar_tWW_UU_2AW4_1782B_2586B.pdb for header information
Classifying interactions
Found     1 bases within 15 Angstroms from 3d structure
Found     1 pairs that are possibly interacting
Classification took 0.01 minutes, or   85 classifications per minute
Saved PairExemplar_tWW_UU_2AW4_1782B_2586B.mat
>> clf
>> zDisplayNT('PairExemplar_tWW_UU_2AW4_1782B_2586B')
Loaded PairExemplar_tWW_UU_2AW4_1782B_2586B.mat
>> zWriteExemplarPDB(1)
>> zDisplayNT('PairExemplar_tWW_GC_1h4q_15T_48T')
Read PairExemplar_tWW_GC_1h4q_15T_48T.pdb for header information
Classifying interactions
Found     1 bases within 15 Angstroms from 3d structure
Found     1 pairs that are possibly interacting
Classification took 0.01 minutes, or   85 classifications per minute
Saved PairExemplar_tWW_GC_1h4q_15T_48T.mat
>> zWriteExemplarPDB(1)
>> zDisplayNT('PairExemplar_tWW_CU_1s72_13940_14320')
Read PairExemplar_tWW_CU_1s72_13940_14320.pdb for header information
Classifying interactions
Found     1 bases within 15 Angstroms from 3d structure
Found     1 pairs that are possibly interacting
Classification took 0.01 minutes, or   85 classifications per minute
Saved PairExemplar_tWW_CU_1s72_13940_14320.mat
>> clf
>> zDisplayNT('PairExemplar_tWW_CU_1s72_13940_14320')
Loaded PairExemplar_tWW_CU_1s72_13940_14320.mat
>> File = zGetNTData('HighResolution_list')
??? Input argument "ReadCode" is undefined.

Error in ==> zGetNTData at 41
  if (ReadCode < 4) & (exist(strcat(Filename,'.mat'),'file') > 0),

>> File = zGetNTData('HighResolution_list',0)
Read HighResolution_list.pdb for header information
??? Subscript indices must either be real positive integers or logicals.

Error in ==> mSynList at 16
synlist(1,length(File.NT))=0;

Error in ==> zReadandAnalyze at 283
SynList = mSynList(File);

Error in ==> zGetNTData at 54
      File = zReadandAnalyze(Filename);

>> File = zGetNTData('HighResolution_list',0)
Read HighResolution_list.pdb for header information
??? Output argument "Files" (and maybe others) not assigned during call to "C:\Documents and Settings\zirbel\My Documents\FR3D\zGetNTData.m (zGetNTData)".

Error in ==> zGetNTData at 16
path(path,pwd);

>> File = zAddNTData('HighResolution_list')
Read HighResolution_list.pdb
Loaded 1s72.mat
Loaded 1j5e.mat
Loaded 2aw4.mat
Loaded 2avy.mat

File = 

1x4 struct array with fields:
    CI
    ClassVersion
    Comment
    Distance
    Edge
    Filename
    HandClass
    Header
    Inter
    Modified
    NT
    NumNT
    Pair

>> File = zAddNTData('HighResolution_list',0,File)
Read HighResolution_list.pdb
Read 1q96.pdb for header information
Classifying interactions
Found   598 bases within 15 Angstroms from 3d structure
Found   238 pairs that are possibly interacting
Classification took 0.03 minutes, or 8308 classifications per minute
Saved 1q96.mat
Loaded 157d.mat
Loaded 165d.mat
Loaded 1a9n.mat
Loaded 1bmv.mat
Loaded 1d96.mat
Loaded 1dk1.mat
Loaded 1dpl.mat
Loaded 1fir.mat
Loaded 1hmh.mat
Loaded 1i5l.mat
Loaded 1i2y.mat
Loaded 1id9.mat
Loaded 1idw.mat
Loaded 1iha.mat
Loaded 1kuq.mat
Classifying interactions
Found   405 bases within 15 Angstroms from 3d structure
Found   165 pairs that are possibly interacting
Classification took 0.02 minutes, or 7119 classifications per minute
Saved 1kuq.mat
Loaded 1l2x.mat
Classifying interactions
Found   196 bases within 15 Angstroms from 3d structure
Found    92 pairs that are possibly interacting
Classification took 0.02 minutes, or 5047 classifications per minute
Saved 1l2x.mat
Loaded 1l3d.mat
Classifying interactions
Found   195 bases within 15 Angstroms from 3d structure
Found    89 pairs that are possibly interacting
Classification took 0.02 minutes, or 4953 classifications per minute
Saved 1l3d.mat
Read 1mme.pdb for header information
Classifying interactions
Found   628 bases within 15 Angstroms from 3d structure
Found   232 pairs that are possibly interacting
Classification took 0.03 minutes, or 8405 classifications per minute
Saved 1mme.mat
Read 1mwl.pdb for header information
Classifying interactions
Found   262 bases within 15 Angstroms from 3d structure
Found    98 pairs that are possibly interacting
Classification took 0.02 minutes, or 5376 classifications per minute
Saved 1mwl.mat
Read 1nlc.pdb for header information
Classifying interactions
Found   342 bases within 15 Angstroms from 3d structure
Found   115 pairs that are possibly interacting
Classification took 0.02 minutes, or 5968 classifications per minute
Saved 1nlc.mat
Read 1nuj.pdb for header information
Classifying interactions
Found   675 bases within 15 Angstroms from 3d structure
Found   270 pairs that are possibly interacting
Classification took 0.03 minutes, or 8569 classifications per minute
Saved 1nuj.mat
Read 1nuv.pdb for header information
Classifying interactions
Found   677 bases within 15 Angstroms from 3d structure
Found   269 pairs that are possibly interacting
Classification took 0.03 minutes, or 8537 classifications per minute
Saved 1nuv.mat
Read 1o3z.pdb for header information
Classifying interactions
Found   345 bases within 15 Angstroms from 3d structure
Found   113 pairs that are possibly interacting
Classification took 0.02 minutes, or 5786 classifications per minute
Saved 1o3z.mat
Read 1osu.pdb for header information
Classifying interactions
Found    13 bases within 15 Angstroms from 3d structure
Found     5 pairs that are possibly interacting
Classification took 0.01 minutes, or  427 classifications per minute
Saved 1osu.mat
Read 1q9a.pdb for header information
Classifying interactions
Found   173 bases within 15 Angstroms from 3d structure
Found    67 pairs that are possibly interacting
Classification took 0.02 minutes, or 4150 classifications per minute
Saved 1q9a.mat
Read 1rpu.pdb for header information
Classifying interactions
Found   274 bases within 15 Angstroms from 3d structure
Found   103 pairs that are possibly interacting
Classification took 0.02 minutes, or 5732 classifications per minute
Saved 1rpu.mat
Read 1sdr.pdb for header information
Classifying interactions
Found   365 bases within 15 Angstroms from 3d structure
Found   142 pairs that are possibly interacting
Classification took 0.02 minutes, or 6570 classifications per minute
Saved 1sdr.mat
Read 1ser.pdb for header information
Classifying interactions
Found   432 bases within 15 Angstroms from 3d structure
Found   162 pairs that are possibly interacting
Classification took 0.02 minutes, or 7069 classifications per minute
Saved 1ser.mat
Read 1t03.pdb for header information
Classifying interactions
Found   199 bases within 15 Angstroms from 3d structure
Found    63 pairs that are possibly interacting
Classification took 0.02 minutes, or 4100 classifications per minute
Saved 1t03.mat
Read 1t0e.pdb for header information
Classifying interactions
Found   229 bases within 15 Angstroms from 3d structure
Found    81 pairs that are possibly interacting
Classification took 0.02 minutes, or 4713 classifications per minute
Saved 1t0e.mat
Read 1tn2.pdb for header information
Classifying interactions
Found   410 bases within 15 Angstroms from 3d structure
Found   150 pairs that are possibly interacting
Classification took 0.02 minutes, or 6698 classifications per minute
Saved 1tn2.mat
Read 1tra.pdb for header information
Classifying interactions
Found   404 bases within 15 Angstroms from 3d structure
Found   149 pairs that are possibly interacting
Classification took 0.02 minutes, or 6653 classifications per minute
Saved 1tra.mat
Read 1wvd.pdb for header information
Classifying interactions
Found   346 bases within 15 Angstroms from 3d structure
Found   117 pairs that are possibly interacting
Classification took 0.02 minutes, or 5990 classifications per minute
Saved 1wvd.mat
Read 1y6s.pdb for header information
Classifying interactions
Found   341 bases within 15 Angstroms from 3d structure
Found   113 pairs that are possibly interacting
Classification took 0.02 minutes, or 5786 classifications per minute
Saved 1y6s.mat
Read 1y6t.pdb for header information
Classifying interactions
Found   344 bases within 15 Angstroms from 3d structure
Found   115 pairs that are possibly interacting
Classification took 0.02 minutes, or 5811 classifications per minute
Saved 1y6t.mat
Read 1y73.pdb for header information
Classifying interactions
Found   346 bases within 15 Angstroms from 3d structure
Found   112 pairs that are possibly interacting
Classification took 0.02 minutes, or 5734 classifications per minute
Saved 1y73.mat
Read 1y90.pdb for header information
Classifying interactions
Found   322 bases within 15 Angstroms from 3d structure
Found   103 pairs that are possibly interacting
Classification took 0.02 minutes, or 5650 classifications per minute
Saved 1y90.mat
Read 1y95.pdb for header information
Classifying interactions
Found   342 bases within 15 Angstroms from 3d structure
Found   114 pairs that are possibly interacting
Classification took 0.02 minutes, or 5837 classifications per minute
Saved 1y95.mat
Read 1y99.pdb for header information
Classifying interactions
Found   348 bases within 15 Angstroms from 3d structure
Found   115 pairs that are possibly interacting
Classification took 0.02 minutes, or 5811 classifications per minute
Saved 1y99.mat
Read 205d.pdb for header information
Classifying interactions
Found   160 bases within 15 Angstroms from 3d structure
Found    64 pairs that are possibly interacting
Classification took 0.02 minutes, or 4029 classifications per minute
Saved 205d.mat
Read 247d.pdb for header information
Classifying interactions
Found    16 bases within 15 Angstroms from 3d structure
Found     6 pairs that are possibly interacting
Classification took 0.01 minutes, or  501 classifications per minute
Saved 247d.mat
Read 248d.pdb for header information
Classifying interactions
Found    17 bases within 15 Angstroms from 3d structure
Found     6 pairs that are possibly interacting
Classification took 0.01 minutes, or  512 classifications per minute
Saved 248d.mat
Read 255d.pdb for header information
Classifying interactions
Found    32 bases within 15 Angstroms from 3d structure
Found    11 pairs that are possibly interacting
Classification took 0.01 minutes, or  899 classifications per minute
Saved 255d.mat
Read 2tra.pdb for header information
Classifying interactions
Found   471 bases within 15 Angstroms from 3d structure
Found   176 pairs that are possibly interacting
Classification took 0.02 minutes, or 7427 classifications per minute
Saved 2tra.mat
Read 3tra.pdb for header information
Classifying interactions
Found   465 bases within 15 Angstroms from 3d structure
Found   167 pairs that are possibly interacting
Classification took 0.02 minutes, or 7287 classifications per minute
Saved 3tra.mat
Read 413d.pdb for header information
Classifying interactions
Found    35 bases within 15 Angstroms from 3d structure
Found    12 pairs that are possibly interacting
Classification took 0.01 minutes, or 1002 classifications per minute
Saved 413d.mat
Read 437d.pdb for header information
Classifying interactions
Found   196 bases within 15 Angstroms from 3d structure
Found    92 pairs that are possibly interacting
Classification took 0.02 minutes, or 4976 classifications per minute
Saved 437d.mat
Read 462d.pdb for header information
Classifying interactions
Found   344 bases within 15 Angstroms from 3d structure
Found   114 pairs that are possibly interacting
Classification took 0.02 minutes, or 5837 classifications per minute
Saved 462d.mat
Read 483d.pdb for header information
Classifying interactions
Found   173 bases within 15 Angstroms from 3d structure
Found    67 pairs that are possibly interacting
Classification took 0.02 minutes, or 4150 classifications per minute
Saved 483d.mat
Read 4tna.pdb for header information
Classifying interactions
Found   410 bases within 15 Angstroms from 3d structure
Found   150 pairs that are possibly interacting
Classification took 0.02 minutes, or 6776 classifications per minute
Saved 4tna.mat
Read 4tra.pdb for header information
Classifying interactions
Found   405 bases within 15 Angstroms from 3d structure
Found   148 pairs that are possibly interacting
Classification took 0.02 minutes, or 6686 classifications per minute
Saved 4tra.mat
Read 6tna.pdb for header information
Classifying interactions
Found   408 bases within 15 Angstroms from 3d structure
Found   152 pairs that are possibly interacting
Classification took 0.02 minutes, or 6709 classifications per minute
Saved 6tna.mat
Loaded 1a34.mat
Loaded 1apg.mat
Loaded 1aq3.mat
Loaded 1aq4.mat
Loaded 1asy.mat
Loaded 1asz.mat
Loaded 1av6.mat
Loaded 1b23.mat
Loaded 1b7f.mat
Loaded 1br3.mat
Loaded 1by4.mat
Loaded 1c04.mat
Loaded 1c0a.mat
Loaded 1csl.mat
Loaded 1cvj.mat
Loaded 1cwp.mat
Loaded 1cx0.mat
Loaded 1d4r.mat
Loaded 1d9h.mat
Loaded 1ddl.mat
Loaded 1ddy.mat
Loaded 1dfu.mat
Loaded 1di2.mat
Loaded 1dno.mat
Loaded 1dnt.mat
Loaded 1dnx.mat
Loaded 1dqf.mat
Loaded 1dqh.mat
Loaded 1drz.mat
Loaded 1duh.mat
Loaded 1dul.mat
Loaded 1duq.mat
Loaded 1dzs.mat
Loaded 1e6t.mat
Loaded 1e7k.mat
Loaded 1e7x.mat
Loaded 1e8o.mat
Loaded 1ec6.mat
Loaded 1efo.mat
Loaded 1efw.mat
Loaded 1egk.mat
Loaded 1ehz.mat
Loaded 1eiy.mat
Loaded 1et4.mat
Loaded 1euq.mat
Loaded 1euy.mat
Loaded 1evp.mat
Loaded 1evv.mat
Loaded 1exd.mat
Loaded 1f1t.mat
Loaded 1f27.mat
Loaded 1f7u.mat
Loaded 1f7v.mat
Loaded 1f8v.mat
Loaded 1feu.mat
Loaded 1fix.mat
Loaded 1fuf.mat
Loaded 1fxl.mat
Loaded 1g1x.mat
Loaded 1g2e.mat
Loaded 1g2j.mat
Loaded 1g4q.mat
Loaded 1g59.mat
Loaded 1gax.mat
Loaded 1gid.mat
Loaded 1gkv.mat
Loaded 1gkw.mat
Loaded 1grz.mat
Loaded 1gtf.mat
Loaded 1gtn.mat
Loaded 1gtr.mat
Loaded 1gts.mat
Loaded 1h1k.mat
Loaded 1h2d.mat
Loaded 1h38.mat
Loaded 1h3e.mat
Loaded 1h4q.mat
Loaded 1h4s.mat
Loaded 1h8j.mat
Loaded 1hc8.mat
Loaded 1hdw.mat
Loaded 1he0.mat
Loaded 1he6.mat
Loaded 1hq1.mat
Loaded 1hr2.mat
Loaded 1hvu.mat
Loaded 1hys.mat
Loaded 1i0f.mat
Loaded 1i0g.mat
Loaded 1i0j.mat
Loaded 1i0k.mat
Loaded 1i0m.mat
Loaded 1i0n.mat
Loaded 1i0o.mat
Loaded 1i0p.mat
Loaded 1i0q.mat
Loaded 1i2x.mat
Loaded 1i6h.mat
Loaded 1i6u.mat
Loaded 1i9v.mat
Loaded 1i9x.mat
Loaded 1icg.mat
Loaded 1ik5.mat
Loaded 1il2.mat
Loaded 1ivs.mat
Loaded 1j1u.mat
Loaded 1j2b.mat
Loaded 1j6s.mat
Loaded 1j7t.mat
Loaded 1j8g.mat
Loaded 1j9h.mat
Loaded 1jb8.mat
Loaded 1jbr.mat
Loaded 1jbs.mat
Loaded 1jbt.mat
Loaded 1jid.mat
Loaded 1jjm.mat
Loaded 1jjn.mat
Loaded 1jo2.mat
Loaded 1jzv.mat
Loaded 1k8w.mat
Loaded 1k9w.mat
Loaded 1kd3.mat
Loaded 1kd4.mat
Loaded 1kd5.mat
Loaded 1kfo.mat
Loaded 1kh6.mat
Loaded 1knz.mat
Loaded 1kog.mat
Loaded 1kq2.mat
Loaded 1kuo.mat
Classifying interactions
Found   109 bases within 15 Angstroms from 3d structure
Found    47 pairs that are possibly interacting
Classification took 0.01 minutes, or 3281 classifications per minute
Saved 1kuo.mat
Loaded 1kxk.mat
Classifying interactions
Found   513 bases within 15 Angstroms from 3d structure
Found   185 pairs that are possibly interacting
Classification took 0.02 minutes, or 7722 classifications per minute
Saved 1kxk.mat
Loaded 1l3z.mat
Classifying interactions
Found    49 bases within 15 Angstroms from 3d structure
Found    19 pairs that are possibly interacting
Classification took 0.01 minutes, or 1459 classifications per minute
Saved 1l3z.mat
Loaded 1l8v.mat
Classifying interactions
Found  2910 bases within 15 Angstroms from 3d structure
Found  1198 pairs that are possibly interacting
Classification took 0.10 minutes, or 11887 classifications per minute
Saved 1l8v.mat
Loaded 1l9a.mat
Classifying interactions
Found   984 bases within 15 Angstroms from 3d structure
Found   408 pairs that are possibly interacting
Classification took 0.04 minutes, or 10043 classifications per minute
Saved 1l9a.mat
Loaded 1laj.mat
Classifying interactions
Found     3 bases within 15 Angstroms from 3d structure
Found     1 pairs that are possibly interacting
Classification took 0.01 minutes, or   89 classifications per minute
Saved 1laj.mat
Loaded 1lc4.mat
Classifying interactions
Found   272 bases within 15 Angstroms from 3d structure
Found   100 pairs that are possibly interacting
Classification took 0.02 minutes, or 5408 classifications per minute
Saved 1lc4.mat
Loaded 1lng.mat
Classifying interactions
Found   744 bases within 15 Angstroms from 3d structure
Found   301 pairs that are possibly interacting
Classification took 0.03 minutes, or 9173 classifications per minute
Saved 1lng.mat
Loaded 1lnt.mat
Classifying interactions
Found   132 bases within 15 Angstroms from 3d structure
Found    49 pairs that are possibly interacting
Classification took 0.01 minutes, or 3360 classifications per minute
Saved 1lnt.mat
Read 1m5k.pdb for header information
Classifying interactions
Found  1810 bases within 15 Angstroms from 3d structure
Found   678 pairs that are possibly interacting
Classification took 0.06 minutes, or 11079 classifications per minute
Saved 1m5k.mat
Read 1m5o.pdb for header information
Classifying interactions
Found  1772 bases within 15 Angstroms from 3d structure
Bases C    3(D) and C   90(E) fall into categories   1.00  23.10 
Found   668 pairs that are possibly interacting
Classification took 0.06 minutes, or 10962 classifications per minute
Saved 1m5o.mat
Read 1m5p.pdb for header information
Classifying interactions
Found  1772 bases within 15 Angstroms from 3d structure
Found   678 pairs that are possibly interacting
Classification took 0.06 minutes, or 11032 classifications per minute
Saved 1m5p.mat
Read 1m5v.pdb for header information
Classifying interactions
Found  1777 bases within 15 Angstroms from 3d structure
Found   672 pairs that are possibly interacting
Classification took 0.06 minutes, or 10981 classifications per minute
Saved 1m5v.mat
Read 1m8v.pdb for header information
Classifying interactions
Found   112 bases within 15 Angstroms from 3d structure
Found    29 pairs that are possibly interacting
Classification took 0.01 minutes, or 2101 classifications per minute
Saved 1m8v.mat
Read 1m8x.pdb for header information
Classifying interactions
Found    24 bases within 15 Angstroms from 3d structure
Found     3 pairs that are possibly interacting
Classification took 0.01 minutes, or  250 classifications per minute
Saved 1m8x.mat
Read 1m8y.pdb for header information
Classifying interactions
Found    41 bases within 15 Angstroms from 3d structure
Found     8 pairs that are possibly interacting
Classification took 0.01 minutes, or  654 classifications per minute
Saved 1m8y.mat
Read 1mdg.pdb for header information
Classifying interactions
Found     7 bases within 15 Angstroms from 3d structure
Found     3 pairs that are possibly interacting
Classification took 0.01 minutes, or  256 classifications per minute
Saved 1mdg.mat
Read 1mfq.pdb for header information
Classifying interactions
Found  1016 bases within 15 Angstroms from 3d structure
Found   393 pairs that are possibly interacting
Classification took 0.04 minutes, or 10061 classifications per minute
Saved 1mfq.mat
Read 1mhk.pdb for header information
Classifying interactions
Found   155 bases within 15 Angstroms from 3d structure
Found    63 pairs that are possibly interacting
Classification took 0.02 minutes, or 4032 classifications per minute
Saved 1mhk.mat
Read 1mji.pdb for header information
Classifying interactions
Found   519 bases within 15 Angstroms from 3d structure
Found   206 pairs that are possibly interacting
Classification took 0.03 minutes, or 7832 classifications per minute
Saved 1mji.mat
Read 1mms.pdb for header information
Classifying interactions
Found  1094 bases within 15 Angstroms from 3d structure
Found   457 pairs that are possibly interacting
Classification took 0.05 minutes, or 10144 classifications per minute
Saved 1mms.mat
Read 1msw.pdb for header information
Classifying interactions
Found   171 bases within 15 Angstroms from 3d structure
Found    70 pairs that are possibly interacting
Classification took 0.02 minutes, or 4407 classifications per minute
Saved 1msw.mat
Read 1msy.pdb for header information
Classifying interactions
Found   174 bases within 15 Angstroms from 3d structure
Found    68 pairs that are possibly interacting
Classification took 0.02 minutes, or 4212 classifications per minute
Saved 1msy.mat
Read 1mzp.pdb for header information
Classifying interactions
Found   422 bases within 15 Angstroms from 3d structure
Found   187 pairs that are possibly interacting
Classification took 0.03 minutes, or 7480 classifications per minute
Saved 1mzp.mat
Read 1n1h.pdb for header information
Classifying interactions
Found     5 bases within 15 Angstroms from 3d structure
Found     2 pairs that are possibly interacting
Classification took 0.01 minutes, or  179 classifications per minute
Saved 1n1h.mat
Read 1n35.pdb for header information
Classifying interactions
Found    52 bases within 15 Angstroms from 3d structure
Found    24 pairs that are possibly interacting
Classification took 0.01 minutes, or 1807 classifications per minute
Saved 1n35.mat
Read 1n38.pdb for header information
Classifying interactions
Found    13 bases within 15 Angstroms from 3d structure
Found     9 pairs that are possibly interacting
Classification took 0.01 minutes, or  768 classifications per minute
Saved 1n38.mat
Read 1n77.pdb for header information
Classifying interactions
Found  1178 bases within 15 Angstroms from 3d structure
Found   439 pairs that are possibly interacting
Classification took 0.04 minutes, or 10406 classifications per minute
Saved 1n77.mat
Read 1n78.pdb for header information
Classifying interactions
Found  1176 bases within 15 Angstroms from 3d structure
Found   438 pairs that are possibly interacting
Classification took 0.04 minutes, or 10447 classifications per minute
Saved 1n78.mat
Read 1n7a.pdb for header information
Classifying interactions
Found   411 bases within 15 Angstroms from 3d structure
Found   201 pairs that are possibly interacting
Classification took 0.03 minutes, or 7494 classifications per minute
Saved 1n7a.mat
Read 1n7b.pdb for header information
Classifying interactions
Found   411 bases within 15 Angstroms from 3d structure
Found   201 pairs that are possibly interacting
Classification took 0.03 minutes, or 7494 classifications per minute
Saved 1n7b.mat
Read 1nbs.pdb for header information
Classifying interactions
Found  2361 bases within 15 Angstroms from 3d structure
Found   979 pairs that are possibly interacting
Classification took 0.09 minutes, or 11427 classifications per minute
Saved 1nbs.mat
Read 1nh3.pdb for header information
Classifying interactions
Found   247 bases within 15 Angstroms from 3d structure
Found    83 pairs that are possibly interacting
Classification took 0.02 minutes, or 4757 classifications per minute
Saved 1nh3.mat
Read 1nb7.pdb for header information
Classifying interactions
Found     8 bases within 15 Angstroms from 3d structure
Found     3 pairs that are possibly interacting
Classification took 0.01 minutes, or  268 classifications per minute
Saved 1nb7.mat
Read 1nta.pdb for header information
Classifying interactions
Found   317 bases within 15 Angstroms from 3d structure
Found   141 pairs that are possibly interacting
Classification took 0.02 minutes, or 6523 classifications per minute
Saved 1nta.mat
Read 1ntb.pdb for header information
Classifying interactions
Found   321 bases within 15 Angstroms from 3d structure
Found   139 pairs that are possibly interacting
Classification took 0.02 minutes, or 6431 classifications per minute
Saved 1ntb.mat
Read 1nyi.pdb for header information
Classifying interactions
Found   270 bases within 15 Angstroms from 3d structure
Found   108 pairs that are possibly interacting
Classification took 0.02 minutes, or 5841 classifications per minute
Saved 1nyi.mat
Read 1o0b.pdb for header information
Classifying interactions
Found   541 bases within 15 Angstroms from 3d structure
Found   209 pairs that are possibly interacting
Classification took 0.03 minutes, or 7946 classifications per minute
Saved 1o0b.mat
Read 1o0c.pdb for header information
Classifying interactions
Found   529 bases within 15 Angstroms from 3d structure
Found   203 pairs that are possibly interacting
Classification took 0.03 minutes, or 8036 classifications per minute
Saved 1o0c.mat
Read 1o9m.pdb for header information
Classifying interactions
Found   276 bases within 15 Angstroms from 3d structure
Found   104 pairs that are possibly interacting
Classification took 0.02 minutes, or 5625 classifications per minute
Saved 1o9m.mat
Read 1ob2.pdb for header information
Classifying interactions
Found   410 bases within 15 Angstroms from 3d structure
Found   147 pairs that are possibly interacting
Classification took 0.02 minutes, or 6720 classifications per minute
Saved 1ob2.mat
Read 1ofx.pdb for header information
Classifying interactions
Found    70 bases within 15 Angstroms from 3d structure
Found    29 pairs that are possibly interacting
Classification took 0.01 minutes, or 2142 classifications per minute
Saved 1ofx.mat
Read 1ooa.pdb for header information
Classifying interactions
Found   344 bases within 15 Angstroms from 3d structure
Found   142 pairs that are possibly interacting
Classification took 0.02 minutes, or 6902 classifications per minute
Saved 1ooa.mat
Read 1p6v.pdb for header information
Classifying interactions
Found   414 bases within 15 Angstroms from 3d structure
Found   171 pairs that are possibly interacting
Classification took 0.02 minutes, or 7061 classifications per minute
Saved 1p6v.mat
Read 1p79.pdb for header information
Classifying interactions
Found    10 bases within 15 Angstroms from 3d structure
Found     6 pairs that are possibly interacting
Classification took 0.01 minutes, or  501 classifications per minute
Saved 1p79.mat
Read 1pgl.pdb for header information
Classifying interactions
Found    14 bases within 15 Angstroms from 3d structure
Found     5 pairs that are possibly interacting
Classification took 0.01 minutes, or  427 classifications per minute
Saved 1pgl.mat
Read 1pvo.pdb for header information
Classifying interactions
Found     5 bases within 15 Angstroms from 3d structure
Found     2 pairs that are possibly interacting
Classification took 0.01 minutes, or  183 classifications per minute
Saved 1pvo.mat
Read 1pwf.pdb for header information
Classifying interactions
Found    49 bases within 15 Angstroms from 3d structure
Found    23 pairs that are possibly interacting
Classification took 0.01 minutes, or 1802 classifications per minute
Saved 1pwf.mat
Read 1q2r.pdb for header information
Classifying interactions
Found   225 bases within 15 Angstroms from 3d structure
Found    85 pairs that are possibly interacting
Classification took 0.02 minutes, or 4800 classifications per minute
Saved 1q2r.mat
Read 1q2s.pdb for header information
Classifying interactions
Found   221 bases within 15 Angstroms from 3d structure
Found    93 pairs that are possibly interacting
Classification took 0.02 minutes, or 5030 classifications per minute
Saved 1q2s.mat
Read 1q93.pdb for header information
Classifying interactions
Found   599 bases within 15 Angstroms from 3d structure
Found   236 pairs that are possibly interacting
Classification took 0.03 minutes, or 8239 classifications per minute
Saved 1q93.mat
Read 1qa6.pdb for header information
Classifying interactions
Found  1040 bases within 15 Angstroms from 3d structure
Found   436 pairs that are possibly interacting
Classification took 0.04 minutes, or 9848 classifications per minute
Saved 1qa6.mat
Read 1qbp.pdb for header information
Classifying interactions
Found   536 bases within 15 Angstroms from 3d structure
Found   164 pairs that are possibly interacting
Classification took 0.02 minutes, or 7323 classifications per minute
Saved 1qbp.mat
Read 1qc0.pdb for header information
Classifying interactions
Found   423 bases within 15 Angstroms from 3d structure
Found   146 pairs that are possibly interacting
Classification took 0.02 minutes, or 6837 classifications per minute
Saved 1qc0.mat
Read 1qf6.pdb for header information
Classifying interactions
Found   509 bases within 15 Angstroms from 3d structure
Found   189 pairs that are possibly interacting
Classification took 0.03 minutes, or 7560 classifications per minute
Saved 1qf6.mat
Read 1qln.pdb for header information
Classifying interactions
Found   103 bases within 15 Angstroms from 3d structure
Found    36 pairs that are possibly interacting
Classification took 0.01 minutes, or 2658 classifications per minute
Saved 1qln.mat
Read 1qrs.pdb for header information
Classifying interactions
Found   532 bases within 15 Angstroms from 3d structure
Found   218 pairs that are possibly interacting
Classification took 0.03 minutes, or 8127 classifications per minute
Saved 1qrs.mat
Read 1qrt.pdb for header information
Classifying interactions
Found   529 bases within 15 Angstroms from 3d structure
Found   215 pairs that are possibly interacting
Classification took 0.03 minutes, or 7938 classifications per minute
Saved 1qrt.mat
Read 1qru.pdb for header information
Classifying interactions
Found   534 bases within 15 Angstroms from 3d structure
Found   216 pairs that are possibly interacting
Classification took 0.03 minutes, or 7975 classifications per minute
Saved 1qru.mat
Read 1qtq.pdb for header information
Classifying interactions
Found   533 bases within 15 Angstroms from 3d structure
Found   208 pairs that are possibly interacting
Classification took 0.03 minutes, or 8068 classifications per minute
Saved 1qtq.mat
Read 1qu2.pdb for header information
Classifying interactions
Found   564 bases within 15 Angstroms from 3d structure
Found   210 pairs that are possibly interacting
Classification took 0.03 minutes, or 7906 classifications per minute
Saved 1qu2.mat
Read 1qu3.pdb for header information
Classifying interactions
Found   569 bases within 15 Angstroms from 3d structure
Found   207 pairs that are possibly interacting
Classification took 0.03 minutes, or 8029 classifications per minute
Saved 1qu3.mat
Read 1qzw.pdb for header information
Classifying interactions
Found  1164 bases within 15 Angstroms from 3d structure
Found   456 pairs that are possibly interacting
Classification took 0.04 minutes, or 10612 classifications per minute
Saved 1qzw.mat
Read 1r3e.pdb for header information
Classifying interactions
Found   363 bases within 15 Angstroms from 3d structure
Found   137 pairs that are possibly interacting
Classification took 0.02 minutes, or 6416 classifications per minute
Saved 1r3e.mat
Read 1r3g.pdb for header information
Classifying interactions
Found    75 bases within 15 Angstroms from 3d structure
Found    29 pairs that are possibly interacting
Classification took 0.01 minutes, or 2184 classifications per minute
Saved 1r3g.mat
Read 1r9f.pdb for header information
Classifying interactions
Found   262 bases within 15 Angstroms from 3d structure
Found    93 pairs that are possibly interacting
Classification took 0.02 minutes, or 5176 classifications per minute
Saved 1r9f.mat
Read 1r9s.pdb for header information
Classifying interactions
Found    80 bases within 15 Angstroms from 3d structure
Found    31 pairs that are possibly interacting
Classification took 0.01 minutes, or 2246 classifications per minute
Saved 1r9s.mat
Read 1r9t.pdb for header information
Classifying interactions
Found   188 bases within 15 Angstroms from 3d structure
Found    65 pairs that are possibly interacting
Classification took 0.02 minutes, or 4160 classifications per minute
Saved 1r9t.mat
Read 1rc7.pdb for header information
Classifying interactions
Found   276 bases within 15 Angstroms from 3d structure
Found    99 pairs that are possibly interacting
Classification took 0.02 minutes, or 5510 classifications per minute
Saved 1rc7.mat
Read 1rlg.pdb for header information
Classifying interactions
Found   234 bases within 15 Angstroms from 3d structure
Found   111 pairs that are possibly interacting
Classification took 0.02 minutes, or 5683 classifications per minute
Saved 1rlg.mat
Read 1rmv.pdb for header information
Classifying interactions
Found     3 bases within 15 Angstroms from 3d structure
Found     2 pairs that are possibly interacting
Classification took 0.01 minutes, or  171 classifications per minute
Saved 1rmv.mat
Read 1rna.pdb for header information
Classifying interactions
Found   205 bases within 15 Angstroms from 3d structure
Found    73 pairs that are possibly interacting
Classification took 0.02 minutes, or 4313 classifications per minute
Saved 1rna.mat
Read 1rxa.pdb for header information
Classifying interactions
Found    20 bases within 15 Angstroms from 3d structure
Found     7 pairs that are possibly interacting
Classification took 0.01 minutes, or  584 classifications per minute
Saved 1rxa.mat
Read 1rxb.pdb for header information
Classifying interactions
Found    89 bases within 15 Angstroms from 3d structure
Found    39 pairs that are possibly interacting
Classification took 0.01 minutes, or 2773 classifications per minute
Saved 1rxb.mat
Read 1s03.pdb for header information
Classifying interactions
Found   692 bases within 15 Angstroms from 3d structure
Found   272 pairs that are possibly interacting
Classification took 0.03 minutes, or 8777 classifications per minute
Saved 1s03.mat
Read 1s0v.pdb for header information
Classifying interactions
Found   542 bases within 15 Angstroms from 3d structure
Found   214 pairs that are possibly interacting
Classification took 0.03 minutes, or 8385 classifications per minute
Saved 1s0v.mat
Read 1s76.pdb for header information
Classifying interactions
Found   189 bases within 15 Angstroms from 3d structure
Found    78 pairs that are possibly interacting
Classification took 0.02 minutes, or 4538 classifications per minute
Saved 1s76.mat
Read 1s77.pdb for header information
Classifying interactions
Found   188 bases within 15 Angstroms from 3d structure
Found    75 pairs that are possibly interacting
Classification took 0.02 minutes, or 4431 classifications per minute
Saved 1s77.mat
Read 1si2.pdb for header information
Classifying interactions
Found    19 bases within 15 Angstroms from 3d structure
Found     7 pairs that are possibly interacting
Classification took 0.01 minutes, or  584 classifications per minute
Saved 1si2.mat
Read 1si3.pdb for header information
Classifying interactions
Found    22 bases within 15 Angstroms from 3d structure
Found    10 pairs that are possibly interacting
Classification took 0.01 minutes, or  817 classifications per minute
Saved 1si3.mat
Read 1sz1.pdb for header information
Classifying interactions
Found   820 bases within 15 Angstroms from 3d structure
Found   297 pairs that are possibly interacting
Classification took 0.03 minutes, or 9051 classifications per minute
Saved 1sz1.mat
Read 1t0d.pdb for header information
Classifying interactions
Found   456 bases within 15 Angstroms from 3d structure
Found   168 pairs that are possibly interacting
Classification took 0.02 minutes, or 7168 classifications per minute
Saved 1t0d.mat
Read 1t0k.pdb for header information
Classifying interactions
Found   170 bases within 15 Angstroms from 3d structure
Found    83 pairs that are possibly interacting
Classification took 0.02 minutes, or 4980 classifications per minute
Saved 1t0k.mat
Read 1tfw.pdb for header information
Classifying interactions
Found   489 bases within 15 Angstroms from 3d structure
Found   178 pairs that are possibly interacting
Classification took 0.02 minutes, or 7430 classifications per minute
Saved 1tfw.mat
Read 1tfy.pdb for header information
Classifying interactions
Found   466 bases within 15 Angstroms from 3d structure
Found   179 pairs that are possibly interacting
Classification took 0.02 minutes, or 7553 classifications per minute
Saved 1tfy.mat
Read 1ttt.pdb for header information
Classifying interactions
Found  1354 bases within 15 Angstroms from 3d structure
Found   506 pairs that are possibly interacting
Classification took 0.05 minutes, or 10503 classifications per minute
Saved 1ttt.mat
Read 1u0b.pdb for header information
Classifying interactions
Found   575 bases within 15 Angstroms from 3d structure
Found   220 pairs that are possibly interacting
Classification took 0.03 minutes, or 7970 classifications per minute
Saved 1u0b.mat
Read 1u1y.pdb for header information
Classifying interactions
Found   164 bases within 15 Angstroms from 3d structure
Found    64 pairs that are possibly interacting
Classification took 0.02 minutes, or 4029 classifications per minute
Saved 1u1y.mat
Read 1u63.pdb for header information
Classifying interactions
Found   723 bases within 15 Angstroms from 3d structure
Found   335 pairs that are possibly interacting
Classification took 0.04 minutes, or 9189 classifications per minute
Saved 1u63.mat
Read 1u6b.pdb for header information
Classifying interactions
Found  1972 bases within 15 Angstroms from 3d structure
Found   779 pairs that are possibly interacting
Classification took 0.07 minutes, or 11162 classifications per minute
Saved 1u6b.mat
Read 1u8d.pdb for header information
Classifying interactions
Found   601 bases within 15 Angstroms from 3d structure
Found   232 pairs that are possibly interacting
Classification took 0.03 minutes, or 8099 classifications per minute
Saved 1u8d.mat
Read 1u9s.pdb for header information
Classifying interactions
Found  1370 bases within 15 Angstroms from 3d structure
Found   538 pairs that are possibly interacting
Classification took 0.05 minutes, or 10594 classifications per minute
Saved 1u9s.mat
Read 1un6.pdb for header information
Classifying interactions
Found   889 bases within 15 Angstroms from 3d structure
Found   338 pairs that are possibly interacting
Classification took 0.04 minutes, or 9544 classifications per minute
Saved 1un6.mat
Read 1urn.pdb for header information
Classifying interactions
Found   227 bases within 15 Angstroms from 3d structure
Found   101 pairs that are possibly interacting
Classification took 0.02 minutes, or 5387 classifications per minute
Saved 1urn.mat
Read 1utd.pdb for header information
Classifying interactions
Found    41 bases within 15 Angstroms from 3d structure
Found    26 pairs that are possibly interacting
Classification took 0.01 minutes, or 1884 classifications per minute
Saved 1utd.mat
Read 1utf.pdb for header information
Classifying interactions
Found    19 bases within 15 Angstroms from 3d structure
Found    13 pairs that are possibly interacting
Classification took 0.01 minutes, or 1040 classifications per minute
Saved 1utf.mat
Read 1utv.pdb for header information
Classifying interactions
Found    15 bases within 15 Angstroms from 3d structure
Found    12 pairs that are possibly interacting
Classification took 0.01 minutes, or 1002 classifications per minute
Saved 1utv.mat
Read 1uvi.pdb for header information
Classifying interactions
Found    12 bases within 15 Angstroms from 3d structure
Found     6 pairs that are possibly interacting
Classification took 0.01 minutes, or  512 classifications per minute
Saved 1uvi.mat
Read 1uvj.pdb for header information
Classifying interactions
Found    15 bases within 15 Angstroms from 3d structure
Found     6 pairs that are possibly interacting
Classification took 0.01 minutes, or  501 classifications per minute
Saved 1uvj.mat
Read 1uvm.pdb for header information
Classifying interactions
Found    21 bases within 15 Angstroms from 3d structure
Found    18 pairs that are possibly interacting
Classification took 0.01 minutes, or 1471 classifications per minute
Saved 1uvm.mat
Read 1vfg.pdb for header information
Classifying interactions
Found   431 bases within 15 Angstroms from 3d structure
Found   188 pairs that are possibly interacting
Classification took 0.03 minutes, or 7520 classifications per minute
Saved 1vfg.mat
Read 1wmq.pdb for header information
Classifying interactions
Found    20 bases within 15 Angstroms from 3d structure
Found     8 pairs that are possibly interacting
Classification took 0.01 minutes, or  683 classifications per minute
Saved 1wmq.mat
Read 1wsu.pdb for header information
Classifying interactions
Found   451 bases within 15 Angstroms from 3d structure
Found   174 pairs that are possibly interacting
Classification took 0.02 minutes, or 7424 classifications per minute
Saved 1wsu.mat
Read 1x8w.pdb for header information
Classifying interactions
Found  9379 bases within 15 Angstroms from 3d structure
Found  3818 pairs that are possibly interacting
Classification took 0.31 minutes, or 12446 classifications per minute
Saved 1x8w.mat
Read 1xjr.pdb for header information
Classifying interactions
Found   357 bases within 15 Angstroms from 3d structure
Found   144 pairs that are possibly interacting
Classification took 0.02 minutes, or 6505 classifications per minute
Saved 1xjr.mat
Read 1xok.pdb for header information
Classifying interactions
Found   183 bases within 15 Angstroms from 3d structure
Found    81 pairs that are possibly interacting
Classification took 0.02 minutes, or 4860 classifications per minute
Saved 1xok.mat
Read 1xpo.pdb for header information
Classifying interactions
Found     6 bases within 15 Angstroms from 3d structure
Found     4 pairs that are possibly interacting
Classification took 0.01 minutes, or  341 classifications per minute
Saved 1xpo.mat
Read 1xpr.pdb for header information
Classifying interactions
Found     5 bases within 15 Angstroms from 3d structure
Found     4 pairs that are possibly interacting
Classification took 0.01 minutes, or  341 classifications per minute
Saved 1xpr.mat
Read 1xpu.pdb for header information
Classifying interactions
Found     6 bases within 15 Angstroms from 3d structure
Found     2 pairs that are possibly interacting
Classification took 0.01 minutes, or  171 classifications per minute
Saved 1xpu.mat
Read 1y0q.pdb for header information
Classifying interactions
Found  2173 bases within 15 Angstroms from 3d structure
Found   927 pairs that are possibly interacting
Classification took 0.08 minutes, or 11409 classifications per minute
Saved 1y0q.mat
Read 1y1w.pdb for header information
Classifying interactions
Found   123 bases within 15 Angstroms from 3d structure
Found    49 pairs that are possibly interacting
Classification took 0.01 minutes, or 3301 classifications per minute
Saved 1y1w.mat
Read 1y26.pdb for header information
Classifying interactions
Found   618 bases within 15 Angstroms from 3d structure
Found   238 pairs that are possibly interacting
Classification took 0.03 minutes, or 8308 classifications per minute
Saved 1y26.mat
Read 1y27.pdb for header information
Classifying interactions
Found   594 bases within 15 Angstroms from 3d structure
Found   228 pairs that are possibly interacting
Classification took 0.03 minutes, or 8260 classifications per minute
Saved 1y27.mat
Read 1y39.pdb for header information
Classifying interactions
Found  1039 bases within 15 Angstroms from 3d structure
Found   452 pairs that are possibly interacting
Classification took 0.04 minutes, or 10091 classifications per minute
Saved 1y39.mat
Read 1y77.pdb for header information
Classifying interactions
Found   121 bases within 15 Angstroms from 3d structure
Found    48 pairs that are possibly interacting
Classification took 0.01 minutes, or 3234 classifications per minute
Saved 1y77.mat
Read 1ykq.pdb for header information
Classifying interactions
Found   782 bases within 15 Angstroms from 3d structure
Found   335 pairs that are possibly interacting
Classification took 0.04 minutes, or 9390 classifications per minute
Saved 1ykq.mat
Read 1ykv.pdb for header information
Classifying interactions
Found   779 bases within 15 Angstroms from 3d structure
Found   323 pairs that are possibly interacting
Classification took 0.03 minutes, or 9396 classifications per minute
Saved 1ykv.mat
Read 1yls.pdb for header information
Classifying interactions
Found   651 bases within 15 Angstroms from 3d structure
Found   268 pairs that are possibly interacting
Classification took 0.03 minutes, or 8796 classifications per minute
Saved 1yls.mat
Read 1yrj.pdb for header information
Classifying interactions
Found   286 bases within 15 Angstroms from 3d structure
Found   110 pairs that are possibly interacting
Classification took 0.02 minutes, or 5786 classifications per minute
Saved 1yrj.mat
Read 1ytu.pdb for header information
Classifying interactions
Found    67 bases within 15 Angstroms from 3d structure
Found    37 pairs that are possibly interacting
Classification took 0.01 minutes, or 2732 classifications per minute
Saved 1ytu.mat
Read 1yvp.pdb for header information
Classifying interactions
Found   155 bases within 15 Angstroms from 3d structure
Found    67 pairs that are possibly interacting
Classification took 0.02 minutes, or 4218 classifications per minute
Saved 1yvp.mat
Read 1zdh.pdb for header information
Classifying interactions
Found   111 bases within 15 Angstroms from 3d structure
Found    44 pairs that are possibly interacting
Classification took 0.01 minutes, or 3072 classifications per minute
Saved 1zdh.mat
Read 1zdi.pdb for header information
Classifying interactions
Found   136 bases within 15 Angstroms from 3d structure
Found    51 pairs that are possibly interacting
Classification took 0.01 minutes, or 3436 classifications per minute
Saved 1zdi.mat
Read 1zdj.pdb for header information
Classifying interactions
Found    52 bases within 15 Angstroms from 3d structure
Found    24 pairs that are possibly interacting
Classification took 0.01 minutes, or 1807 classifications per minute
Saved 1zdj.mat
Read 1zdk.pdb for header information
Classifying interactions
Found   211 bases within 15 Angstroms from 3d structure
Found    75 pairs that are possibly interacting
Classification took 0.02 minutes, or 4431 classifications per minute
Saved 1zdk.mat
Read 1ze2.pdb for header information
Classifying interactions
Found   277 bases within 15 Angstroms from 3d structure
Found   112 pairs that are possibly interacting
Classification took 0.02 minutes, or 5892 classifications per minute
Saved 1ze2.mat
Read 1zjw.pdb for header information
Classifying interactions
Found   535 bases within 15 Angstroms from 3d structure
Found   208 pairs that are possibly interacting
Classification took 0.03 minutes, or 8150 classifications per minute
Saved 1zjw.mat
Read 246d.pdb for header information
Classifying interactions
Found    97 bases within 15 Angstroms from 3d structure
Found    38 pairs that are possibly interacting
Classification took 0.01 minutes, or 2702 classifications per minute
Saved 246d.mat
Read 259d.pdb for header information
Classifying interactions
Found    90 bases within 15 Angstroms from 3d structure
Found    39 pairs that are possibly interacting
Classification took 0.01 minutes, or 2880 classifications per minute
Saved 259d.mat
Read 280d.pdb for header information
Classifying interactions
Found   337 bases within 15 Angstroms from 3d structure
Found   120 pairs that are possibly interacting
Classification took 0.02 minutes, or 6227 classifications per minute
Saved 280d.mat
Read 283d.pdb for header information
Classifying interactions
Found    34 bases within 15 Angstroms from 3d structure
Found    11 pairs that are possibly interacting
Classification took 0.01 minutes, or  899 classifications per minute
Saved 283d.mat
Read 299d.pdb for header information
Classifying interactions
Found   285 bases within 15 Angstroms from 3d structure
Found   112 pairs that are possibly interacting
Classification took 0.02 minutes, or 5892 classifications per minute
Saved 299d.mat
Read 2a8v.pdb for header information
Classifying interactions
Found     8 bases within 15 Angstroms from 3d structure
Found     2 pairs that are possibly interacting
Classification took 0.01 minutes, or  175 classifications per minute
Saved 2a8v.mat
Read 2bbv.pdb for header information
Classifying interactions
Found    26 bases within 15 Angstroms from 3d structure
Found    10 pairs that are possibly interacting
Classification took 0.01 minutes, or  873 classifications per minute
Saved 2bbv.mat
Read 2bgg.pdb for header information
Classifying interactions
Found   130 bases within 15 Angstroms from 3d structure
Found    57 pairs that are possibly interacting
Classification took 0.02 minutes, or 3710 classifications per minute
Saved 2bgg.mat
Read 2bh2.pdb for header information
Classifying interactions
Found   269 bases within 15 Angstroms from 3d structure
Found   137 pairs that are possibly interacting
Classification took 0.02 minutes, or 6416 classifications per minute
Saved 2bh2.mat
Read 2bj6.pdb for header information
Classifying interactions
Found   114 bases within 15 Angstroms from 3d structure
Found    37 pairs that are possibly interacting
Classification took 0.01 minutes, or 2631 classifications per minute
Saved 2bj6.mat
Read 2fmt.pdb for header information
Classifying interactions
Found  1032 bases within 15 Angstroms from 3d structure
Found   381 pairs that are possibly interacting
Classification took 0.04 minutes, or 9885 classifications per minute
Saved 2fmt.mat
Read 300d.pdb for header information
Classifying interactions
Found   287 bases within 15 Angstroms from 3d structure
Found   114 pairs that are possibly interacting
Classification took 0.02 minutes, or 5837 classifications per minute
Saved 300d.mat
Read 301d.pdb for header information
Classifying interactions
Found   285 bases within 15 Angstroms from 3d structure
Found   112 pairs that are possibly interacting
Classification took 0.02 minutes, or 5812 classifications per minute
Saved 301d.mat
Read 315d.pdb for header information
Classifying interactions
Found    88 bases within 15 Angstroms from 3d structure
Found    37 pairs that are possibly interacting
Classification took 0.01 minutes, or 2631 classifications per minute
Saved 315d.mat
Read 332d.pdb for header information
Classifying interactions
Found    85 bases within 15 Angstroms from 3d structure
Found    36 pairs that are possibly interacting
Classification took 0.01 minutes, or 2658 classifications per minute
Saved 332d.mat
Read 333d.pdb for header information
Classifying interactions
Found    15 bases within 15 Angstroms from 3d structure
Found     5 pairs that are possibly interacting
Classification took 0.01 minutes, or  417 classifications per minute
Saved 333d.mat
Read 353d.pdb for header information
Classifying interactions
Found   164 bases within 15 Angstroms from 3d structure
Found    60 pairs that are possibly interacting
Classification took 0.02 minutes, or 3840 classifications per minute
Saved 353d.mat
Read 354d.pdb for header information
Classifying interactions
Found   139 bases within 15 Angstroms from 3d structure
Found    54 pairs that are possibly interacting
Classification took 0.02 minutes, or 3575 classifications per minute
Saved 354d.mat
Read 357d.pdb for header information
Classifying interactions
Found   404 bases within 15 Angstroms from 3d structure
Found   153 pairs that are possibly interacting
Classification took 0.02 minutes, or 6912 classifications per minute
Saved 357d.mat
Read 359d.pdb for header information
Classifying interactions
Found   280 bases within 15 Angstroms from 3d structure
Found   112 pairs that are possibly interacting
Classification took 0.02 minutes, or 5812 classifications per minute
Saved 359d.mat
Read 361d.pdb for header information
Classifying interactions
Found   269 bases within 15 Angstroms from 3d structure
Found   103 pairs that are possibly interacting
Classification took 0.02 minutes, or 5650 classifications per minute
Saved 361d.mat
Read 364d.pdb for header information
Classifying interactions
Found   391 bases within 15 Angstroms from 3d structure
Found   148 pairs that are possibly interacting
Classification took 0.02 minutes, or 7016 classifications per minute
Saved 364d.mat
Read 377d.pdb for header information
Classifying interactions
Found   127 bases within 15 Angstroms from 3d structure
Found    58 pairs that are possibly interacting
Classification took 0.02 minutes, or 3712 classifications per minute
Saved 377d.mat
Read 379d.pdb for header information
Classifying interactions
Found   283 bases within 15 Angstroms from 3d structure
Found   118 pairs that are possibly interacting
Classification took 0.02 minutes, or 6123 classifications per minute
Saved 379d.mat
Read 387d.pdb for header information
Classifying interactions
Found   122 bases within 15 Angstroms from 3d structure
Found    46 pairs that are possibly interacting
Classification took 0.01 minutes, or 3099 classifications per minute
Saved 387d.mat
Read 397d.pdb for header information
Classifying interactions
Found   177 bases within 15 Angstroms from 3d structure
Found    70 pairs that are possibly interacting
Classification took 0.02 minutes, or 4267 classifications per minute
Saved 397d.mat
Read 398d.pdb for header information
Classifying interactions
Found   193 bases within 15 Angstroms from 3d structure
Found    72 pairs that are possibly interacting
Classification took 0.02 minutes, or 4320 classifications per minute
Saved 398d.mat
Read 402d.pdb for header information
Classifying interactions
Found    92 bases within 15 Angstroms from 3d structure
Found    37 pairs that are possibly interacting
Classification took 0.01 minutes, or 2732 classifications per minute
Saved 402d.mat
Read 404d.pdb for header information
Classifying interactions
Found    61 bases within 15 Angstroms from 3d structure
Found    23 pairs that are possibly interacting
Classification took 0.01 minutes, or 1732 classifications per minute
Saved 404d.mat
Read 405d.pdb for header information
Classifying interactions
Found   224 bases within 15 Angstroms from 3d structure
Found    77 pairs that are possibly interacting
Classification took 0.02 minutes, or 4693 classifications per minute
Saved 405d.mat
Read 409d.pdb for header information
Classifying interactions
Found   210 bases within 15 Angstroms from 3d structure
Found   100 pairs that are possibly interacting
Classification took 0.02 minutes, or 5565 classifications per minute
Saved 409d.mat
Read 418d.pdb for header information
Classifying interactions
Found    86 bases within 15 Angstroms from 3d structure
Found    36 pairs that are possibly interacting
Classification took 0.01 minutes, or 2658 classifications per minute
Saved 418d.mat
Read 419d.pdb for header information
Classifying interactions
Found   205 bases within 15 Angstroms from 3d structure
Found   103 pairs that are possibly interacting
Classification took 0.02 minutes, or 5493 classifications per minute
Saved 419d.mat
Read 420d.pdb for header information
Classifying interactions
Found   229 bases within 15 Angstroms from 3d structure
Found    78 pairs that are possibly interacting
Classification took 0.02 minutes, or 4608 classifications per minute
Saved 420d.mat
Read 421d.pdb for header information
Classifying interactions
Found    38 bases within 15 Angstroms from 3d structure
Found    14 pairs that are possibly interacting
Classification took 0.01 minutes, or 1169 classifications per minute
Saved 421d.mat
Read 422d.pdb for header information
Classifying interactions
Found   132 bases within 15 Angstroms from 3d structure
Found    47 pairs that are possibly interacting
Classification took 0.01 minutes, or 3223 classifications per minute
Saved 422d.mat
Read 429d.pdb for header information
Classifying interactions
Found   305 bases within 15 Angstroms from 3d structure
Found   124 pairs that are possibly interacting
Classification took 0.02 minutes, or 6027 classifications per minute
Saved 429d.mat
Read 430d.pdb for header information
Classifying interactions
Found   180 bases within 15 Angstroms from 3d structure
Found    68 pairs that are possibly interacting
Classification took 0.02 minutes, or 4352 classifications per minute
Saved 430d.mat
Read 433d.pdb for header information
Classifying interactions
Found   187 bases within 15 Angstroms from 3d structure
Found    71 pairs that are possibly interacting
Classification took 0.02 minutes, or 4397 classifications per minute
Saved 433d.mat
Read 434d.pdb for header information
Classifying interactions
Found   180 bases within 15 Angstroms from 3d structure
Found    67 pairs that are possibly interacting
Classification took 0.02 minutes, or 4361 classifications per minute
Saved 434d.mat
Read 435d.pdb for header information
Classifying interactions
Found   185 bases within 15 Angstroms from 3d structure
Found    66 pairs that are possibly interacting
Classification took 0.02 minutes, or 4155 classifications per minute
Saved 435d.mat
Read 438d.pdb for header information
Classifying interactions
Found    94 bases within 15 Angstroms from 3d structure
Found    39 pairs that are possibly interacting
Classification took 0.01 minutes, or 2773 classifications per minute
Saved 438d.mat
Read 439d.pdb for header information
Classifying interactions
Found    90 bases within 15 Angstroms from 3d structure
Found    38 pairs that are possibly interacting
Classification took 0.01 minutes, or 2702 classifications per minute
Saved 439d.mat
Read 464d.pdb for header information
Classifying interactions
Found   179 bases within 15 Angstroms from 3d structure
Found    67 pairs that are possibly interacting
Classification took 0.02 minutes, or 4288 classifications per minute
Saved 464d.mat
Read 466d.pdb for header information
Classifying interactions
Found   182 bases within 15 Angstroms from 3d structure
Found    67 pairs that are possibly interacting
Classification took 0.02 minutes, or 4218 classifications per minute
Saved 466d.mat
Read 472d.pdb for header information
Classifying interactions
Found    98 bases within 15 Angstroms from 3d structure
Found    40 pairs that are possibly interacting
Classification took 0.01 minutes, or 2898 classifications per minute
Saved 472d.mat
Read 479d.pdb for header information
Classifying interactions
Found    51 bases within 15 Angstroms from 3d structure
Found    18 pairs that are possibly interacting
Classification took 0.01 minutes, or 1411 classifications per minute
Saved 479d.mat
Read 480d.pdb for header information
Classifying interactions
Found   169 bases within 15 Angstroms from 3d structure
Found    67 pairs that are possibly interacting
Classification took 0.02 minutes, or 4150 classifications per minute
Saved 480d.mat
Read 485d.pdb for header information
Classifying interactions
Found    89 bases within 15 Angstroms from 3d structure
Found    38 pairs that are possibly interacting
Classification took 0.01 minutes, or 2753 classifications per minute
Saved 485d.mat
Read 488d.pdb for header information
Classifying interactions
Found   706 bases within 15 Angstroms from 3d structure
Found   297 pairs that are possibly interacting
Classification took 0.03 minutes, or 9197 classifications per minute
Saved 488d.mat
Read 5msf.pdb for header information
Classifying interactions
Found   175 bases within 15 Angstroms from 3d structure
Found    64 pairs that are possibly interacting
Classification took 0.02 minutes, or 4165 classifications per minute
Saved 5msf.mat
Read 6msf.pdb for header information
Classifying interactions
Found    99 bases within 15 Angstroms from 3d structure
Found    46 pairs that are possibly interacting
Classification took 0.01 minutes, or 3099 classifications per minute
Saved 6msf.mat
Could not open file 7msf    .pdb
??? Output argument "Files" (and maybe others) not assigned during call to "C:\Documents and Settings\zirbel\My Documents\FR3D\zGetNTData.m (zGetNTData)".

Error in ==> zGetNTData at 16
if nargin < 2,

Error in ==> zAddNTData at 46
    File(F+1) = zGetNTData(FullList{f},ReadCode); %   load it

>> File

File = 

1x4 struct array with fields:
    CI
    ClassVersion
    Comment
    Distance
    Edge
    Filename
    HandClass
    Header
    Inter
    Modified
    NT
    NumNT
    Pair

>> zGetNTData('7msf   ')
Could not open file 7msf   .pdb
>> zGetNTData('7msf   ')
>> File = zAddNTData('HighResolution_list',0,File)
Read HighResolution_list.pdb
Loaded 1q96  .mat
Loaded 157d.mat
Loaded 165d.mat
Loaded 1a9n.mat
Loaded 1bmv.mat
Loaded 1d96.mat
Loaded 1dk1.mat
Loaded 1dpl.mat
Loaded 1fir.mat
Loaded 1hmh.mat
Loaded 1i5l.mat
Loaded 1i2y.mat
Loaded 1id9.mat
Loaded 1idw.mat
Loaded 1iha.mat
Loaded 1kuq.mat
Loaded 1l2x.mat
Loaded 1l3d.mat
Error in ==> zGetNTData at 46
      load(strcat(Filename,'.mat'),'File','-mat');

Error in ==> zAddNTData at 46
    File(F+1) = zGetNTData(FullList{f},ReadCode); %   load it


>> File = zAddNTData('HighResolution_list',0,File)
Read HighResolution_list.pdb
Loaded 1q96.mat
Loaded 157d.mat
Loaded 165d.mat
Loaded 1a9n.mat
Loaded 1bmv.mat
Loaded 1d96.mat
Loaded 1dk1.mat
Loaded 1dpl.mat
Loaded 1fir.mat
Loaded 1hmh.mat
Loaded 1i5l.mat
Loaded 1i2y.mat
Loaded 1id9.mat
Loaded 1idw.mat
Loaded 1iha.mat
Loaded 1kuq.mat
Loaded 1l2x.mat
Loaded 1l3d.mat
Loaded 1mme.mat
Loaded 1mwl.mat
Loaded 1nlc.mat
Loaded 1nuj.mat
Loaded 1nuv.mat
Loaded 1o3z.mat
Loaded 1osu.mat
Loaded 1q9a.mat
Loaded 1rpu.mat
Loaded 1sdr.mat
Loaded 1ser.mat
Loaded 1t03.mat
Loaded 1t0e.mat
Loaded 1tn2.mat
Loaded 1tra.mat
Loaded 1wvd.mat
Loaded 1y6s.mat
Loaded 1y6t.mat
Loaded 1y73.mat
Loaded 1y90.mat
Loaded 1y95.mat
Loaded 1y99.mat
Loaded 205d.mat
Loaded 247d.mat
Loaded 248d.mat
Loaded 255d.mat
Loaded 2tra.mat
Loaded 3tra.mat
Loaded 413d.mat
Loaded 437d.mat
Loaded 462d.mat
Loaded 483d.mat
Loaded 4tna.mat
Loaded 4tra.mat
Loaded 6tna.mat
Loaded 1a34.mat
Loaded 1apg.mat
Loaded 1aq3.mat
Loaded 1aq4.mat
Loaded 1asy.mat
Loaded 1asz.mat
Loaded 1av6.mat
Loaded 1b23.mat
Loaded 1b7f.mat
Loaded 1br3.mat
Loaded 1by4.mat
Loaded 1c04.mat
Loaded 1c0a.mat
Loaded 1csl.mat
Loaded 1cvj.mat
Loaded 1cwp.mat
Loaded 1cx0.mat
Loaded 1d4r.mat
Loaded 1d9h.mat
Loaded 1ddl.mat
Loaded 1ddy.mat
Loaded 1dfu.mat
Loaded 1di2.mat
Loaded 1dno.mat
Loaded 1dnt.mat
Loaded 1dnx.mat
Loaded 1dqf.mat
Loaded 1dqh.mat
Loaded 1drz.mat
Loaded 1duh.mat
Loaded 1dul.mat
Loaded 1duq.mat
Loaded 1dzs.mat
Loaded 1e6t.mat
Loaded 1e7k.mat
Loaded 1e7x.mat
Loaded 1e8o.mat
Loaded 1ec6.mat
Loaded 1efo.mat
Loaded 1efw.mat
Loaded 1egk.mat
Loaded 1ehz.mat
Loaded 1eiy.mat
Loaded 1et4.mat
Loaded 1euq.mat
Loaded 1euy.mat
Loaded 1evp.mat
Loaded 1evv.mat
Loaded 1exd.mat
Loaded 1f1t.mat
Loaded 1f27.mat
Loaded 1f7u.mat
Loaded 1f7v.mat
Loaded 1f8v.mat
Loaded 1feu.mat
Loaded 1fix.mat
Loaded 1fuf.mat
Loaded 1fxl.mat
Loaded 1g1x.mat
Loaded 1g2e.mat
Loaded 1g2j.mat
Loaded 1g4q.mat
Loaded 1g59.mat
Loaded 1gax.mat
Loaded 1gid.mat
Loaded 1gkv.mat
Loaded 1gkw.mat
Loaded 1grz.mat
Loaded 1gtf.mat
Loaded 1gtn.mat
Loaded 1gtr.mat
Loaded 1gts.mat
Loaded 1h1k.mat
Loaded 1h2d.mat
Loaded 1h38.mat
Loaded 1h3e.mat
Loaded 1h4q.mat
Loaded 1h4s.mat
Loaded 1h8j.mat
Loaded 1hc8.mat
Loaded 1hdw.mat
Loaded 1he0.mat
Loaded 1he6.mat
Loaded 1hq1.mat
Loaded 1hr2.mat
Loaded 1hvu.mat
Loaded 1hys.mat
Loaded 1i0f.mat
Loaded 1i0g.mat
Loaded 1i0j.mat
Loaded 1i0k.mat
Loaded 1i0m.mat
Loaded 1i0n.mat
Loaded 1i0o.mat
Loaded 1i0p.mat
Loaded 1i0q.mat
Loaded 1i2x.mat
Loaded 1i6h.mat
Loaded 1i6u.mat
Loaded 1i9v.mat
Loaded 1i9x.mat
Loaded 1icg.mat
Loaded 1ik5.mat
Loaded 1il2.mat
Loaded 1ivs.mat
Loaded 1j1u.mat
Loaded 1j2b.mat
Loaded 1j6s.mat
Loaded 1j7t.mat
Loaded 1j8g.mat
Loaded 1j9h.mat
Loaded 1jb8.mat
Loaded 1jbr.mat
Loaded 1jbs.mat
Loaded 1jbt.mat
Loaded 1jid.mat
Loaded 1jjm.mat
Loaded 1jjn.mat
Loaded 1jo2.mat
Loaded 1jzv.mat
Loaded 1k8w.mat
Loaded 1k9w.mat
Loaded 1kd3.mat
Loaded 1kd4.mat
Loaded 1kd5.mat
Loaded 1kfo.mat
Loaded 1kh6.mat
Loaded 1knz.mat
Loaded 1kog.mat
Loaded 1kq2.mat
Loaded 1kuo.mat
Loaded 1kxk.mat
Loaded 1l3z.mat
Loaded 1l8v.mat
Loaded 1l9a.mat
Loaded 1laj.mat
Loaded 1lc4.mat
Loaded 1lng.mat
Loaded 1lnt.mat
Loaded 1m5k.mat
Loaded 1m5o.mat
Loaded 1m5p.mat
Loaded 1m5v.mat
Loaded 1m8v.mat
Loaded 1m8x.mat
Loaded 1m8y.mat
Loaded 1mdg.mat
Loaded 1mfq.mat
Loaded 1mhk.mat
Loaded 1mji.mat
Loaded 1mms.mat
Loaded 1msw.mat
Loaded 1msy.mat
Loaded 1mzp.mat
Loaded 1n1h.mat
Loaded 1n35.mat
Loaded 1n38.mat
Loaded 1n77.mat
Loaded 1n78.mat
Loaded 1n7a.mat
Loaded 1n7b.mat
Loaded 1nbs.mat
Loaded 1nh3.mat
Loaded 1nb7.mat
Loaded 1nta.mat
Loaded 1ntb.mat
Loaded 1nyi.mat
Loaded 1o0b.mat
Loaded 1o0c.mat
Loaded 1o9m.mat
Loaded 1ob2.mat
Loaded 1ofx.mat
Loaded 1ooa.mat
Loaded 1p6v.mat
Loaded 1p79.mat
Loaded 1pgl.mat
Loaded 1pvo.mat
Loaded 1pwf.mat
Loaded 1q2r.mat
Loaded 1q2s.mat
Loaded 1q93.mat
Loaded 1qa6.mat
Loaded 1qbp.mat
Loaded 1qc0.mat
Loaded 1qf6.mat
Loaded 1qln.mat
Loaded 1qrs.mat
Loaded 1qrt.mat
Loaded 1qru.mat
Loaded 1qtq.mat
Loaded 1qu2.mat
Loaded 1qu3.mat
Loaded 1qzw.mat
Loaded 1r3e.mat
Loaded 1r3g.mat
Loaded 1r9f.mat
Loaded 1r9s.mat
Loaded 1r9t.mat
Loaded 1rc7.mat
Loaded 1rlg.mat
Loaded 1rmv.mat
Loaded 1rna.mat
Loaded 1rxa.mat
Loaded 1rxb.mat
Loaded 1s03.mat
Loaded 1s0v.mat
Loaded 1s76.mat
Loaded 1s77.mat
Loaded 1si2.mat
Loaded 1si3.mat
Loaded 1sz1.mat
Loaded 1t0d.mat
Loaded 1t0k.mat
Loaded 1tfw.mat
Loaded 1tfy.mat
Loaded 1ttt.mat
Loaded 1u0b.mat
Loaded 1u1y.mat
Loaded 1u63.mat
Loaded 1u6b.mat
Loaded 1u8d.mat
Loaded 1u9s.mat
Loaded 1un6.mat
Loaded 1urn.mat
Loaded 1utd.mat
Loaded 1utf.mat
Loaded 1utv.mat
Loaded 1uvi.mat
Loaded 1uvj.mat
Loaded 1uvm.mat
Loaded 1vfg.mat
Loaded 1wmq.mat
Loaded 1wsu.mat
Loaded 1x8w.mat
Loaded 1xjr.mat
Loaded 1xok.mat
Loaded 1xpo.mat
Loaded 1xpr.mat
Loaded 1xpu.mat
Loaded 1y0q.mat
Loaded 1y1w.mat
Loaded 1y26.mat
Loaded 1y27.mat
Loaded 1y39.mat
Loaded 1y77.mat
Loaded 1ykq.mat
Loaded 1ykv.mat
Loaded 1yls.mat
Loaded 1yrj.mat
Loaded 1ytu.mat
Loaded 1yvp.mat
Loaded 1zdh.mat
Loaded 1zdi.mat
Loaded 1zdj.mat
Loaded 1zdk.mat
Loaded 1ze2.mat
Loaded 1zjw.mat
Loaded 246d.mat
Loaded 259d.mat
Loaded 280d.mat
Loaded 283d.mat
Loaded 299d.mat
Loaded 2a8v.mat
Loaded 2bbv.mat
Loaded 2bgg.mat
Loaded 2bh2.mat
Loaded 2bj6.mat
Loaded 2fmt.mat
Loaded 300d.mat
Loaded 301d.mat
Loaded 315d.mat
Loaded 332d.mat
Loaded 333d.mat
Loaded 353d.mat
Loaded 354d.mat
Loaded 357d.mat
Loaded 359d.mat
Loaded 361d.mat
Loaded 364d.mat
Loaded 377d.mat
Loaded 379d.mat
Loaded 387d.mat
Loaded 397d.mat
Loaded 398d.mat
Loaded 402d.mat
Loaded 404d.mat
Loaded 405d.mat
Loaded 409d.mat
Loaded 418d.mat
Loaded 419d.mat
Loaded 420d.mat
Loaded 421d.mat
Loaded 422d.mat
Loaded 429d.mat
Loaded 430d.mat
Loaded 433d.mat
Loaded 434d.mat
Loaded 435d.mat
Loaded 438d.mat
Loaded 439d.mat
Loaded 464d.mat
Loaded 466d.mat
Loaded 472d.mat
Loaded 479d.mat
Loaded 480d.mat
Loaded 485d.mat
Loaded 488d.mat
Loaded 5msf.mat
Loaded 6msf.mat
Read 7msf.pdb for header information
Classifying interactions
Found   114 bases within 15 Angstroms from 3d structure
Found    48 pairs that are possibly interacting
Classification took 0.01 minutes, or 3234 classifications per minute
Saved 7msf.mat
Loaded 1f7y.mat

File = 

1x368 struct array with fields:
    CI
    ClassVersion
    Comment
    Distance
    Edge
    Filename
    HandClass
    Header
    Inter
    Modified
    NT
    NumNT
    Pair

>> PairViewer(File)

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [ 13] 1

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [1 ] 11

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
Found     12  AA  cSS cis Sugar edge / Sugar edge (second base dominant)  pairs.

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 1
Scheme 1 - Blue - cutoff matches; Green - hand matches; Orange - both agree
Scheme 2 - Cutoff classification -12 to 12
Scheme 3 - Cutoff classification 15 to 18
Scheme 4 - Cutoff classification -12 to 30
Scheme 5 - Pair code
Scheme 6 - Appropriate angle of rotation
Scheme 7 - Gap
Scheme 8 - Blue - cutoff matches; Green - exemplar matches; Orange - both agree
Scheme 9 - Distance to nearest exemplar
Scheme 10 - Exemplar classification
Scheme 11 - Degree of stacking overlap
Scheme 12 - Subcategory of first category
Scheme 13 - Distance to specified pair
Enter color scheme [9] 9
Enter 1 to show normal vector of second base [0] 
Enter 1 to show class limits [1] 

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 4
Pair   1 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.92  -1.75  -0.43  -0.82   34.7  -78.9  -0.84  11.00  158.9 
Hydrogen bonds:  158.9deg  2.46A |
q
Saved 7msf.mat

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [4] 9

 1 - File number
 2 - Paircode
 3 - Displacement(1)
 4 - Displacement(2)
 5 - Displacement(3)
 6 - Normal(3)
 7 - Angle between planes
 8 - Angle of rotation
 9 - Absolute value of Gap
10 - Computer classification
11 - Hand classification
12 - Minimum distance
13 - Index of first nucleotide
14 - Lowest index in pair
15 - Angle in the plane of the first base
16 - Smallest hydrogen bond angle
17 - Largest hydrogen bond angle
18 - Smallest hydrogen bond length
19 - Largest hydrogen bond length
20 - Distance to nearest exemplar
21 - Classification of nearest exemplar
22 - Degree of stacking overlap
23 - Discrepancy with specified pair

Enter sort key(s) in brackets, negative for reversed order
20

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [4] 2
       Filename  Nucl1  Nucl2  Disp1  Disp2  Disp3  Norm3   Ang    Gap MinDist Cutoff   Hand   Exem  EDistPairDisc
  1        1et4 A  434 A  535   6.57  -1.69   1.61  -0.93 -82.4   0.75    2.51  11.00   0.00  11.00   0.91   0.00
  2        1et4 A  234 A  335   6.58  -2.05   1.75  -0.92 -82.8   0.78    2.66  11.00   0.00  11.00   1.23   0.00
  3        1j5e A   16 A 1080   6.79  -1.27   1.17  -0.99 -71.1   1.39    2.57  11.00   0.00  11.00   1.60   0.00
  4        1s72 A 1909 A 2266   6.68  -1.95   0.47  -0.75 -83.5  -0.86    2.47  11.00   0.00  11.00   1.78   0.00
  5        2AVY A   16 A 1080   6.51  -1.60   0.62  -0.95 -78.0   1.25    2.55  11.00   0.00  11.00   2.08   0.00
  6        2AVY A  909 A 1413   7.28  -1.92   1.29  -0.67 -83.7  -1.84    2.89  11.00   0.00  11.00   2.17   0.00
  7        2AVY A   66 A  172   6.96  -2.68   0.92  -0.73 -80.7  -1.96    2.96  11.00   0.00  11.00   2.30   0.00
  8        1j5e A  130 A  263   6.51  -2.26   0.63  -0.94 -76.5   1.40    2.80  11.00   0.00  11.00   2.51   0.00
  9        1j5e A  909 A 1413   6.71  -1.78   1.95  -0.65 269.4   0.10    2.89  11.00   0.00  11.00   2.61   0.00
 10        1s72 A 1132 A 2521   6.92  -1.75  -0.43  -0.82 -78.9  -0.84    2.46  11.00   0.00  11.00   2.73   0.00
 11        2AVY A  130 A  263   6.64  -2.47   0.55  -0.89 -77.6   1.53    3.12  11.00   0.00  11.00   2.99   0.00
 12        2AW4 A 2077 A 2434   7.07  -2.65  -0.23  -0.91 -74.8   1.58    3.08  11.00   0.00  11.00   3.86   0.00

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [2] 4
Pair   1 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.57  -1.69   1.61  -0.93   21.0  -82.4   0.75  11.00  161.4 
Hydrogen bonds:  161.4deg  2.51A |
s
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.57  -1.69   1.61  -0.93   21.0  -82.4   0.75  11.00  161.4 

Pair   2 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.58  -2.05   1.75  -0.92   23.6  -82.8   0.78  11.00  157.1 
Hydrogen bonds:  157.1deg  2.66A |
rr
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.58  -2.05   1.75  -0.92   23.6  -82.8   0.78  11.00  157.1 

Pair   3 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.79  -1.27   1.17  -0.99    6.6  -71.1   1.39  11.00  146.7 
Hydrogen bonds:  146.7deg  2.57A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.79  -1.27   1.17  -0.99    6.6  -71.1   1.39  11.00  146.7 

Pair   4 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.68  -1.95   0.47  -0.75   41.4  -83.5  -0.86  11.00  149.8 
Hydrogen bonds:  149.8deg  2.47A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.68  -1.95   0.47  -0.75   41.4  -83.5  -0.86  11.00  149.8 

Pair   5 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.51  -1.60   0.62  -0.95   18.0  -78.0   1.25  11.00  150.4 
Hydrogen bonds:  150.4deg  2.55A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.51  -1.60   0.62  -0.95   18.0  -78.0   1.25  11.00  150.4 

Pair   6 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.28  -1.92   1.29  -0.67   47.6  -83.7  -1.84  11.00  157.7 
Hydrogen bonds:  157.7deg  2.89A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.28  -1.92   1.29  -0.67   47.6  -83.7  -1.84  11.00  157.7 

Pair   7 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.96  -2.68   0.92  -0.73   43.3  -80.7  -1.96  11.00  143.3 
Hydrogen bonds:  143.3deg  3.06A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.96  -2.68   0.92  -0.73   43.3  -80.7  -1.96  11.00  143.3 

Pair   8 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.51  -2.26   0.63  -0.94   19.6  -76.5   1.40  11.00  144.0 
Hydrogen bonds:  144.0deg  2.80A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.51  -2.26   0.63  -0.94   19.6  -76.5   1.40  11.00  144.0 

Pair   9 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.71  -1.78   1.95  -0.65   49.4  269.4   0.10  11.00  159.2 
Hydrogen bonds:  159.2deg  2.89A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.71  -1.78   1.95  -0.65   49.4  269.4   0.10  11.00  159.2 

Pair  10 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.92  -1.75  -0.43  -0.82   34.7  -78.9  -0.84  11.00  158.9 
Hydrogen bonds:  158.9deg  2.46A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.92  -1.75  -0.43  -0.82   34.7  -78.9  -0.84  11.00  158.9 

Pair  11 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.64  -2.47   0.55  -0.89   26.7  -77.6   1.53  11.00  142.7 
Hydrogen bonds:  142.7deg  3.12A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.64  -2.47   0.55  -0.89   26.7  -77.6   1.53  11.00  142.7 

Pair  12 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.07  -2.65  -0.23  -0.91   24.1  -74.8   1.58  11.00  151.4 
Hydrogen bonds:  151.4deg  3.08A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.07  -2.65  -0.23  -0.91   24.1  -74.8   1.58  11.00  151.4 

Wrote 1s72.hand
Saved 1s72.mat
Saved 1j5e.mat
Saved 2AW4.mat
Saved 2AVY.mat
Saved 1et4.mat

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [4] 1
Scheme 1 - Blue - cutoff matches; Green - hand matches; Orange - both agree
Scheme 2 - Cutoff classification -12 to 12
Scheme 3 - Cutoff classification 15 to 18
Scheme 4 - Cutoff classification -12 to 30
Scheme 5 - Pair code
Scheme 6 - Appropriate angle of rotation
Scheme 7 - Gap
Scheme 8 - Blue - cutoff matches; Green - exemplar matches; Orange - both agree
Scheme 9 - Distance to nearest exemplar
Scheme 10 - Exemplar classification
Scheme 11 - Degree of stacking overlap
Scheme 12 - Subcategory of first category
Scheme 13 - Distance to specified pair
Enter color scheme [9] 6
Enter 1 to show normal vector of second base [0] 
Enter 1 to show class limits [1] 

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 2
       Filename  Nucl1  Nucl2  Disp1  Disp2  Disp3  Norm3   Ang    Gap MinDist Cutoff   Hand   Exem  EDistPairDisc
  1        1et4 A  535 A  434   6.57  -1.69   1.61  -0.93 -82.4   0.75    2.51  11.00   0.00  11.00   0.00   0.00
  2        1et4 A  335 A  234   6.58  -2.05   1.75  -0.92 -82.8   0.78    2.66  11.00   0.00  11.00   0.53   0.00
  3        1j5e A 1080 A   16   6.79  -1.27   1.17  -0.99 -71.1   1.39    2.57  11.00   0.00  11.00   2.01   0.00
  4        1s72 A 2266 A 1909   6.68  -1.95   0.47  -0.75 -83.5  -0.86    2.47  11.00   0.00  11.00   2.36   0.00
  5        2AVY A 1080 A   16   6.51  -1.60   0.62  -0.95 -78.0   1.25    2.55  11.00   0.00  11.00   2.79   0.00
  6        2AVY A 1413 A  909   7.28  -1.92   1.29  -0.67 -83.7  -1.84    2.89  11.00   0.00  11.00   2.24   0.00
  7        2AVY A  172 A   66   6.96  -2.68   0.92  -0.73 -80.7  -1.96    2.96  11.00   0.00  11.00   2.51   0.00
  8        1j5e A  263 A  130   6.51  -2.26   0.63  -0.94 -76.5   1.40    2.80  11.00   0.00  11.00   3.09   0.00
  9        1j5e A 1413 A  909   6.71  -1.78   1.95  -0.65 269.4   0.10    2.89  11.00   0.00  11.00   1.87   0.00
 10        1s72 A 2521 A 1132   6.92  -1.75  -0.43  -0.82 -78.9  -0.84    2.46  11.00   0.00  11.00   3.59   0.00
 11        2AVY A  263 A  130   6.64  -2.47   0.55  -0.89 -77.6   1.53    3.12  11.00   0.00  11.00   3.58   0.00
 12        2AW4 A 2434 A 2077   7.07  -2.65  -0.23  -0.91 -74.8   1.58    3.08  11.00   0.00  11.00   4.42   0.00

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [2] 4
Pair   1 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.57  -1.69   1.61  -0.93   21.0  -82.4   0.75  11.00  161.4 
Hydrogen bonds:  161.4deg  2.51A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.31  -5.98   1.36  -0.93   21.0  -82.4   0.60 111.00   88.9 

Pair   2 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.58  -2.05   1.75  -0.92   23.6  -82.8   0.78  11.00  157.1 
Hydrogen bonds:  157.1deg  2.66A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.75  -5.92   1.19  -0.92   23.6  -82.8   0.41 111.00   93.7 

Pair   3 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.79  -1.27   1.17  -0.99    6.6  -71.1   1.39  11.00  146.7 
Hydrogen bonds:  146.7deg  2.57A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.27  -5.92   1.82  -0.99    6.6  -71.1   1.57 111.00   84.9 

Pair   4 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.68  -1.95   0.47  -0.75   41.4  -83.5  -0.86  11.00  149.8 
Hydrogen bonds:  149.8deg  2.47A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  4.19  -5.40   1.40  -0.75   41.4  -83.5  -0.75 111.00   93.1 

Pair   5 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.51  -1.60   0.62  -0.95   18.0  -78.0   1.25  11.00  150.4 
Hydrogen bonds:  150.4deg  2.55A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  2.69  -5.76   2.20  -0.95   18.0  -78.0   1.73 111.00   81.3 

Pair   6 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.28  -1.92   1.29  -0.67   47.6  -83.7  -1.84  11.00  157.7 
Hydrogen bonds:  157.7deg  2.89A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  5.02  -5.66   1.04  -0.67   47.6  -83.7  -0.79 111.00  101.5 

Pair   7 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.96  -2.68   0.92  -0.73   43.3  -80.7  -1.96  11.00  143.3 
Hydrogen bonds:  143.3deg  3.06A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  5.29  -5.32   0.29  -0.73   43.3  -80.7  -1.42 111.00  105.1 

Pair   8 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.51  -2.26   0.63  -0.94   19.6  -76.5   1.40  11.00  144.0 
Hydrogen bonds:  144.0deg  2.80A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.46  -5.58   2.21  -0.94   19.6  -76.5   2.00 111.00   88.0 

Pair   9 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.71  -1.78   1.95  -0.65   49.4  269.4   0.10  11.00  159.2 
Hydrogen bonds:  159.2deg  2.89A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.73  -6.15  -0.51  -0.65   49.4  269.4  -1.23  30.00   92.9 

Pair  10 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.92  -1.75  -0.43  -0.82   34.7  -78.9  -0.84  11.00  158.9 
Hydrogen bonds:  158.9deg  2.46A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.89  -5.31   2.81  -0.82   34.7  -78.9   0.28 111.00   91.1 

Pair  11 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.64  -2.47   0.55  -0.89   26.7  -77.6   1.53  11.00  142.7 
Hydrogen bonds:  142.7deg  3.12A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.43  -5.44   3.02  -0.89   26.7  -77.6   2.49 111.00   87.9 

Pair  12 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.07  -2.65  -0.23  -0.91   24.1  -74.8   1.58  11.00  151.4 
Hydrogen bonds:  151.4deg  3.08A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  4.69  -5.88   0.59  -0.91   24.1  -74.8   1.29 111.00  100.5 

Wrote 1s72.hand
Saved 1s72.mat
Saved 1j5e.mat
Saved 2AW4.mat
Saved 2AVY.mat
Saved 1et4.mat

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [4] 2
       Filename  Nucl1  Nucl2  Disp1  Disp2  Disp3  Norm3   Ang    Gap MinDist Cutoff   Hand   Exem  EDistPairDisc
  1        1et4 A  434 A  535   3.31  -5.98   1.36  -0.93 -82.4   0.60    2.51 111.00   0.00  11.00   6.40   0.00
  2        1et4 A  234 A  335   3.75  -5.92   1.19  -0.92 -82.8   0.41    2.66 111.00   0.00  11.00   6.21   0.00
  3        1j5e A   16 A 1080   3.27  -5.92   1.82  -0.99 -71.1   1.57    2.57 111.00   0.00  11.00   6.45   0.00
  4        1s72 A 1909 A 2266   4.19  -5.40   1.40  -0.75 -83.5  -0.75    2.47 111.00   0.00  11.00   5.80   0.00
  5        2AVY A   16 A 1080   2.69  -5.76   2.20  -0.95 -78.0   1.73    2.55 111.00   0.00  11.00   6.87   0.00
  6        2AVY A  909 A 1413   5.02  -5.66   1.04  -0.67 -83.7  -0.79    2.89 111.00   0.00  11.00   6.19   0.00
  7        2AVY A   66 A  172   5.29  -5.32   0.29  -0.73 -80.7  -1.42    2.96 111.00   0.00  11.00   5.77   0.00
  8        1j5e A  130 A  263   3.46  -5.58   2.21  -0.94 -76.5   2.00    2.80 111.00   0.00  11.00   6.46   0.00
  9        1j5e A  909 A 1413   3.73  -6.15  -0.51  -0.65 269.4  -1.23    2.89  30.00   0.00  11.00   8.07   0.00
 10        1s72 A 1132 A 2521   3.89  -5.31   2.81  -0.82 -78.9   0.28    2.46 111.00   0.00  11.00   5.45   0.00
 11        2AVY A  130 A  263   3.43  -5.44   3.02  -0.89 -77.6   2.49    3.12 111.00   0.00  11.00   6.69   0.00
 12        2AW4 A 2077 A 2434   4.69  -5.88   0.59  -0.91 -74.8   1.29    3.08 111.00   0.00  11.00   6.71   0.00

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [2] 1
Scheme 1 - Blue - cutoff matches; Green - hand matches; Orange - both agree
Scheme 2 - Cutoff classification -12 to 12
Scheme 3 - Cutoff classification 15 to 18
Scheme 4 - Cutoff classification -12 to 30
Scheme 5 - Pair code
Scheme 6 - Appropriate angle of rotation
Scheme 7 - Gap
Scheme 8 - Blue - cutoff matches; Green - exemplar matches; Orange - both agree
Scheme 9 - Distance to nearest exemplar
Scheme 10 - Exemplar classification
Scheme 11 - Degree of stacking overlap
Scheme 12 - Subcategory of first category
Scheme 13 - Distance to specified pair
Enter color scheme [6] 
Enter 1 to show normal vector of second base [0] 
Enter 1 to show class limits [1] 

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] Error in ==> zEnterViewMode at 46
inp = input('');

Error in ==> PairViewer at 37
        ViewParam = zEnterViewMode(Param,ViewParam);


>> zGetNTData('7msf   ')
>> PairViewer(File)

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [ 13] 1

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [1 ] 11

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 1
Found     12  AA  cSS cis Sugar edge / Sugar edge (second base dominant)  pairs.

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 1
Scheme 1 - Blue - cutoff matches; Green - hand matches; Orange - both agree
Scheme 2 - Cutoff classification -12 to 12
Scheme 3 - Cutoff classification 15 to 18
Scheme 4 - Cutoff classification -12 to 30
Scheme 5 - Pair code
Scheme 6 - Appropriate angle of rotation
Scheme 7 - Gap
Scheme 8 - Blue - cutoff matches; Green - exemplar matches; Orange - both agree
Scheme 9 - Distance to nearest exemplar
Scheme 10 - Exemplar classification
Scheme 11 - Degree of stacking overlap
Scheme 12 - Subcategory of first category
Scheme 13 - Distance to specified pair
Enter color scheme [9] 
Enter 1 to show normal vector of second base [0] 
Enter 1 to show class limits [1] 

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 0

ans = 

1x368 struct array with fields:
    CI
    ClassVersion
    Comment
    Distance
    Edge
    Filename
    HandClass
    Header
    Inter
    Modified
    NT
    NumNT
    Pair

>> File = PairViewer(File)

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [ 13] 1

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [1 ] 11

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
Found     12  AA  cSS cis Sugar edge / Sugar edge (second base dominant)  pairs.

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 4
Pair   1 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.92  -1.75  -0.43  -0.82   34.7  -78.9  -0.84  11.00  158.9 
Hydrogen bonds:  158.9deg  2.46A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.92  -1.75  -0.43  -0.82   34.7  -78.9  -0.84 111.00  158.9 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.89  -5.31   2.81  -0.82   34.7  -78.9   0.28 111.00   91.1 

Pair   2 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.68  -1.95   0.47  -0.75   41.4  -83.5  -0.86  11.00  149.8 
Hydrogen bonds:  149.8deg  2.47A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.68  -1.95   0.47  -0.75   41.4  -83.5  -0.86 111.00  149.8 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  4.19  -5.40   1.40  -0.75   41.4  -83.5  -0.75 111.00   93.1 

Pair   3 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.79  -1.27   1.17  -0.99    6.6  -71.1   1.39  11.00  146.7 
Hydrogen bonds:  146.7deg  2.57A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.79  -1.27   1.17  -0.99    6.6  -71.1   1.39 111.00  146.7 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.27  -5.92   1.82  -0.99    6.6  -71.1   1.57 111.00   84.9 

Pair   4 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.51  -2.26   0.63  -0.94   19.6  -76.5   1.40  11.00  144.0 
Hydrogen bonds:  144.0deg  2.80A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.51  -2.26   0.63  -0.94   19.6  -76.5   1.40 111.00  144.0 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.46  -5.58   2.21  -0.94   19.6  -76.5   2.00 111.00   88.0 

Pair   5 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.71  -1.78   1.95  -0.65   49.4  269.4   0.10  11.00  159.2 
Hydrogen bonds:  159.2deg  2.89A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.71  -1.78   1.95  -0.65   49.4  269.4   0.10 111.00  159.2 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.73  -6.15  -0.51  -0.65   49.4  269.4  -1.23  30.40   92.9 

Pair   6 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.07  -2.65  -0.23  -0.91   24.1  -74.8   1.58  11.00  151.4 
Hydrogen bonds:  151.4deg  3.08A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.07  -2.65  -0.23  -0.91   24.1  -74.8   1.58 111.00  151.4 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  4.69  -5.88   0.59  -0.91   24.1  -74.8   1.29 111.00  100.5 

Pair   7 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.51  -1.60   0.62  -0.95   18.0  -78.0   1.25  11.00  150.4 
Hydrogen bonds:  150.4deg  2.55A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.51  -1.60   0.62  -0.95   18.0  -78.0   1.25 111.00  150.4 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  2.69  -5.76   2.20  -0.95   18.0  -78.0   1.73 111.00   81.3 

Pair   8 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.96  -2.68   0.92  -0.73   43.3  -80.7  -1.96  11.00  143.3 
Hydrogen bonds:  143.3deg  3.06A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.96  -2.68   0.92  -0.73   43.3  -80.7  -1.96 111.00  143.3 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  5.29  -5.32   0.29  -0.73   43.3  -80.7  -1.42 111.00  105.1 

Pair   9 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.64  -2.47   0.55  -0.89   26.7  -77.6   1.53  11.00  142.7 
Hydrogen bonds:  142.7deg  3.12A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.64  -2.47   0.55  -0.89   26.7  -77.6   1.53 111.00  142.7 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.43  -5.44   3.02  -0.89   26.7  -77.6   2.49 111.00   87.9 

Pair  10 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.28  -1.92   1.29  -0.67   47.6  -83.7  -1.84  11.00  157.7 
Hydrogen bonds:  157.7deg  2.89A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.28  -1.92   1.29  -0.67   47.6  -83.7  -1.84 111.00  157.7 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  5.02  -5.66   1.04  -0.67   47.6  -83.7  -0.79 111.00  101.5 

Pair  11 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.58  -2.05   1.75  -0.92   23.6  -82.8   0.78  11.00  157.1 
Hydrogen bonds:  157.1deg  2.66A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.58  -2.05   1.75  -0.92   23.6  -82.8   0.78 111.00  157.1 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.75  -5.92   1.19  -0.92   23.6  -82.8   0.41 111.00   93.7 

Pair  12 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.57  -1.69   1.61  -0.93   21.0  -82.4   0.75  11.00  161.4 
Hydrogen bonds:  161.4deg  2.51A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.57  -1.69   1.61  -0.93   21.0  -82.4   0.75 111.00  161.4 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  3.31  -5.98   1.36  -0.93   21.0  -82.4   0.60 111.00   88.9 

Wrote 1s72.hand
Saved 1s72.mat
Saved 1j5e.mat
Saved 2AW4.mat
Saved 2AVY.mat
Saved 1et4.mat
Saved 7msf.mat

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 1
Scheme 1 - Blue - cutoff matches; Green - hand matches; Orange - both agree
Scheme 2 - Cutoff classification -12 to 12
Scheme 3 - Cutoff classification 15 to 18
Scheme 4 - Cutoff classification -12 to 30
Scheme 5 - Pair code
Scheme 6 - Appropriate angle of rotation
Scheme 7 - Gap
Scheme 8 - Blue - cutoff matches; Green - exemplar matches; Orange - both agree
Scheme 9 - Distance to nearest exemplar
Scheme 10 - Exemplar classification
Scheme 11 - Degree of stacking overlap
Scheme 12 - Subcategory of first category
Scheme 13 - Distance to specified pair
Enter color scheme [9] 
Enter 1 to show normal vector of second base [0] 
Enter 1 to show class limits [1] 

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 8

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [  1] 11

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [11 ] 11

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
Found      1  GG  cSS cis Sugar edge / Sugar edge (second base dominant)  pairs.

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 2
       Filename  Nucl1  Nucl2  Disp1  Disp2  Disp3  Norm3   Ang    Gap MinDist Cutoff   Hand   Exem  EDistPairDisc
  1        1s72 G  885 G 2113   5.14  -3.40   1.69  -0.86 252.9   2.09    2.36  11.00   0.00  11.00   0.03   0.00

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [2] 4
Pair   1 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  5.14  -3.40   1.69  -0.86   31.0  252.9   2.09  11.00  124.6 
Hydrogen bonds:  124.6deg  2.62A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  5.14  -3.40   1.69  -0.86   31.0  252.9   2.09  11.00  124.6 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  2.06  -5.89  -1.36  -0.86   31.0  252.9   1.13 111.00  138.4 

Wrote 1s72.hand
Saved 1s72.mat

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [4] 8

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [ 11] 1

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [11 ] 12

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
Found     11  AA  tSS trans Sugar edge / Sugar edge (second base dominant)  pairs.

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 2
       Filename  Nucl1  Nucl2  Disp1  Disp2  Disp3  Norm3   Ang    Gap MinDist Cutoff   Hand   Exem  EDistPairDisc
  1        1s72 A  340 A  306   7.52  -1.85   1.06   0.84 156.0  -0.73    2.44  12.00   0.00  12.00   3.36   0.00
  2        1s72 A 2311 A 1007   7.06  -2.45  -1.18   0.82 144.8   1.38    2.49  12.00   0.00  12.00   4.90   0.00
  3        1s72 A 2604 A 1233   7.40  -2.29   0.76   0.74 153.7  -1.04    2.85  12.00   0.00  12.00   3.13   0.00
  4        1s72 A 1242 A 1247   7.53  -3.79   0.71   0.96 132.9  -0.31    2.49  12.00   0.00  12.00   0.06   0.00
  5        1s72 A 1313 A 1341   7.58  -3.30  -0.45   0.76 139.2  -1.05    2.49  12.00   0.00  12.00   2.84   0.00
  6        1s72 A 1829 A 2019   7.53  -2.20   1.86   0.97 138.0   2.55    2.90  12.00   0.00  12.00   3.46   0.00
  7        1j5e A  919 A 1080   7.48  -3.98   1.26   0.95 128.7   0.10    2.52  12.00   0.00  12.00   0.83   0.00
  8        2AW4 A 1260 A  515   7.28  -2.52  -0.38   0.96 149.4   0.87    2.20  12.00   0.00  12.00   3.38   0.00
  9        2AW4 A 1773 A 1978   7.37  -2.36   0.91   0.92 137.4   2.18    2.64  12.00   0.00  12.00   3.37   0.00
 10        2AVY A   19 A  572   7.61  -1.92   1.17   0.75 152.0  -0.68    2.40  12.00   0.00  12.00   3.45   0.00
 11        1mfq A  127 A  173   7.99  -2.72  -0.20   0.79 148.7  -0.63    2.51  12.00   0.00  12.00   0.00   0.00

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [2] 4
Pair   1 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.52  -1.85   1.06   0.84   33.2  156.0  -0.73  12.00  122.3 
Hydrogen bonds:  122.3deg  2.44A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.52  -1.85   1.06   0.84   33.2  156.0  -0.73  12.00  122.3 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.21   1.69   2.50   0.84   33.2  204.0   1.58  30.00  121.6 

Pair   2 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.06  -2.45  -1.18   0.82   35.0  144.8   1.38  12.00  114.5 
Hydrogen bonds:  114.5deg  2.49A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.06  -2.45  -1.18   0.82   35.0  144.8   1.38  12.00  114.5 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.46   2.14  -3.31   0.82   35.0  215.2  -2.15  30.00 

Pair   3 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.40  -2.29   0.76   0.74   42.3  153.7  -1.04  12.00  110.9 
Hydrogen bonds:  110.9deg  2.85A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.40  -2.29   0.76   0.74   42.3  153.7  -1.04  12.00  110.9 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.91   1.99   2.97   0.74   42.3  206.3   2.54  30.00  117.1 

Pair   4 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.53  -3.79   0.71   0.96   16.0  132.9  -0.31  12.00  157.0 
Hydrogen bonds:  157.0deg  2.49A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.53  -3.79   0.71   0.96   16.0  132.9  -0.31  12.00  157.0 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.81   3.02   1.15   0.96   16.0  227.1   0.97  30.00   89.9   56.0 

Pair   5 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.58  -3.30  -0.45   0.76   40.4  139.2  -1.05  12.00  128.8 
Hydrogen bonds:  128.8deg  2.49A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.58  -3.30  -0.45   0.76   40.4  139.2  -1.05  12.00  128.8 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.47   3.36   1.19   0.76   40.4  220.8   2.71  30.00   90.1   57.2 

Pair   6 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.53  -2.20   1.86   0.97   14.1  138.0   2.55  12.00  104.2 
Hydrogen bonds:  104.2deg  2.90A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.53  -2.20   1.86   0.97   14.1  138.0   2.55  12.00  104.2 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.57   3.76  -2.78   0.97   14.1  222.0  -2.83 102.00   85.2   58.7 

Pair   7 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.48  -3.98   1.26   0.95   18.7  128.7   0.10  12.00  162.0 
Hydrogen bonds:  162.0deg  2.52A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.48  -3.98   1.26   0.95   18.7  128.7   0.10  12.00  162.0 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.81   3.42   0.81   0.95   18.7  231.3   0.68  30.00   93.3   56.1 

Pair   8 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.28  -2.52  -0.38   0.96   16.7  149.4   0.87  12.00  135.7 
Hydrogen bonds:  135.7deg  2.20A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.28  -2.52  -0.38   0.96   16.7  149.4   0.87  12.00  135.7 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.37   1.48  -1.73   0.96   16.7  210.6  -0.50  30.00  115.6 

Pair   9 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.37  -2.36   0.91   0.92   22.4  137.4   2.18  12.00  102.4 
Hydrogen bonds:  102.4deg  2.64A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.37  -2.36   0.91   0.92   22.4  137.4   2.18  12.00  102.4 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  6.32   3.65  -2.73   0.92   22.4  222.6  -2.58  30.00   83.5   57.0 

Pair  10 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.61  -1.92   1.17   0.75   41.0  152.0  -0.68  12.00  121.2 
Hydrogen bonds:  121.2deg  2.40A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.61  -1.92   1.17   0.75   41.0  152.0  -0.68  12.00  121.2 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.16   2.40   2.46   0.75   41.0  208.0   1.99  30.00  124.1 

Pair  11 | ret: next | q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify | r: reverse
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.99  -2.72  -0.20   0.79   38.3  148.7  -0.63  12.00  138.2 
Hydrogen bonds:  138.2deg  2.51A |
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.99  -2.72  -0.20   0.79   38.3  148.7  -0.63  12.00  138.2 
r
 Sh(1)  Sh(2)  Sh(3) Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles
  7.96   2.53   1.22   0.79   38.3  211.3   2.36  30.00  124.0 

Wrote 1s72.hand
Saved 1s72.mat
Saved 1j5e.mat
Saved 2AW4.mat
Saved 2AVY.mat
Saved 1mfq.mat

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [4] 2
       Filename  Nucl1  Nucl2  Disp1  Disp2  Disp3  Norm3   Ang    Gap MinDist Cutoff   Hand   Exem  EDistPairDisc
  1        1s72 A  340 A  306   7.21   1.69   2.50   0.84 204.0   1.58    2.44  30.00   0.00  12.00   8.29   0.00
  2        1s72 A 2311 A 1007   6.46   2.14  -3.31   0.82 215.2  -2.15    2.49  30.00   0.00  23.10   8.76   0.00
  3        1s72 A 2604 A 1233   6.91   1.99   2.97   0.74 206.3   2.54    2.85  30.00   0.00  12.00   8.84   0.00
  4        1s72 A 1242 A 1247   7.81   3.02   1.15   0.96 227.1   0.97    2.49  30.00   0.00   2.00   9.33   0.00
  5        1s72 A 1313 A 1341   7.47   3.36   1.19   0.76 220.8   2.71    2.49  30.00   0.00   2.00   8.98   0.00
  6        1s72 A 1829 A 2019   6.57   3.76  -2.78   0.97 222.0  -2.83    2.90 102.00   0.00   2.00   7.87   0.00
  7        1j5e A  919 A 1080   7.81   3.42   0.81   0.95 231.3   0.68    2.52  30.00   0.00   2.00   9.06   0.00
  8        2AW4 A 1260 A  515   7.37   1.48  -1.73   0.96 210.6  -0.50    2.20  30.00   0.00  12.00   8.22   0.00
  9        2AW4 A 1773 A 1978   6.32   3.65  -2.73   0.92 222.6  -2.58    2.64  30.00   0.00   2.00   8.14   0.00
 10        2AVY A   19 A  572   7.16   2.40   2.46   0.75 208.0   1.99    2.40  30.00   0.00  12.00   8.98   0.00
 11        1mfq A  127 A  173   7.96   2.53   1.22   0.79 211.3   2.36    2.51  30.00   0.00  12.00   8.56   0.00

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [2] 1
Scheme 1 - Blue - cutoff matches; Green - hand matches; Orange - both agree
Scheme 2 - Cutoff classification -12 to 12
Scheme 3 - Cutoff classification 15 to 18
Scheme 4 - Cutoff classification -12 to 30
Scheme 5 - Pair code
Scheme 6 - Appropriate angle of rotation
Scheme 7 - Gap
Scheme 8 - Blue - cutoff matches; Green - exemplar matches; Orange - both agree
Scheme 9 - Distance to nearest exemplar
Scheme 10 - Exemplar classification
Scheme 11 - Degree of stacking overlap
Scheme 12 - Subcategory of first category
Scheme 13 - Distance to specified pair
Enter color scheme [9] 
Enter 1 to show normal vector of second base [0] 
Enter 1 to show class limits [1] 

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 0

File = 

1x368 struct array with fields:
    CI
    ClassVersion
    Comment
    Distance
    Edge
    Filename
    HandClass
    Header
    Inter
    Modified
    NT
    NumNT
    Pair

>> File = PairViewer(File)
>> File = PairViewer(File)

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [ 13] 1

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [1 ] 12

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
No  AA  tSS trans Sugar edge / Sugar edge (second base dominant)  pairs were found.

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [  1] 

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [12 ] -12

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 

No  AA  tSS trans Sugar edge / Sugar edge (first base dominant)  pairs were found.

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [  1] 
  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [-12 ] Error in ==> zEnterPairSelection at 55
inp = lower(input('','s'));

Error in ==> PairViewer at 23
      Param = zEnterPairSelection(Param);

>> 11
>> File = PairViewer(File)

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [ 13] 11

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [1 ] 11

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
No  GG  cSS cis Sugar edge / Sugar edge (second base dominant)  pairs were found.

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [ 11] 1

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [11 ] 11

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
No  AA  cSS cis Sugar edge / Sugar edge (second base dominant)  pairs were found.

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [  1] Error in ==> zEnterPairSelection at 32
inp = input('');

Error in ==> PairViewer at 23
      Param = zEnterPairSelection(Param);

>> 
>> File = zGetNTData({'1s72','1j53','2aw4','2avy','1mfq','1et4','7msf'},1)
Loaded 1s72.mat
Classifying interactions
Found 29027 bases within 15 Angstroms from 3d structure
Found 12235 pairs that are possibly interacting
Classification took 1.09 minutes, or 11242 classifications per minute
Saved 1s72.mat
Could not open file 1j53.pdb
Loaded 2aw4.mat
Classifying interactions
Found 29431 bases within 15 Angstroms from 3d structure
Found 12185 pairs that are possibly interacting
Classification took 1.08 minutes, or 11269 classifications per minute
Saved 2AW4.mat
Loaded 2avy.mat
Classifying interactions
Found 15141 bases within 15 Angstroms from 3d structure
Bases C  948(A) and C 1234(A) fall into categories   1.00  23.10 
Found  6348 pairs that are possibly interacting
Classification took 0.52 minutes, or 12231 classifications per minute
Saved 2AVY.mat
Loaded 1mfq.mat
Classifying interactions
Found  1016 bases within 15 Angstroms from 3d structure
Found   393 pairs that are possibly interacting
Classification took 0.04 minutes, or 10061 classifications per minute
Saved 1mfq.mat
Loaded 1et4.mat
Classifying interactions
Found  1419 bases within 15 Angstroms from 3d structure
Found   624 pairs that are possibly interacting
Classification took 0.06 minutes, or 11303 classifications per minute
Saved 1et4.mat
Loaded 7msf.mat
Classifying interactions
Found   114 bases within 15 Angstroms from 3d structure
Found    48 pairs that are possibly interacting
Classification took 0.01 minutes, or 3291 classifications per minute
Saved 7msf.mat

File = 

1x7 struct array with fields:
    CI
    ClassVersion
    Comment
    Distance
    Edge
    Filename
    HandClass
    Header
    Inter
    Modified
    NT
    NumNT
    Pair

>> PairViewer(File)

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [ 13] 1

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [1 ] 11

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
No  AA  cSS cis Sugar edge / Sugar edge (second base dominant)  pairs were found.

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [  1] 11

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [11 ] 11

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
No  GG  cSS cis Sugar edge / Sugar edge (second base dominant)  pairs were found.

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [ 11] 1

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [11 ] 12

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
Found      3  AA  tSS trans Sugar edge / Sugar edge (second base dominant)  pairs.

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 2
       Filename  Nucl1  Nucl2  Disp1  Disp2  Disp3  Norm3   Ang    Gap MinDist Cutoff   Hand   Exem  EDistPairDisc
  1        1s72 A  306 A  340   7.21   1.69   2.50   0.84 204.0   1.58    2.44  12.00   0.00  12.00   8.29   0.00
  2        1s72 A 1233 A 2604   6.91   1.99   2.97   0.74 206.3   2.54    2.85  12.00   0.00  12.00   8.84   0.00
  3        2AW4 A  515 A 1260   7.37   1.48  -1.73   0.96 210.6  -0.50    2.20  12.00   0.00  12.00   8.22   0.00

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 0

ans = 

1x7 struct array with fields:
    CI
    ClassVersion
    Comment
    Distance
    Edge
    Filename
    HandClass
    Header
    Inter
    Modified
    NT
    NumNT
    Pair

>> File

File = 

1x7 struct array with fields:
    CI
    ClassVersion
    Comment
    Distance
    Edge
    Filename
    HandClass
    Header
    Inter
    Modified
    NT
    NumNT
    Pair

>> File.Filename

ans =

1s72


ans =

     []


ans =

2AW4


ans =

2AVY


ans =

1mfq


ans =

1et4


ans =

7msf

>> 
>> File = zGetNTData({'1s72','1j5e','2aw4','2avy','1mfq','1et4','7msf'},1)
Error in ==> zGetNTData at 46
      load(strcat(Filename,'.mat'),'File','-mat');


>> zGetNTData({'1j5e'},1);
Loaded 1j5e.mat
Classifying interactions
Found 15003 bases within 15 Angstroms from 3d structure
Bases C 1437(A) and C 1465(A) fall into categories   1.00  23.10 
Found  6352 pairs that are possibly interacting
Classification took 0.52 minutes, or 12129 classifications per minute
Saved 1j5e.mat
>> File = zGetNTData({'1s72','1j5e','2aw4','2avy','1mfq','1et4','7msf'},0)
Loaded 1s72.mat
Loaded 1j5e.mat
Loaded 2aw4.mat
Loaded 2avy.mat
Loaded 1mfq.mat
Loaded 1et4.mat
Loaded 7msf.mat

File = 

1x7 struct array with fields:
    CI
    ClassVersion
    Comment
    Distance
    Edge
    Filename
    HandClass
    Header
    Inter
    Modified
    NT
    NumNT
    Pair

>> File = zGetNTData({'1s72','1j5e','2aw4','2avy','1mfq','1et4','7msf'},1)
Loaded 1s72.mat
Classifying interactions
Found 29027 bases within 15 Angstroms from 3d structure
Bases C 2002(0) and C 1987(0) fall into categories   1.00  23.10 
Found 12235 pairs that are possibly interacting
Classification took 1.13 minutes, or 10796 classifications per minute
Saved 1s72.mat
Loaded 1j5e.mat
Classifying interactions
Found 15003 bases within 15 Angstroms from 3d structure
Error in ==> zAxisAngleRadians at 6
[v,d] = eig(R);                % get eigenvectors and values of R

Error in ==> zPairDiscrepancy at 9
[ax,ang] = zAxisAngleRadians(r);              % rotation angle without a flip

Error in ==> zDistanceToExemplars at 13
    d(k) = zPairDiscrepancy(Exemplar(j,pc),Pair);

Error in ==> zAnalyzePair at 129
    [c,d,h] = zDistanceToExemplars(Exemplar,Pair);

Error in ==> zClassifyPair at 40
    Pair = zAnalyzePair(M1,M2,CL,Exemplar,sh); % analyze and classify pair

Error in ==> zClassifyPairs at 47
  [Pair,s] = zClassifyPair(Ni,Nj,CL,Exemplar);

Error in ==> zGetNTData at 94
      File = zClassifyPairs(File);


>> File = zGetNTData({'1s72','1j5e','2aw4','2avy','1mfq','1et4','7msf'},1)
Loaded 1s72.mat
Classifying interactions
Found 29027 bases within 15 Angstroms from 3d structure
Bases C 2002(0) and C 1987(0) fall into categories   1.00  23.10 
Found 12235 pairs that are possibly interacting
Classification took 1.13 minutes, or 10788 classifications per minute
Saved 1s72.mat
Loaded 1j5e.mat
Classifying interactions
Found 15003 bases within 15 Angstroms from 3d structure
Bases C 1437(A) and C 1465(A) fall into categories   1.00  23.10 
Found  6352 pairs that are possibly interacting
Classification took 0.54 minutes, or 11665 classifications per minute
Saved 1j5e.mat
Loaded 2aw4.mat
Classifying interactions
Found 29431 bases within 15 Angstroms from 3d structure
Found 12185 pairs that are possibly interacting
Classification took 1.13 minutes, or 10779 classifications per minute
Saved 2AW4.mat
Loaded 2avy.mat
Classifying interactions
Found 15141 bases within 15 Angstroms from 3d structure
Bases C  948(A) and C 1234(A) fall into categories   1.00  23.10 
Found  6348 pairs that are possibly interacting
Classification took 0.54 minutes, or 11776 classifications per minute
Saved 2AVY.mat
Loaded 1mfq.mat
Classifying interactions
Found  1016 bases within 15 Angstroms from 3d structure
Found   393 pairs that are possibly interacting
Classification took 0.04 minutes, or 9736 classifications per minute
Saved 1mfq.mat
Loaded 1et4.mat
Classifying interactions
Found  1419 bases within 15 Angstroms from 3d structure
Found   624 pairs that are possibly interacting
Classification took 0.06 minutes, or 10941 classifications per minute
Saved 1et4.mat
Loaded 7msf.mat
Classifying interactions
Found   114 bases within 15 Angstroms from 3d structure
Found    48 pairs that are possibly interacting
Classification took 0.01 minutes, or 3234 classifications per minute
Saved 7msf.mat

File = 

1x7 struct array with fields:
    CI
    ClassVersion
    Comment
    Distance
    Edge
    Filename
    HandClass
    Header
    Inter
    Modified
    NT
    NumNT
    Pair

>> PairViewer(File)

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [ 13] 1

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [1 ] 11

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
No  AA  cSS cis Sugar edge / Sugar edge (second base dominant)  pairs were found.

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [  1] 11

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [11 ] 11

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
No  GG  cSS cis Sugar edge / Sugar edge (second base dominant)  pairs were found.

1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
Enter paircode(s) to view        [ 11] 1

  0  Display all categories together
  1  cWW cis Watson-Crick / Watson-Crick
  2  tWW trans Watson-Crick / Watson-Crick
  3  cWH cis Watson-Crick / Hoogsteen
 -3  cHW cis Hoogsteen / Watson-Crick
  4  tWH trans Watson-Crick / Hoogsteen
 -4  tHW trans Hoogsteen / Watson-Crick
  5  cWS cis Watson-Crick / Sugar edge
 -5  cSW cis Sugar edge / Watson-Crick
  6  tWS trans Watson-Crick / Sugar edge
 -6  tSW trans Sugar edge / Watson-Crick
  7  cHH cis Hoogsteen / Hoogsteen
  8  tHH trans Hoogsteen / Hoogsteen
  9  cHS cis Hoogsteen / Sugar edge
 -9  cSH cis Sugar edge / Hoogsteen
 10  tHS trans Hoogsteen / Sugar edge
-10  tSH trans Sugar edge / Hoogsteen
 11  cSS cis Sugar edge / Sugar edge (second base dominant)
-11  cSS cis Sugar edge / Sugar edge (first base dominant)
 12  tSS trans Sugar edge / Sugar edge (second base dominant)
-12  tSS trans Sugar edge / Sugar edge (first base dominant)
 14  mis miscellaneous hand classification
 21  s35 second base faces up, above first
-21  s53 second base faces up, below first
 22  s33 second base faces down, above first
 23  s55 second base faces down, below first
 25  part of a motif (not implemented)
 26  sugar stacked on base (not implemented)
 30  no interaction
 40  potential stacking above, second base facing up
 41  potential stacking above, second base facing down
 42  potential stacking below, second base facing up
 43  potential stacking below, second base facing down
 44  potential pairing
 45  potential pairing, computer classification 30
 50  in given box
 51  near a pair you specify
Enter category(ies) to view     [11 ] 12

Group 1 - Computer classification matches category
Group 2 - Hand classification matches category
Group 3 - Either cutoff or hand classification matches
Group 4 - Computer classification matches but hand differs
Group 5 - Hand matches but computer differs
Group 6 - Computer classification matches and pair is sequential
Group 7 - Nearest exemplar matches
Group 8 - Nearest exemplar matches but cutoff classification differs
Group 9 - Computer classification matches but nearest exemplar differs
Group 10 - Computer or nearest exemplar match category
Enter group to view        [1] 
Found      9  AA  tSS trans Sugar edge / Sugar edge (second base dominant)  pairs.

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 2
       Filename  Nucl1  Nucl2  Disp1  Disp2  Disp3  Norm3   Ang    Gap MinDist Cutoff   Hand   Exem  EDistPairDisc
  1        1s72 A  306 A  340   7.21   1.69   2.50   0.84 204.0   1.58    2.44  12.00   0.00  12.00   8.29   0.00
  2        1s72 A 1233 A 2604   6.91   1.99   2.97   0.74 206.3   2.54    2.85  12.00   0.00  12.00   8.84   0.00
  3        1s72 A 1247 A 1242   7.81   3.02   1.15   0.96 227.1   0.97    2.49  12.00   0.00   2.00   9.33   0.00
  4        1s72 A 1341 A 1313   7.47   3.36   1.19   0.76 220.8   2.71    2.49  12.00   0.00   2.00   8.98   0.00
  5        1s72 A 2019 A 1829   6.57   3.76  -2.78   0.97 222.0  -2.83    2.90  12.00   0.00   2.00   7.87   0.00
  6        2AW4 A  515 A 1260   7.37   1.48  -1.73   0.96 210.6  -0.50    2.20  12.00   0.00  12.00   8.22   0.00
  7        2AW4 A 1978 A 1773   6.32   3.65  -2.73   0.92 222.6  -2.58    2.64  12.00   0.00   2.00   8.14   0.00
  8        2AVY A  572 A   19   7.16   2.40   2.46   0.75 208.0   1.99    2.40  12.00   0.00  12.00   8.98   0.00
  9        1mfq A  173 A  127   7.96   2.53   1.22   0.79 211.3   2.36    2.51  12.00   0.00  12.00   8.56   0.00

Mode 1 - Scatterplots of pairs
Mode 2 - List of pair parameters
Mode 3 - Pair parameter statistics
Mode 4 - View individual pairs and hand classify
Mode 5 - View context for these pairs
Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [1] 