
function [Disc,R,MM,CM,A] = xTripleDiscrepancy(File1,i1,File2,i2)

% if File1 is a text string (filename), load the file

if strcmp(class(File1),'char'),
  File1name = File1;
  File1 = zAddNTData(File1name,2);
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(i1),'char'),
  i1 = {i1};
end

if strcmp(class(i1),'cell'),
  i1 = File1.NT(zIndexLookup(File1,i1));
else
  i1 = File1.NT(i1);
end

% if File2 is a text string (filename), load the file and display

if strcmp(class(File2),'char'),
  File2name = File2;
  File2 = zAddNTData(File2name,2);
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(i2),'char'),
  i2 = {i2};
end

if strcmp(class(i2),'cell'),
  i2 = File2.NT(zIndexLookup(File2,i2));
else
  i2 = File2.NT(i2);
end

if length(i1) ~= length(i2),           % motif sizes must be the same
  Disc = [];
  R    = eye(3);
else

L = length(i2);

A = zeros(L,1);                             % rotation angles for bases

  LocationWeight = ones(1,L);
  AngleWeight    = ones(1,L);

% ------------------------------ Calculate discrepancy


  ModelCenters = cat(1,File1.NT(i1).Fit(1,:));

  ModelWeightedCenter = LocationWeight * ModelCenters / L;
  MCC                 = ModelCenters-ones(L,1)*ModelWeightedCenter;

  CandiCenters = cat(1,Cand.Center);
  CMean = LocationWeight * CandiCenters / L;
  CC = CandiCenters - ones(L,1) * CMean;  % subtract mean
  
  R = zBestRotation(CC, diag(LocationWeight)*MCC);      % candidate onto model
  
  S = LocationWeight * sum(((MCC - CC*R').^2)')';  % distances between centers

  n = 1;                                    % nucleotide number for angles
  v = 4 * AngleWeight.^2;                   % precompute a little
  
  while (n <= L),
    angbytwo = acos(min(1,sqrt(trace(R*Cand(n).Rot*(Model(n).Rot)')+1)/2));
    A(n) = angbytwo;
    S   = S + (angbytwo^2)*v(n);
    n   = n + 1;
  end

  Disc = sqrt(S)/L;

  MM = ModelWeightedCenter;
  CM = CMean;


R = R';
