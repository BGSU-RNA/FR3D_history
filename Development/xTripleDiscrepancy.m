
function [Disc,R,MM,CM,A] = xTripleDiscrepancy(File1,Triple1,File2,Triple2)

% if File1 is a text string (filename), load the file

if strcmp(class(File1),'char'),
  File1name = File1;
  File1 = zAddNTData(File1name,2);
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(Triple1),'char'),
  Triple1 = {Triple1};
end

if strcmp(class(Triple1),'cell'),
  Triple1 = File1.NT(zIndexLookup(File1,Triple1));
else
  Triple1 = File1.NT(Triple1);
end

% if File2 is a text string (filename), load the file and display

if strcmp(class(File2),'char'),
  File2name = File2;
  File2 = zAddNTData(File2name,2);
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(Triple2),'char'),
  Triple2 = {Triple2};
end

if strcmp(class(Triple2),'cell'),
  Triple2 = File2.NT(zIndexLookup(File2,Triple2));
else
  Triple2 = File2.NT(Triple2);
end

if length(Triple1) ~= length(Triple2),           % motif sizes must be the same
  Disc = [];
  R    = eye(3);
else

L = length(Triple2);

A = zeros(L,1);                             % rotation angles for bases

  LocationWeight = ones(1,L);
  AngleWeight    = ones(1,L);

% ------------------------------ Calculate discrepancy

  Triple1Centers = [];
  for i = 1:L,
    Triple1Centers = [Triple1Centers; Triple1(i).Fit(1,:)];
  end

  Triple1Centers


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

end

R = R';
