
% xFindPattern hijacks FR3D to search for a pattern in a field of points

function [Candidates] = xFindPattern(field,pattern,D,Verbose)

if nargin < 1,
  pattern = rand(5,3);                 % pattern to look for
  field = rand(20,3);                  % where to look
  field = [field; pattern];
  D = 0.01;                             % maximum discrepancy
end

if nargin < 4,
  Verbose = 2;
end

if Verbose > 1,
  figure(1)
  clf
  subplot(3,1,1);
  plot(field(:,1),field(:,2),'.k');
  axis([min(field(:,1)) max(field(:,1)) min(field(:,2)) max(field(:,2))]);
  subplot(3,1,2);
  plot(pattern(:,1),pattern(:,2),'.r');
  axis([min(field(:,1)) max(field(:,1)) min(field(:,2)) max(field(:,2))]);
end

File.Distance = zDistance(field);
File.Filename = 'findingpattern';
File.NumNT    = length(field(:,1));
for n = 1:File.NumNT,
  File.NT(n).Code = 1;
end


Query.Geometric = 1;
Query.Diameter = 8;
Query.Distance = zDistance(pattern);
Query.NumNT = length(pattern(:,1));
Query.LocWeight = ones(1,Query.NumNT);
Query.AngleWeight = ones(1,Query.NumNT);
Query.DiscCutoff = D;
Query.RelCutoff = Query.DiscCutoff;
Query.SSCutoff  =(Query.NumNT^2)*(Query.RelCutoff^2)*cumsum(Query.LocWeight);


Candidates = xFindCandidates(File,Query,2);

L = length(Candidates(:,1));                     % number of candidates

tic

if Verbose > 0,
  fprintf('Calculating discrepancies: ');
end

MCC = pattern - ones(Query.NumNT,1)*mean(pattern);

Discrepancy = [];

for i=1:L,
  C  = field(Candidates(i,1:Query.NumNT),:);         % coordinates of candidate i
  CC = C - ones(Query.NumNT,1)*mean(C);

  R = zBestRotation(CC, diag(Query.LocWeight)*MCC);      % candidate onto model
  
  S = Query.LocWeight * sum(((MCC - CC*R').^2)')';  % distances between centers

  Discrepancy(i) = S;

  if Verbose > 0,
    if (mod(i,round(L/10)) == 0) && (Verbose > 0)
      fprintf(' %d', fix((L-i)*toc/i)); 
      drawnow
    end
  end
end

if Verbose > 0,
  fprintf('\n');
end

Discrepancy = sqrt(Discrepancy)/Query.NumNT;

[y,i]       = sort(Discrepancy);                    % sort by discrepancy
Candidates  = Candidates(i,:);
Discrepancy = Discrepancy(i);

if Verbose > 1,
  fprintf('Calculating discrepancy took        %8.3f seconds\n',toc);
end

if Verbose > 1,
  figure(1)
  subplot(3,1,3);
  for j = 1:length(Candidates(:,1)),
    i = Candidates(j,1:Query.NumNT);
    plot(field(i,1),field(i,2),'.k');
    axis([min(field(:,1)) max(field(:,1)) min(field(:,2)) max(field(:,2))]);
    pause
  end
end
