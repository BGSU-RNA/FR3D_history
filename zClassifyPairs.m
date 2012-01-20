% zClassifyPairs(File) calculates the rotation matrix, axis, angle, and shift
% between bases in File that are close enough to possibly be interacting, then
% classifies the interaction

function [File] = zClassifyPairs(File)

if isfield(File,'Pair'),
  File = rmfield(File,'Pair');                    % remove previous pair info
end

if File.NumNT > 0,

fprintf('Classifying interactions\n')

t = cputime;

CL = zClassLimits;                              % read ClassLimits matrix

% -------- First screening of base pairs ------------------------------------ 

DistCutoff = 15;
[i,j] = find((File.Distance < DistCutoff).*(File.Distance > 0)); 
                                                % screen by C-C distance
k = find(i<j);                                  % look at each pair only once
i = i(k);                                       % reduce list of indices
j = j(k);                                       % reduce list of indices

fprintf('Found %5d bases within %2d Angstroms from 3d structure\n', length(i), DistCutoff);

% -------- Screen and analyze base pairs ------------------------------------ 

pc = 1;                                         % index for pairs

for k = 1:length(i),                            % loop through possible pairs
  N1 = File.NT(i(k));                           % first nucleotide
  N2 = File.NT(j(k));                           % second nucleotide
  paircode = 4*(N2.Code-1) + N1.Code;           % AA is 1, CA is 2, etc.
  switch paircode
    case {2, 3, 4, 8, 10, 12},                  % put N2 at the origin
      M  = N1;
      N1 = N2;
      N2 = M;
      s  = -1;                                  % bases in reversed order
    otherwise
      s  = 1;                                   % bases in original order
  end

  sh = N1.Rot'*(N2.Fit(1,:)-N1.Fit(1,:))';   % vector shift from 1 to 2
                                             % between glycosidic atoms,
                                             % relative to the plane of base 1

  ci = File.CI(i(k),j(k));                   % Comment index


  if (abs(sh(3)) < 5) || (ci > 0)            % if hand classified 
                                             % or small vertical shift
    Pair = zAnalyzePairFast(N1,N2,CL);       % analyze and classify pair

    if ci > 0,                                   % Pair is in hand file
      HandClass = File.HandClass(ci);         % Use the hand class
    else
      HandClass = 0;                             % Hand class 0
    end

    if (Pair.Class == 30) & (N1.Code == N2.Code),  % re-analyze AA, CC, ...
      Pair2 = zAnalyzePairFast(N2,N1,CL);
      if (Pair2.Class == 30) & (ci > 0),     % didn't match, but hand class'd
        if HandClass < 0,                    % higher-numbered base at origin
          Pair = Pair2;                      % put in reversed order
          s = -1;
        end
      elseif abs(Pair2.Class) < 15,          % some sort of base pairing
        Pair = Pair2;                        % matched with N2 at origin
        s = -1;
      else
        Pair.Class = Pair2.Class;            % original order, use 2nd class
      end
    end

   Pair.Classes   = zeros(1,3);
   Pair.Distances = zeros(1,3);
   Pair.ExemIndex = zeros(1,3);

   if (length(Pair.Hydrogen) == 0) & (abs(Pair.Classes(1)) < 14),
     Pair.Hydrogen = zCheckHydrogen(N1,N2,Pair.Classes(1));
   end

   File.Inter(i(k),j(k))  = Pair.Class;            % record this interaction
   File.Inter(j(k),i(k))  = Pair.Class;

   if (s == 1),
     Pair.Base1Index = i(k);                       % bases in original order
     Pair.Base2Index = j(k);
   else
     Pair.Base1Index = j(k);                       % bases in reversed order
     Pair.Base2Index = i(k);
   end

   File.Pair(pc) = Pair;

   pc = pc + 1;                                    % increment pair counter

  end
end

if pc > 1,
  A = cat(1, File.Pair(:).Base1Index);
  B = cat(1, File.Pair(:).Base2Index);
  C = min(A,B);

  [y,i] = sort(C);

  File.Pair = File.Pair(i);                          % order by lower base index
else
  File.Pair = [];
end

fprintf('Found %5d pairs that are possibly interacting\n', pc-1);

fprintf('Classification took %4.2f minutes, or %4.0f classifications per minute\n', (cputime-t)/60, 60*(pc-1)/(cputime-t));

else
  File.Pair = [];
end
