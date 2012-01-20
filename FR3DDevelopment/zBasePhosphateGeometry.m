
function [Distance, Score, MinBPh] = zBasePhosphateGeometry(BPh)

centeroption = 1;

% When BPh = 7 or BPh = 12, two oxygens interact simultaneously with two H's.
% When BPh = 18 or BPh = 19, one oxygen interacts with two hydrogens,
%  18 is called 7BPh, 19 is called 4BPh

if BPh >= 1 && BPh <= 19,

% ------------------------------------- Calculate optimal oxygen locations

zStandardBases

for Code = 1:4,

  L = Lim(2,Code);
  Q = StandardLoc(1:L,:,Code);        % locations of atoms in ideal geometry

  if centeroption == 1,
    M = Lim(1,Code);
    Q = Q - ones(L,1)*mean(Q(1:M,:));
  end

    switch Code
      case 1,                         % Base A
              h   = [11 12 14 15];    % rows of the base hydrogens
              hn  = {'H2','H8','1H6','2H6'}; % names of the base hydrogens
              m   = [ 9  7  6  6];    % rows of the corresponding massive atoms
              e   = [ 1  4  2  3];    % code for location of the interaction
      case 2,                         % Base C
              h   = [10 11 12 13];
              hn  = {'H6','H5','1H4','2H4'}; % names of the base hydrogens
              m   = [ 7  8  6  6];
              e   = [ 9  8  6  5];
      case 3,                         % Base G
              h   = [12 13 15 16];
              hn  = {'H1','H8','1H2','2H2'}; % names of the base hydrogens
              m   = [ 4  7 11 11];
              e   = [13 14 10 11];
      case 4,                         % Base U
              h   = [ 9 11 12];
              hn  = {'H5','H3','H6'}; % names of the base hydrogens
              m   = [ 8  4  7];
              e   = [16 15 17];
    end

  for z = 1:length(h),                % loop through hydrogens
    dx = Q(h(z),1)-Q(m(z),1);
    dy = Q(h(z),2)-Q(m(z),2);

    d = [dx dy];
    d = d / norm(d);

    if AtomNames{m(z),Code}(1) == 'N',
      r = 2.9;
    elseif AtomNames{m(z),Code}(1) == 'C',
      r = 3.3;
    else
      fprintf('Unrecognized heavy atom\n');
      Code 
      m(z)
      AtomNames{m(z),Code}
      r = 6;
    end

    OxygenLocations(e(z),:) = Q(m(z),1:2) + r*d;

  end
end



BPhCodes{1} = [1 2 3 4];                     % codes for A
BPhCodes{2} = [5 6 7 8 9];                   % codes for C, omit 18
BPhCodes{3} = [10 11 12 13 14];              % codes for G, omit 19
BPhCodes{4} = [15 16 17];                    % codes for U


OxygenLocations( 7,:) = [1000 1000];         % don't allow these to be close
OxygenLocations(12,:) = [1000 1000];

OxygenLocations(18,:) = (OxygenLocations( 6,:)+OxygenLocations( 8,:))/2;      % average these
OxygenLocations(19,:) = (OxygenLocations(11,:)+OxygenLocations(13,:))/2;      % average these

BPhDist = zDistance(OxygenLocations,OxygenLocations);

if any(BPh == [7 12]),
  for c = 1:4,
    [Distance(c) i] = min([BPhDist(BPh-1,BPhCodes{c}) BPhDist(BPh+1,BPhCodes{c})]);
    if i < length(BPhCodes{c}),
      MinBPh(c) = BPhCodes{c}(i);
    else
      MinBPh(c) = BPhCodes{c}(i-length(BPhCodes{c}));
    end
  end
else
  for c = 1:4,
    [Distance(c) i] = min(BPhDist(BPh,BPhCodes{c}));
    MinBPh(c) = BPhCodes{c}(i);
  end
end

Score = 1./(1+2*Distance.^2);                  % convert to probabilities
Score = Score / sum(Score);                    % normalize

Lett = 'ACGU';

if 0 > 1,
  BPhCodes{1} = [1 2 3 4];                     % codes for A
  BPhCodes{2} = [5 6 7 8 9 18];                  % codes for C
  BPhCodes{3} = [10 11 12 13 14 19];           % codes for G
  BPhCodes{4} = [15 16 17];                    % codes for U
  for BPh = 1:19,
    for c = 1:4,
      if any(BPh == BPhCodes{c}),
        Letter = Lett(c);
      end
    end
    [D,S] = zBasePhosphateGeometry(BPh);
    fprintf('BPh code %d is %s made by %s\n', BPh, zBasePhosphateText(BPh), Letter);
    D
    S
  end
end

else

Distance = [];
Score = [];

end
