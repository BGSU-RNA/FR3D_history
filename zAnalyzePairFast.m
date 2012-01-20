% zAnalyzePairFast(N1,N2,CL) computes distances, angles, and classification
% codes.  Fast means that the calling program passes CL, the class limits.

function [Pair] = zAnalyzePairFast(N1,N2,CL)

  Pair.Paircode = 4*(N2.Code-1) + N1.Code;     % AA is 1, CA is 2, etc.

  Pair.Displ = (N2.Fit(1,:)-N1.Fit(1,:))*N1.Rot; % vector shift (from 1 to 2)

  ro = N1.Rot'*N2.Rot;                       % rotation matrix from 1 to 2
  Pair.Normal = ro(:,3)';                    % normal to second plane

  if ro(3,3) > 0,                            % depending on orientation of 2,
    [ax,ang] = zAxisAngle(ro);              % rotation angle without a flip
  else
    [ax,ang] = zAxisAngle(ro*diag([-1 1 -1])); % flip base 2 first
  end

  Pair.Rot      = ro;
  Pair.RotAx    = ax';
  Pair.Ang      = ang;

  Pair.PlaneAng = acos(abs(N1.Rot(:,3)'*N2.Rot(:,3)))*57.29577951308232; 
                                             % angle between planes

  Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

  d = zDistance(N2.Fit(1:Lim(2,N2.Code),:), N1.Center); 
                                           % distances to base 1 center
  [y,m] = min(d);                          % identify the closest atom
  m = m(1);                                % in case of a tie, use the first
  Pair.Gap = N1.Rot(:,3)'*(N2.Fit(m,:)-N1.Center)';% height above plane of 1

  Pair.MinDist = min(min(zDistance(N1.Fit,N2.Fit)));

  a = zCheckCutoffs(Pair.Displ,Pair.Normal,Pair.Ang,Pair.Gap,CL(:,:,Pair.Paircode));

  % ---------- Notify and remove multiple classifications

  if length(a) > 1,
    fprintf('Bases %1s%5s(%1s) and %1s%5s(%1s) fall into categories ', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain);
    for k=1:length(a),
      fprintf('%6.2f ',a(k));
    end
    fprintf('\n');
    a = a(1);
  end

  % ---------- Calculate hydrogen bonds for base pairing interactions

  if (abs(a) < 14) ,
    Pair.Hydrogen = zCheckHydrogen(N1,N2,a);
  else
    Pair.Hydrogen = [];
  end

  % ---------- Eliminate out of plane interactions and bad hydrogen bonds

  if ((abs(a) < 11) | (abs(a) == 13)) & (abs(Pair.Gap) > 1.0),
    if length(Pair.Hydrogen) > 0,
      if (min(cat(1,Pair.Hydrogen(:).Angle)) < 110) | (abs(Pair.Gap) > 2.0),
        a = 30.1;
      end
    elseif abs(Pair.Gap) > 1.6,
      a = 30.2;
    end
  elseif ((abs(a) == 11) | (abs(a) == 12)) & (abs(Pair.Gap) > 3.0),
                            % disallow immense gaps in cases 11 and 12
    a = 30.3;
  end

  if (abs(a) < 14) & (length(Pair.Hydrogen) > 0),
    if min(cat(1,Pair.Hydrogen(:).Distance)) > 4,
      a = 30.4;
    end
  end


  % ---------- Measure stacking overlap

  SO1 = zStackingOverlap(N1,N2);
  SO2 = zStackingOverlap(N2,N1);

  if (SO1 > 0) & (SO2 > 0),
    Pair.StackingOverlap = (SO1+SO2)/2;
  else
    Pair.StackingOverlap = 0;
  end
  
  % ----------- Is an unclassified pair really stacked?

  if (a == 30) & (Pair.StackingOverlap > 0) & (Pair.MinDist < 4),
    if Pair.Displ(3) > 0,
      if Pair.Normal(3) > 0,
        a = 15;
      else
        a = 16;
      end
    else
      if Pair.Normal(3) > 0,
        a = 17;
      else
        a = 18;
      end
    end
  end

  Pair.Class    = a;

