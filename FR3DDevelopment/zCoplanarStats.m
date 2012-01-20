
% File = zAddNTData('1s72');

Stats = [];

for f = 1:length(File),
  E = File(f).Edge;
  [i,j] = find( (abs(E) > 0) .* (abs(E) < 14) );
  for k = 1:length(i),
    N1 = File(f).NT(i(k));
    N2 = File(f).NT(j(k));


  Pair.PlaneAng = acos(abs(N1.Rot(:,3)'*N2.Rot(:,3)))*57.29577951308232; 
                                             % angle between planes

  Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

  d = zDistance(N2.Fit(1:Lim(2,N2.Code),:), N1.Center); 
                                           % distances to base 1 center
  [y,m] = min(d);                          % identify the closest atom
  m = m(1);                                % in case of a tie, use the first
  Pair.Gap = N1.Rot(:,3)'*(N2.Fit(m,:)-N1.Center)';% height above plane of 1

  Pair.MinDist = min(min(zDistance(N1.Fit,N2.Fit)));

  % ------------------------ check for coplanarity

  Pair.Coplanar = 0;                      % default value, not coplanar

  % Criteria for possibly being coplanar:
  %   Pair.Gap must be < 2 Angstroms
  %   Pair.MinDist must be < 4.5 Angstroms
  %   Angle between center-center vector and each normal must be > 60 degrees
  %   Angle between normal vectors must be < 45 degrees


    v  = N1.Center - N2.Center;           % vector from center to center
    v  = v / norm(v);                     % normalize

    dot1 = abs(v * N1.Rot(:,3));          % to calculate angle: v and normal
    dot2 = abs(v * N2.Rot(:,3));
    dot3 = abs(N1.Rot(:,3)' * N2.Rot(:,3));

    yy = 0.5;                             % cos(60) = 0.5
    yyy = 1/sqrt(2);                      % cos(45) = 1/sqrt(2)

      d = zDistance(N1.Fit(1:Lim(2,N1.Code),:), N2.Center); 
                                           % distances to base 2 center
      [y,m] = min(d);                      % identify the closest atom
      m = m(1);                            % in case of a tie, use the first
      Gap2 = N2.Rot(:,3)'*(N1.Fit(m,:)-N2.Center)';% height above plane of 1

        Pair.Coplanar = min([(2-abs(Pair.Gap))/2 (2-abs(Gap2))/2 (yy-dot1)/yy (yy-dot2)/yy (yyy-dot3)/yyy min(1,4.5-Pair.MinDist)]);


%      Stats = [Stats; [(2-abs(Pair.Gap))/2 (2-abs(Gap2))/2 (yy-dot1)/yy (yy-dot2)/yy (dot3-yyy)/yyy min(1,4.5-Pair.MinDist)]];

      Stats = [Stats; [abs(Pair.Gap) abs(Gap2) dot1 dot2 dot3 Pair.MinDist]];

  end
end

names{1} = 'abs(Pair.Gap)';
names{2} = 'abs(Gap2)';
names{3} = 'dot1';
names{4} = 'dot2';
names{5} = 'dot3';
names{6} = 'Pair.MinDist';

for v = 1:6,
  figure(v)
  clf
  if v == 6,
    i = find(Stats(:,v)<4);
  else
    i = 1:length(Stats(:,1));
  end
  hist(Stats(i,v),30)
  title(names{v});
  hold on
  plot(quantile(Stats(i,v),0.10),0,'r*');
  plot(quantile(Stats(i,v),0.90),0,'r*');
  plot(quantile(Stats(i,v),0.50),0,'m*');
  plot(quantile(Stats(i,v),0.70),0,'g*');
  plot(quantile(Stats(i,v),0.30),0,'g*');

  if v ~= 5,
    fprintf('1+(%s-%7.4f)*(%7.4f)\n', names{v}, quantile(Stats(i,v),0.70), 0.5/(quantile(Stats(i,v),0.70)-quantile(Stats(i,v),0.90)));

    fprintf('%s <= %7.4f\n', names{v}, (-quantile(Stats(i,v),0.70)+2*quantile(Stats(i,v),0.90)));
  else
    fprintf('1+(%s-%7.4f)*(%7.4f)\n', names{v}, quantile(Stats(i,v),0.30), 0.5/(quantile(Stats(i,v),0.30)-quantile(Stats(i,v),0.10)));
    fprintf('%s >= %7.4f\n', names{v}, (-quantile(Stats(i,v),0.30)+2*quantile(Stats(i,v),0.10)));
  end
end
