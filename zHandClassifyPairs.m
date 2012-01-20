% zHandClassifyPairs displays pairs and allows a reclassification

function [File,SP,ViewParam] = zHandClassifyPairs(File,SP,ViewParam)

load('PairExemplars','Exemplar');

figure(1)

k = 1;
Oldk = 0;

while (k <= length(SP)),
  f = SP(k).Filenum;
  p = SP(k).B1Index;
  q = SP(k).B2Index;
  Pair = File(f).Pair(SP(k).PairIndex);

  if ViewParam.Nearby > 0
    n = ViewParam.Nearby;
    B = find((File(f).Distance(p,:) < n).*(File(f).Distance(p,:) > 0) ...
           .*(File(f).Distance(q,:) < n).*(File(f).Distance(q,:) > 0));
  else
    n = 0;
    B = [];
  end

  if ViewParam.Nearby < 0,             % display sequential neighbors
    n = ViewParam.Nearby;
    B = [(SP(k).B1Index+n:SP(k).B1Index-1) ...
         (SP(k).B1Index+1:SP(k).B1Index-n) ...
         (SP(k).B2Index+n:SP(k).B2Index-1) ...
         (SP(k).B2Index+1:SP(k).B2Index-n)];
    B = max(1,B);
    B = min(B,length(File(f).NT));
  end

  C = [p q B];                     % indices of bases to plot

  if ViewParam.Nearby ~= 0,
    zShowInteractionTable(File(f),sort(C));
  end

  clf

  ViewParam.AtOrigin = 1;
  zDisplayNT(File(f),C,ViewParam);

  R = File(f).NT(C(1)).Rot;                   % Rotation matrix for first base
  S = File(f).NT(C(1)).Fit(1,:);              % Location of glycosidic atom

  zPlotHydrogen(File(f),SP(k).PairIndex,R,S);

  grid on
 
  if ViewParam.Exemplars == 1,
    [c,d,h] = zDistanceToExemplars(Exemplar,Pair);
  fprintf('Distance %5.2f to class %5.2f exemplar (dashed lines)\n',d(1),c(1));
  fprintf('Distance %5.2f to class %5.2f exemplar (dotted lines)\n',d(2),c(2));
  fprintf('Distance %5.2f to class %5.2f exemplar (not shown)\n',d(3),c(3));

    VP = ViewParam;
    VP.LineStyle = '--';
    zPlotExemplar(Exemplar(h(1),Pair.Paircode),VP);
    VP.LineStyle = ':';
    zPlotExemplar(Exemplar(h(2),Pair.Paircode),VP);
  end

  Title = strcat(File(f).NT(p).Base,File(f).NT(p).Number);
  for j=2:length(C),
    Title = strcat(Title,'-',File(f).NT(C(j)).Base,File(f).NT(C(j)).Number);
  end;

  Title = [Title ' CompCat ' num2str(File(f).Pair(SP(k).PairIndex).Class,'%4.2f')];
  Title = [Title ' HandCat ' num2str(SP(k).HandClass,'%4.2f')];

  title(Title);
  axis equal

  if SP(k).CI > 0,                       % if there IS a comment
    xlabel(File(f).Comment(SP(k).CI));
  end
  ylabel(strrep(File(f).Filename,'_','\_'));

  view(ViewParam.az,ViewParam.el);

  axis equal

  if (ViewParam.Nearby == 0),
    switch File(f).Pair(SP(k).PairIndex).Class,
      case {1, 2},  axis([-2 10 -2 12 -5 5]);
      case 15,      axis([-5 3 -3 5 -6 2]);
      case 16,      axis([-4 4 -3 5 -3 5]);
      otherwise,    axis([-6 10 -6 10 -6 10]);
    end
  end

  if k ~= Oldk,
    fprintf('Pair %3d | ',k);
    fprintf('q: quit | b: back | s: sugars | n: nearby | j: jump | e: exemplars | m: modify');
    if (max(Pair.Paircode == [1 6 11 16]) == 1),
      fprintf(' | r: reverse');
    end
    fprintf('\n');
    if length(Pair.Hydrogen) > 0,
      fprintf('Hydrogen bonds: ');
      for i=1:length(Pair.Hydrogen),
        fprintf(' %5.1fdeg %5.2fA |',Pair.Hydrogen(i).Angle,Pair.Hydrogen(i).Distance);
      end
      fprintf('\n');
    end
  end
  Oldk = k;
  inp = lower(input('','s'));
  if isempty(inp),
    k = k + 1;                                   % move to next pair
  else
    if inp(1) == 'b',
      k = max(1,k - 1);
    elseif inp(1) == 'q', 
      k = length(SP) + 1;
    elseif inp(1) == 's', 
      ViewParam.Sugar = 1 - ViewParam.Sugar;
    elseif inp(1) == 'n', 
      if length(inp)>1,
        ViewParam.Nearby = str2num(inp(2:length(inp)));
      elseif ViewParam.Nearby > 0,
        ViewParam.Nearby = 0;
      elseif ViewParam.Nearby == 0,
        ViewParam.Nearby = 8;
      end
    elseif inp(1) == 'j',
      k = abs(str2num(inp(2:length(inp))));
    elseif inp(1) == 'e',
      ViewParam.Exemplars = 1-ViewParam.Exemplars;
    elseif (inp(1) == 'a'),
      N1 = File(f).NT(p);
      N2 = File(f).NT(q);
      RevPair = zAnalyzePair(N1,N2);

      [c,d,h] = zDistanceToExemplars(Exemplar,Pair);
      RevPair.Classes = c(1:3);
      RevPair.Distances = d(1:3);
      RevPair.ExemIndex = h(1:3);

      if (length(Pair.Hydrogen) == 0) & (abs(Pair.Classes(1)) < 14),
        Pair.Hydrogen = zCheckHydrogen(N1,N2,Pair.Classes(1));
      end

      RevPair.Base1Index = p;
      RevPair.Base2Index = q;

      File(f).Inter(p,q)  = RevPair.Class;            % record this interaction
      File(f).Inter(q,p)  = RevPair.Class;

      SP(k).B1Index = p;
      SP(k).B2Index = q;

      zListPairData(RevPair,1);      
      Pair = RevPair;

      File(f).Pair(SP(k).PairIndex) = RevPair;
      File(f).Modified = 0;                          % flag to save .hand
      
    elseif (inp(1) == 'r') & (max(Pair.Paircode == [1 6 11 16]) == 1),
      N1 = File(f).NT(q);
      N2 = File(f).NT(p);
      RevPair = zAnalyzePair(N1,N2);

      [c,d,h] = zDistanceToExemplars(Exemplar,Pair);
      RevPair.Classes = c(1:3);
      RevPair.Distances = d(1:3);
      RevPair.ExemIndex = h(1:3);

      if (length(Pair.Hydrogen) == 0) & (abs(Pair.Classes(1)) < 14),
        Pair.Hydrogen = zCheckHydrogen(N1,N2,Pair.Classes(1));
      end

      RevPair.Base1Index = q;
      RevPair.Base2Index = p;

      File(f).Inter(p,q)  = RevPair.Class;            % record this interaction
      File(f).Inter(q,p)  = RevPair.Class;

      SP(k).B1Index = q;
      SP(k).B2Index = p;

      zListPairData(RevPair,1);      
      Pair = RevPair;

Pair
RevPair

      File(f).Pair(SP(k).PairIndex) = RevPair;
      File(f).Modified = 1;                          % flag to save .hand
    elseif inp(1) == 'm',
      fprintf('Modify current hand classification [%4.2f] ',SP(k).HandClass);
      inp1 = input('','s');

      if SP(k).CI > 0,
        comm = File(f).Comment{SP(k).CI};
      else
        comm = [];
      end
      fprintf('Replace or append (+) current comment [%s] ', comm);
      inp2 = input('','s');
      
      if ~isempty(str2num(inp1)) | ~isempty(inp2),     % something changed
        File(f).Modified = 1;                          % flag to save .hand
        if SP(k).CI > 0,
          m = SP(k).CI;                                % pair has hand class
        else
          m = length(File(f).HandClass)+1;             % new hand pair
          SP(k).CI = m;
          File(f).CI(SP(k).B1Index,SP(k).B2Index) = m;
          File(f).CI(SP(k).B2Index,SP(k).B1Index) = m;
          SP(k).HandClass = 0;                         % default value
          File(f).HandClass(m) = 0;                    % default value
          File(f).Comment{m}   = [];
        end
        
        if ~isempty(str2num(inp1)),
          SP(k).HandClass = str2num(inp1);
          File(f).HandClass(m) = str2num(inp1);
        end

        if ~isempty(inp2),
          if inp2(1) == '+',
            File(f).Comment{m} = [comm ' ' inp2];
          else
            File(f).Comment{m} = inp2;
          end
        end
      end
      k = k + 1;
    end
  end
  [ViewParam.az,ViewParam.el] = view;
end

