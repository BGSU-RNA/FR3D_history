% zContextViewer(File,SP,Param,ViewParam) is in development.  It will
% someday display multiple 3d scatter plots of pairs with context
% information.

function [FigsDone] = zContextViewer(File,SP,Param,ViewParam)

  % --------- Prompt user to choose desired context ---------------

  fprintf('\n');
  fprintf('Enter context to view\n');
  fprintf('1 Pair stacked after current pair\n');
  fprintf('2 Pair stacked before current pair\n');
  fprintf('3 Pairs stacked before and after current pair\n');
  fprintf('Enter context to view        [');
  fprintf(' %1d',Param.Context);
  fprintf('] ');
  inp = input('');
  if ~isempty(inp),
    Param.Context = inp;
  end

  if (length(SP) > 0) & (Param.Sequential == 0),
    zContextAnalyzer(File,SP,Param,ViewParam);
  elseif Param.Sequential == 1,
    zSequentialAnalyzer(File,SP,Param,ViewParam);
  end

return


  ViewParam.Mode  = 1;
  ViewParam = zEnterViewMode(Param,ViewParam);

  SP = zColorPairs(File,SP,Param,ViewParam);

T = cat(1,SP(:).Color);

switch ViewParam.Color,
  case 1, ColorAxis =  [0 8];
          fprintf('Blue - computer matches; Green - hand matches; Orange - both match\n');
  case 2, ColorAxis =  [-12 12];
  case 3, ColorAxis =  [15 22];
  case 4, ColorAxis =  [-12 30];
  case 5, ColorAxis =  [1 16];
  case 8, ColorAxis =  [0 8];
          fprintf('Blue - cutoff matches; Green - exemplar matches; Orange - both match\n');
  case 9, ColorAxis =  [0 4];
  case 10, ColorAxis = [-12 12];
  otherwise, ColorAxis =  [min(T) max(T)+0.01];
end

%-----------------------------------------------------------------------
figure(ViewParam.FigNum)
clf

set(gcf,'Renderer','OpenGL');

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  e = p.Displ;

  scatter3(e(1),e(2),e(3),18,c,'filled')
  hold on


end

if length(Param.Paircode) == 1,             % if only one pair code used
  p = File(SP(1).Filenum).Pair(SP(1).PairIndex);   % use the first pair
  n = File(SP(1).Filenum).NT(p.Base1Index);        % use the first base
  R = n.Rot;                                       % Rotation matrix for first
  S = n.Fit(1,:)';                                 % Location of glycosidic
  L = length(n.Fit(:,1));                          % number of atoms in base
  m.Code = n.Code;
  m.Fit = (R'*(-S*ones(1,L) + n.Fit'))';           % rotated into position
  VP = ViewParam;
  VP.Sugar = 0;
  zPlotOneNT(m,VP);                                 % plot this nucleotide
end

caxis(ColorAxis);
view(2)
axis square

Title =[n.Base ' shown at the origin, N1/N9 atom of second base shown by dots'];
if ViewParam.Normal == 1,
 Title = [Title ', lines are the normal vector of second base'];
end
title(Title);
xlabel('Perpendicular to glycosidic bond');
ylabel('Parallel to glycosidic bond');
zlabel('Vertical with respect to first base');

%----------------------------------------------------------------------
figure(ViewParam.FigNum + 1)
clf

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  scatter3(p.Ang,p.Normal(3),p.Gap,18,c,'filled')
  hold on
end

xlabel('Appropriate angle of rotation');
ylabel('Third component of normal vector');
zlabel('Gap');
title('');
caxis(ColorAxis);
grid on

view(2)

FigsDone = 2;

return

%-----------------------------------------------------------------------
figure(ViewParam.FigNum + 2)
clf

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  scatter3(SP(k).MinDist,p.Displ(3),p.Gap,18,c,'filled')
  hold on
end

xlabel('Minimum distance');
ylabel('Vertical displacement');
zlabel('Gap');
title('');
caxis(ColorAxis);
grid on

FigsDone = 3;

%-----------------------------------------------------------------------
figure(ViewParam.FigNum + 3)
clf

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  scatter3(SP(k).C1pC1p,p.Ang,p.Gap,18,c,'filled')
  hold on
end

xlabel('C1*-C1* distance');
ylabel('Angle of rotation');
zlabel('Gap');
title('');
caxis(ColorAxis);
grid on

FigsDone = 4;

return

%-----------------------------------------------------------------------
figure(ViewParam.FigNum + 4)
clf

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  L = length(p.Hydrogen); 
  if L > 0,
    a2 = p.Hydrogen(min(2,L)).Angle;
    a3 = p.Hydrogen(min(3,L)).Angle;
    scatter3(p.Hydrogen(1).Angle,a2,a3,18,c,'filled')
    hold on
  end
end

xlabel('First hydrogen bond');
ylabel('Second hydrogen bond, if any, otherwise first');
zlabel('Third hydrogen bond, if any, otherwise first');
title('');
caxis(ColorAxis);
grid on

FigsDone = 5;


%-----------------------------------------------------------------------
figure(ViewParam.FigNum + 5)
clf

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  L = length(p.Hydrogen); 
  if L > 0,
    a2 = p.Hydrogen(min(2,L)).Distance;
    a3 = p.Hydrogen(min(3,L)).Distance;
    scatter3(p.Hydrogen(1).Distance,a2,a3,18,c,'filled')
    hold on
  end
end

xlabel('First hydrogen bond length');
ylabel('Second hydrogen bond length, if any, otherwise first');
zlabel('Third hydrogen bond length, if any, otherwise first');
title('');
caxis(ColorAxis);
grid on

FigsDone = 6;


return

%---------------------------------------------------------------------

scatter3(RotAx(G,1),RotAx(G,2),RotAng(G),18,T,'filled')
title('Two components of axis of rotation, vertical is angle')
xlabel('Component 1 of axis of rotation')
ylabel('Component 2 of axis of rotation')
zlabel('Angle of rotation')
axis([-1 1 -1 1 -90 270]);
caxis(ColorAx);

figure(4)
clf
scatter3(MinDist(G),PlaneAng(G),Ang(G),18,T,'filled')
xlabel('Minimum distance');
ylabel('Angle between planes');
zlabel('Appropriate angle of rotation');
caxis(ColorAx);

