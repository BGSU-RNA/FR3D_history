% zDisplayExamplar(File,Paircode,Category) displays the best known
% representative of the given Paircode and Category.  It makes a
% scatterplot of all pairs in that category and marks the exemplar
% specially.

function [void] = zDisplayExemplars(File,Paircode,Category)

% load exemplars from previous session -------------------------------------

 load('PairExemplars','Exemplar');

% specify parameters for viewing -------------------------------------------

  ViewParam.Mode      = 1; 
  ViewParam.Color     = 1;
  ViewParam.Normal    = 1;
  ViewParam.ColorAxis = [-12 30];
  ViewParam.SortKeys  = [];
  ViewParam.Nearby    = 0;
  ViewParam.Sugar     = 0;
  ViewParam.Hydrogen  = 1;
  ViewParam.Sort      = 0;
  ViewParam.az        = 51;
  ViewParam.el        = 14;
  ViewParam.LineStyle = '-';

% loop through computer classifications ----------------------

 for row = 1:length(Exemplar(:,Paircode)),

  E = Exemplar(row,Paircode);

  % specify criteria for selection of pairs ----------------------------------

  Param.Paircode = Paircode;
  Param.Category = E.Class;        
  Param.Decimal  = 1;        % 1 - use 1.0 only; 0 - round to 1
  Param.Group    = 1;        % computer classification matches
  Param.Expert   = 0;
  Param.Sequential = 0;
%  Param.Inbox    = [-5 5 -5 5 -5 5 -1.1 1.1 -95 275];  % if category = 50,
                           % only pairs in this box will be selected

  fprintf('Paircode %2d Class %3.2f ', Paircode, E.Class);

  % select pairs using selection criteria -----------------------------------

  SP = zSelectPairs(File,Param);

  f = zFileNumberLookup(File,E.Filename);
  if (length(SP) > 0) & (f > 0),
    figure(1)
    clf
    zPlotExemplar(File,E,ViewParam)

    figure(2)
    clf

    ViewParam.Normal = 0;
    for k = 1:length(SP),                               % Loop through pairs
      p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
      e = p.Displ;

      if E.Displ == p.Displ,
        scatter3(e(1),e(2),e(3),28,2,'filled')
      else
        scatter3(e(1),e(2),e(3),18,1,'filled')
      end

      hold on

      if ViewParam.Normal == 1,
        v = p.Normal/3;                                  % add normal vector
        plot3([e(1) e(1)+v(1)], [e(2) e(2)+v(2)], [e(3) e(3)+v(3)], 'b');
      end
    end

    p = File(SP(1).Filenum).Pair(SP(1).PairIndex);   % use the first pair
    n = File(SP(1).Filenum).NT(p.Base1Index);        % use the first base
    R = n.Rot;                                       % Rotation matrix for first
    S = n.Fit(1,:)';                                 % Location of glycosidic
    L = length(n.Fit(:,1));                          % number of atoms in base
    m.Code = n.Code;
    m.Fit = (R'*(-S*ones(1,L) + n.Fit'))';           % rotated into position
    zPlotOneNT(m,ViewParam);                         % plot this nucleotide

    caxis([1 4]);
    view(2)

    title('First base shown at origin, N1/N9 atom of second base shown by dots, lines are the normal vector of second base');
    xlabel('Perpendicular to glycosidic bond');
    ylabel('Parallel to glycosidic bond');
    zlabel('Vertical with respect to first base');
 
    figure(3)
    clf

    for k = 1:length(SP),                               % Loop through pairs
      p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
      if E.Displ == p.Displ,
        scatter3(p.PlaneAng,p.Ang,p.Normal(3),28,2,'filled')
      else
        scatter3(p.PlaneAng,p.Ang,p.Normal(3),18,1,'filled')
      end
      hold on
    end

    xlabel('Angle between planes');
    ylabel('Appropriate angle of rotation');
    zlabel('Third component of normal vector');
    title('');
    caxis([1 4]);
    grid on
    view(2)

    pause
  end

 end
end
