% zExamplarTable(Cateogry) displays the best known
% representatives for interactions involving all pairs in
% interaction category(ies) Category

% Here are some ways to run the program:

% zExemplarTable(1:3)

function [void] = zExemplarTable(Category)

% load exemplars -------------------------------------

  load('PairExemplars','Exemplar');

% specify parameters for viewing -------------------------------------------

  ViewParam.Mode      = 1; 
  ViewParam.Normal    = 1;
  ViewParam.ColorAxis = [-12 30];
  ViewParam.SortKeys  = [];
  ViewParam.Nearby    = 0;
  ViewParam.Sugar     = 1;
  ViewParam.ConnectSugar = 0;
  ViewParam.AtOrigin  = 1;
  ViewParam.Hydrogen  = 1;
  ViewParam.Sort      = 0;
  ViewParam.LabelBases= 8;                  % font size

% loop through computer classifications ----------------------

for ca = 1:length(Category),
 figure(fix(Category(ca)))
 clf
 for c1 = 1:4,
  for c2 = 1:4,
   pc = 4*(b2-1)+b1;                         % current paircode

   for r = 1:length(Exemplar(:,1)),          % loop through rows
    E = Exemplar(r,pc);
    if ~isempty(E.Filename),
     if abs(fix(E.Class)) == Category(ca),
       if (E.Class < 0) || ...
          (any(Category(ca) == [1 2 7 8]) && (any(pc == [2 3 4 8 10 12])),


   b1 = c1;
   b2 = c2;


   if any(Category(ca) == [1 2 7 8]),
     


   row = find(abs(fix(cat(1,Exemplar(:,pc).Class))) == Category(ca));
   for su = 1:length(row),

     E = Exemplar(row(su),pc);
     if ~isempty(E.Filename),

       fprintf('%s%s %s %s%s %s Category %3.1f\n',E.NT1.Base,E.NT1.Number,E.Pair.EdgeText,E.NT2.Base,E.NT2.Number,E.Filename,E.Class);
      zListPairData(E.Pair,1);

% display the exemplar pair ---------------------------------------------

       if abs(E.Class - fix(E.Class)) == 0,
         ViewParam.LineStyle = '-';
       elseif abs(E.Class - fix(E.Class)) > 0.29,
         ViewParam.LineStyle = '.';
       elseif abs(E.Class - fix(E.Class)) > 0.19,
         ViewParam.LineStyle = ':';
       elseif abs(E.Class - fix(E.Class)) > 0.09,
         ViewParam.LineStyle = '--';
       end

       subplot(4,4,pc);
       F.NT(1) = E.NT1;
       F.NT(2) = E.NT2;
       F.Filename = E.Filename;
       zDisplayNT(F,[1 2],ViewParam);
       zPlotHydrogenBonds(E.NT1,E.NT2,E.Class,E.NT1.Rot,E.NT1.Fit(1,:));

       view(2)
       grid off
       axis equal

       Title = [E.NT1.Base E.NT2.Base ' ' num2str(E.Class) ' ' strrep(E.Filename,'_','\_') ' '];
       Title = [Title E.NT1.Base E.NT1.Number '-' E.NT2.Base E.NT2.Number];
       CP = norm(E.NT1.Sugar(1,:) - E.NT2.Sugar(1,:));     % c1'-c1' dist
       Title = [Title num2str(CP)];
       title(Title);

%       a = axis;
%       text(0.3*a(1)+0.7*a(2),0.8*a(3)+0.2*a(4),['Count: ' num2str(E.Count)]);
       xlabel(['Count: ' num2str(E.Count)]);

       rotate3d on


      end
    end
   end
  end
end
