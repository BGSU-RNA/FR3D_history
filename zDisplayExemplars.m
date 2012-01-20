% zDisplayExamplars(Paircode,Cateogry) displays the best known
% representatives for interactions involving pairs with the given Paircode
% and interaction Category

% zDisplayExemplars(Paircode,Category), where Paircode and Category can be
% vectors, will loop through all possible Paircode, Category pairs from the
% two vectors.  If a certain category has subcategories, like 1, 1.1, and
% 1.2, it will loop through all of those.

% Here are some ways to run the program:

% zDisplayExemplars(7,15:18) paircode 7, categories 15 to 18 (stacking)
% zDisplayExemplars(1:16,1)   all paircodes, category 1
% zDisplayExemplars(1:16,-12:18) all paircodes, all categories

function [void] = zDisplayExemplars(Paircode,Category)

% load exemplars -------------------------------------

  load('PairExemplars','Exemplar');

% specify parameters for viewing -------------------------------------------

  ViewParam.Mode      = 1; 
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

for pc = 1:length(Paircode),
 for ca = 1:length(Category),
   row = find(fix(cat(1,Exemplar(:,Paircode(pc)).Class)) == Category(ca));
   for su = 1:length(row),

     E = Exemplar(row(su),Paircode(pc));
     if ~isempty(E),
       cla
       zPlotExemplar(E,ViewParam);

       axis equal

       switch fix(E.Class),
         case {1, 2},  axis([-2 10 -2 12 -5 5]);
         case 15,      axis([-5 3 -3 5 -5 5]);
         case 16,      axis([-4 4 -3 5 -3 5]);
         otherwise,    axis([-6 10 -6 10 -6 10]);
       end

       Title = [E.NT1.Base E.NT2.Base ' ' num2str(E.Class) ' ' strrep(E.Filename,'_','\_') ' '];
       Title = [Title E.NT1.Base E.NT1.Number '-' E.NT2.Base E.NT2.Number];
       title(Title);

       fprintf('Press a key to go on\n');
       pause

      end
    end
  end
end
