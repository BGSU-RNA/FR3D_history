function [n] = zGroupNames(group)

switch group,
  case 1, n = 'Cutoff classification matches category';
  case 2, n = 'Hand classification matches category';
  case 3, n = 'Either cutoff or hand classification matches'; 
  case 4, n = 'Cutoff classification matches but hand differs'; 
  case 5, n = 'Hand matches but computer differs';
  case 6, n = 'Cutoff classification matches and pair is sequential';
  case 7, n = 'Nearest exemplar matches';
  case 8, n = 'Nearest exemplar matches but cutoff classification differs';
  case 9, n = 'Cutoff classification matches but nearest exemplar differs';
  case 10, n ='Cutoff or nearest exemplar match';
  case 11, n = 'Berman class matches';
  case 12, n = 'Berman class matches but cutoff classification differs';
  case 13, n = 'Cutoff classification matches but Berman class differs';
  case 14, n ='Cutoff or Berman class match';
end
