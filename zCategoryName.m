% zCategoryName(a) returns the category numbered a

function [n] = zCategoryName(a)

switch fix(a),
  case   0, n = 'no classification'
  case   1, n = 'cis Watson-Crick / Watson-Crick';
  case   2, n = 'trans Watson-Crick / Watson-Crick';
  case   3, n = 'cis Watson-Crick / Hoogsteen';
  case  -3, n = 'cis Hoogsteen / Watson-Crick';
  case   4, n = 'trans Watson-Crick / Hoogsteen';
  case  -4, n = 'trans Hoogsteen / Watson-Crick';
  case   5, n = 'cis Watson-Crick / Sugar edge';
  case  -5, n = 'cis Sugar edge / Watson-Crick';
  case   6, n = 'trans Watson-Crick / Sugar edge';
  case  -6, n = 'trans Sugar edge / Watson-Crick';
  case   7, n = 'cis Hoogsteen / Hoogsteen';
  case   8, n = 'trans Hoogsteen / Hoogsteen';
  case   9, n = 'cis Hoogsteen / Sugar edge';
  case  -9, n = 'cis Sugar edge / Hoogsteen';
  case  10, n = 'trans Hoogsteen / Sugar edge';
  case -10, n = 'trans Sugar edge / Hoogsteen';
  case  11, n = 'cis Sugar edge / Sugar edge (second base dominant)';
  case -11, n = 'cis Sugar edge / Sugar edge (first base dominant)';
  case  12, n = 'trans Sugar edge / Sugar edge (second base dominant)';
  case -12, n = 'trans Sugar edge / Sugar edge (first base dominant)';
  case  13, n = 'bifurcated cis Watson-Crick / Watson-Crick';
  case  14, n = 'miscellaneous hand classification';
  case  15, n = 'second base faces up, above first';
  case  16, n = 'second base faces down, above first';
  case  17, n = 'second base faces up, below first';
  case  18, n = 'second base faces down, below first';
  case  25, n = 'part of a motif';
  case  26, n = 'sugar stacked on base';
  case  30, n = 'no interaction';
  case  40, n = 'potential stacking above, second base facing up';
  case  41, n = 'potential stacking above, second base facing down';
  case  42, n = 'potential stacking below, second base facing up';
  case  43, n = 'potential stacking below, second base facing down';
  case  44, n = 'potential pairing';
  case  45, n = 'potential pairing, computer classification 30';
  case  50, n = 'in given box';
  case  51, n = 'near a pair you specify';
  otherwise, n = 'unknown';
end
