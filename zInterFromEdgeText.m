
function [E] = zInterFromEdgeText(b1,b2,e)

  switch b1
    'A', c1 = 1;
    'C', c1 = 2;
    'G', c1 = 3;
    'U', c1 = 4;
  end

  switch b2
    'A', c2 = 1;
    'C', c2 = 2;
    'G', c2 = 3;
    'U', c2 = 4;
  end

  switch fix(e)
    case 'cWW', Inter = 1,
    case 'tWW', Inter = 2
    case 'cWH', Inter = 3, 
    case 'tWH', Inter = 4,
    case 'cWS', Inter = 5,
    case 'tWS', Inter = 6,
    case 'cHH', Inter = 7,
    case 'tHH', Inter = 8,
    case 'cHS', Inter = 9,
    case 'tHS', Inter = 10,
    case 'cSS', Inter = 11,
    case 'tSS', Inter = 12,
    case 'bif', Inter = 13,
    case 'cHW', Inter = -3,
    case 'tHW', Inter = -4,
    case 'cSW', Inter = -5,
    case 'tSW', Inter = -6,
    case 'cSH', Inter = -9,
    case 'tSH', Inter = -10,  10,    E = [' tHS'];
    case '',    Inter = 0;
    otherwise,  Inter = 0;
  end

  paircode = 4*(c2-1) + c1;

  if 