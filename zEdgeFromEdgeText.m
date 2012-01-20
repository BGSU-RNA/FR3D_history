% zEdgeFromEdgeText converts a text description such as 'cWW' to a number

function [Edge] = zEdgeFromEdgeText(e)

  e = strrep(e,' ','');

  switch e
    case 'cWW', Edge = 1;
    case 'tWW', Edge = 2;
    case 'cWH', Edge = 3;
    case 'tWH', Edge = 4;
    case 'cWS', Edge = 5;
    case 'tWS', Edge = 6;
    case 'cHH', Edge = 7;
    case 'tHH', Edge = 8;
    case 'cHS', Edge = 9;
    case 'tHS', Edge = 10;
    case 'cSS', Edge = 11;
    case 'tSS', Edge = 12;
    case 'bif', Edge = 13;
    case 'cHW', Edge = -3;
    case 'tHW', Edge = -4;
    case 'cSW', Edge = -5;
    case 'tSW', Edge = -6;
    case 'cSH', Edge = -9;
    case 'tSH', Edge = -10;
    otherwise,  
      Edge = 0;
      fprintf('Unknown edge text in zEdgeFromEdgeText\n');
  end

