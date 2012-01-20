
function [E] = zBasePhosphateText(e)

E = [];

for i=1:length(e),
  if e(i) > 100,
    E = [E 'nB' num2str(e(i)-100)];
  elseif e(i) > 0,
    E = [E ' B' num2str(e(i))];
  elseif e(i) < -100,
    E = [E 'nP' num2str(e(i)+100)];
  elseif e(i) < 0,
    E = [E ' P' num2str(e(i))];
  else
    E = [E '   -'];
  end
end

return

for i=1:length(e),
  switch fix(e(i))
    case    1,    E = [E ' B1P'];
    case    2,    E = [E ' B2P'];
    case    3,    E = [E ' B3P'];
    case    4,    E = [E ' B4P'];
    case  101,    E = [E 'nB1P'];
    case  102,    E = [E 'nB2P'];
    case  103,    E = [E 'nB3P'];
    case  104,    E = [E 'nB4P'];
    case   -1,    E = [E ' P1B'];
    case   -2,    E = [E ' P2B'];
    case   -3,    E = [E ' P3B'];
    case   -4,    E = [E ' P4B'];
    case -101,    E = [E 'nP1B'];
    case -102,    E = [E 'nP2B'];
    case -103,    E = [E 'nP3B'];
    case -104,    E = [E 'nP4B'];
    otherwise E = [E '   -'];
  end
end