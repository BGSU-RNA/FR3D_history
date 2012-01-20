
function [E] = zEdgeText(e)

E = [];

for i=1:length(e),
  switch fix(e(i))
    case    1,    E = [E 'cWw'];
    case    2,    E = [E 'tWw'];
    case    3,    E = [E 'cWH'];
    case    4,    E = [E 'tWH'];
    case    5,    E = [E 'cWS'];
    case    6,    E = [E 'tWS'];
    case    7,    E = [E 'cHH'];
    case    8,    E = [E 'tHH'];
    case    9,    E = [E 'cHS'];
    case   10,    E = [E 'tHS'];
    case   11,    E = [E 'cSs'];
    case   12,    E = [E 'tSs'];
    case   13,    E = [E 'bif'];
    case   14,    E = [E 'Rib'];
    case   -1,    E = [E 'cwW'];
    case   -2,    E = [E 'twW'];
    case   -3,    E = [E 'cHW'];
    case   -4,    E = [E 'tHW'];
    case   -5,    E = [E 'cSW'];
    case   -6,    E = [E 'tSW'];
    case   -7,    E = [E 'cHH'];
    case   -8,    E = [E 'tHH'];
    case   -9,    E = [E 'cSH'];
    case  -10,    E = [E 'tSH'];
    case  -11,    E = [E 'csS'];
    case  -12,    E = [E 'tsS'];
    case  -13,    E = [E 'bif'];
    case   14,    E = [E 'Rib'];
    case   21,    E = [E 's35'];
    case  -21,    E = [E 's53'];
    case   22,    E = [E 's33'];
    case  -22,    E = [E 's33'];
    case   23,    E = [E 's55'];
    case  -23,    E = [E 's55'];
    case  101,    E = [E 'ncWW'];
    case  102,    E = [E 'ntWW'];
    case  103,    E = [E 'ncWH'];
    case  104,    E = [E 'ntWH'];
    case  105,    E = [E 'ncWS'];
    case  106,    E = [E 'ntWS'];
    case  107,    E = [E 'ncHH'];
    case  108,    E = [E 'ntHH'];
    case  109,    E = [E 'ncHS'];
    case  110,    E = [E 'ntHS'];
    case  111,    E = [E 'ncSS'];
    case  112,    E = [E 'ntSS'];
    case  113,    E = [E 'nbif'];
    case  114,    E = [E 'nRib'];
    case -101,    E = [E 'ncWW'];
    case -102,    E = [E 'ntWW'];
    case -103,    E = [E 'ncHW'];
    case -104,    E = [E 'ntHW'];
    case -105,    E = [E 'ncSW'];
    case -106,    E = [E 'ntSW'];
    case -107,    E = [E 'ncHH'];
    case -108,    E = [E 'ntHH'];
    case -109,    E = [E 'ncSH'];
    case -110,    E = [E 'ntSH'];
    case -111,    E = [E 'ncSS'];
    case -112,    E = [E 'ntSS'];
    case -113,    E = [E 'nbif'];
    case -114,    E = [E 'nRib'];
    case  121,    E = [E 'ns35'];
    case -121,    E = [E 'ns53'];
    case  122,    E = [E 'ns33'];
    case -122,    E = [E 'ns33'];
    case  123,    E = [E 'ns55'];
    case -123,    E = [E 'ns55'];
    otherwise E = [E '   -'];
  end

if abs(e(i)) < 29,
  if fix(e(i)) == e(i),
    E = [E ' '];
  else
    d = abs(e(i)) - fix(abs(e(i)));
    d = fix(d*10+0.0001);
    d = max(1,d);
    T = 'abcdefghijkl';
    E = [E T(d)];
  end
end

end

