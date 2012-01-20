% zBasePhosphateText(e) converts internal code e for a base-phosphate interaction to text for human use

function [E] = zBasePhosphateText(e)

E = [];

% The following vector converts internal FR3D codes to Base-phosphate categories0 to 9.  To change the categories, change here and in xGetEdgeNums.m

BPCat = [2 6 7 0 6 7 8 9 0 1 3 4 5 0 5 9 0 8 4];  % updated 8-19-2008

for i=1:length(e),
  if e(i) > 100,
    E = [E 'n'];
    a = e(i) - 100;
  elseif e(i) < -100,
    E = [E 'n'];
    a = e(i) + 100;
  else
    E = [E ' '];
    a = e(i);
  end

  if a == 18,
    E = [E '8bBP'];
  elseif a == 19,
    E = [E '4bBP'];
  elseif a > 0,
    E = [E num2str(BPCat(a)) 'BP'];
  elseif a == -18,
    E = [E '8bPB'];
  elseif a == -19,
    E = [E '4bPB'];
  elseif a < 0,
    E = [E num2str(BPCat(-a)) 'PB'];
  else
    E = [E '  -'];
  end
end
