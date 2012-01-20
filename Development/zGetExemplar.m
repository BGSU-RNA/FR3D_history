
function [NT1,NT2,E] = zGetExemplar(Class,Code1,Code2)

load PairExemplars

if strcmp(class(Class),'char'),
  Class = xGetEdgeNums(Class);
  Class = Class(1);
end

if strcmp(class(Code1),'char'),
  Code1 = pLtoN(Code1);
end

if strcmp(class(Code2),'char'),
  Code2 = pLtoN(Code2);
end

paircode = 4*(Code2-1) + Code1;           % AA is 1, CA is 2, etc.

sw = 1;

switch paircode                           % reverse certain pairs
  case {2, 3, 4, 8, 10, 12},                  % put N2 at the origin
    paircode = 4*(Code1-1) + Code2;
    Class = -Class;
    sw  = -1;                                  % bases in reversed order
end

if any(Class == [-1 -2 -7 -8]),
  Class = -Class;
end

% Paircode list
% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

if any(paircode == [1 6 11 16]),               % identical bases
  if Class < 0,
    Class = -Class;
    sw = -sw;
  end
end

[s,t] = size(Exemplar);

E.Filename = '';
NT1.Code = [];
NT2.Code = [];

for r = 1:s,
  if ~isempty(Exemplar(r,paircode).Filename),
    if Exemplar(r,paircode).Class == Class,
      E = Exemplar(r,paircode);
      if sw > 0,
        NT1 = E.NT1;
        NT2 = E.NT2;
      else
        NT1 = E.NT2;
        NT2 = E.NT1;
      end
    end
  end
end

