% pIsoScore returns a 4x4 matrix of scores depending on the specified
% interaction category Edge and the specified Pair having this interaction.
% If an optional 4th argument is used, it interprets this as the deletion
% probability and returns a 1 by 17 vector of normalized probabilities.
% The variable method tells what method to use to calculate scores.
% method = 1;                 Use Isosteric subclasses according to LSW 2002
% method = 2;                 Use IDIs mapped into probabilities
% method = 3;                 Use IDIs times frequencies
% method = 4;                 Use IDI with classes weighted 85, 10, 4, 1

% Note to add:  AG water inserted is self-isosteric.  Call this class 14???

function [S] = pIsoScore(Class,Code1,Code2,method,ExemplarIDI,ExemplarFreq)

if nargin < 4,
  method = 2;
end

if nargin < 5 && method > 1,
  load PairExemplars
end

if strcmp(class(Class),'char'),
  Class = xGetEdgeNums(Class);
end

if strcmp(class(Code1),'char'),
  Code1 = pLtoN(Code1);
end

if strcmp(class(Code2),'char'),
  Code2 = pLtoN(Code2);
end


S = ones(4,4)/16;


if abs(Class) > 13,
  fprintf('Class %8.2f is too large for pIsoScore.  Returning a uniform matrix.\n',Class);
else
  switch method
    case 1,  S = pIsoScore1(Class,Code1,Code2);
    case 2,  S = pIsoScore2(Class,Code1,Code2,ExemplarIDI);
    case 3,  S = pIsoScore3(Class,Code1,Code2,ExemplarIDI,ExemplarFreq);
    case 4,  S = pIsoScore4(Class,Code1,Code2,ExemplarIDI);
    otherwise, fprintf('Invalid method passed to pIsoScore');
  end
end
