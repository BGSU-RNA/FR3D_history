% xGetEdgeNums parses the entries in the GUI concerning base-base interactions

function [ReqEdge,ExEdge,OKPairs,ExPairs,BP1,BP2,EBP1,EBP2,Flank,Range] = xGetEdgeNums(str)

ReqEdge = [];
ExEdge  = [];
OKPairs  = [];
ExPairs  = [];
BP1      = [];
BP2      = [];
EBP1     = [];
EBP2     = [];
Flank    = [];
Range    = [];


EdgeStr{1} = 'cWW';
BPequiv{1} = [1 -1];

EdgeStr{2} = 'tWW';
BPequiv{2} = [2 -2];

EdgeStr{3} = 'cWH';
BPequiv{3} = [3];

EdgeStr{4} = 'cHW';
BPequiv{4} = [-3];

EdgeStr{5} = 'tWH';
BPequiv{5} = [4];

EdgeStr{6} = 'tHW';
BPequiv{6} = [-4];

EdgeStr{7} = 'cWS';
BPequiv{7} = [5];

EdgeStr{8} = 'cSW';
BPequiv{8} = [-5];

EdgeStr{9} = 'tWS';
BPequiv{9} = [6];

EdgeStr{10} = 'tSW';
BPequiv{10} = [-6];

EdgeStr{11} = 'cHH';
BPequiv{11} = [7 -7];

EdgeStr{12} = 'tHH';
BPequiv{12} = [8 -8];

EdgeStr{13} = 'cHS';
BPequiv{13} = [9];

EdgeStr{14} = 'cSH';
BPequiv{14} = [-9];

EdgeStr{15} = 'tHS';
BPequiv{15} = [10];

EdgeStr{16} = 'tSH';
BPequiv{16} = [-10];

EdgeStr{17} = 'cSS';
BPequiv{17} = [11 -11];

EdgeStr{18} = 'cSs';
BPequiv{18} = [11];

EdgeStr{19} = 'csS';
BPequiv{19} = [-11];

EdgeStr{20} = 'tSS';
BPequiv{20} = [12 -12];

EdgeStr{21} = 'tSs';
BPequiv{21} = [12];

EdgeStr{22} = 'tsS';
BPequiv{22} = [-12];

EdgeStr{23} = 'sP';
BPequiv{23} = [21 -21];

EdgeStr{24} = 'sA';
BPequiv{24} = [22 -22 23 -23];

EdgeStr{25} = 's35';
BPequiv{25} = [21];

EdgeStr{26} = 's53';
BPequiv{26} = [-21];

EdgeStr{27} = 's33';
BPequiv{27} = [22 -22];

EdgeStr{28} = 's55';
BPequiv{28} = [23 -23];

EdgeStr{29} = 'NP';            % close enough to interact,
                               % but not classified or near a class
BPequiv{29} = [30 -30];

EdgeStr{30} = 'Pair';
BPequiv{30} = [-12:-1 1:12];

EdgeStr{31} = 'Stack';
BPequiv{31} = [21:23 -23:-21];

EdgeStr{32} = 'cis';
BPequiv{32} = [1 3 5 7 9 11 -1 -3 -5 -7 -9 -11];

EdgeStr{33} = 'trans';
BPequiv{33} = [2 4 6 8 10 12 -2 -4 -6 -8 -10 -12];

EdgeStr{34} = 'bif';
BPequiv{34} = [13 -13];

BPStr{1}    = 'BP';
basephoscode{1}  = 1:15;

BPStr{2}    = 'B1P';
basephoscode{2}  = [1];

BPStr{3}    = 'B2P';
basephoscode{3}  = [2];

BPStr{4}    = 'B3P';
basephoscode{4}  = [3];

BPStr{5}    = 'B4P';
basephoscode{5}  = [4];

BPStr{6}    = 'B5P';
basephoscode{6}  = [5];

BPStr{7}    = 'B6P';
basephoscode{7}  = [6];

BPStr{8}    = 'B7P';
basephoscode{8}  = [7];

BPStr{8}    = 'B9P';
basephoscode{8}  = [9];

BPStr{9}    = 'B10P';
basephoscode{9}  = [10];

BPStr{10}    = 'B11P';
basephoscode{10}  = [11];

BPStr{21}    = 'PB';
basephoscode{21}  = [1 2 3 4];

BPStr{22}    = 'P1B';
basephoscode{22}  = [1];

BPStr{23}    = 'P2B';
basephoscode{23}  = [2];

BPStr{24}    = 'P3B';
basephoscode{24}  = [3];

BPStr{25}    = 'P4B';
basephoscode{25}  = [4];

BPequiv{100} = [];

% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

Pairs = {'AA' 'CA' 'GA' 'UA' 'AC' 'CC' 'GC' 'UC' 'AG' 'CG' 'GG' 'UG' 'AU' 'CU' 'GU' 'UU' ''};

str    = strrep(str,'any','');        % no restriction
  
str    = strrep(str,';',',');       % replace delims by commas
str    = strrep(str,':',',');       % replace delims by commas
str    = strrep(str,'|',',');       % replace delims by commas
str    = strrep(str,' ',',');       % replace delims by commas

while strfind(str,',,'),
  str = strrep(str,',,',',');
end

if str(end) == ',',
  str = str(1:(end-1));
end

if str(1) == ',',
  str = str(2:end);
end

commas = strfind(str,',');               % find locations of commas

if isempty(str),
  lim = [];
else
  lim    = [0 commas length(str)+1];       % locations of edge specifications
end

for i=1:length(lim)-1                    % loop through tokens
    NewStr = str(lim(i)+1:lim(i+1)-1);   % extract next token
    if NewStr(1) == '~',                 % opposite of usual sense
      NewStr = NewStr(2:length(NewStr));
      Reverse = 1;
    else
      Reverse = 0;
    end
    if (NewStr(1) == 'n') | (NewStr(1) == 'N'),  % near
      NewStr = NewStr(2:length(NewStr));
      Near = 1;
    else
      Near = 0;
    end

    EdgeNum = str2num(NewStr);               % convert to number
    PairCode = [];                           % default; nothing
    newBP1   = [];
    newBP2   = [];

    if isempty(EdgeNum)                      % NewStr IS a string
      edg  = find(strcmp(EdgeStr,NewStr));   % case sensitive
      edgi = find(strcmpi(EdgeStr,NewStr));  % case insensitive
      pair = find(strcmpi(Pairs,NewStr));    % case insensitive
      basephos = find(strcmpi(BPStr,NewStr));% case insensitive

      EdgeCode = 99;                         % default; nothing

      if ~isempty(edg),
        EdgeCode = edg(1);
      elseif ~isempty(edgi),
        EdgeCode = edgi(1);
      elseif ~isempty(pair),
        PairCode = pair(1);
      elseif ~isempty(basephos),
        if basephos(1) < 20,
          newBP1 = basephoscode{basephos(1)};
        else
          newBP2 = basephoscode{basephos(1)};          
        end
      end

      if strcmpi(NewStr,'flank'),
        Flank = 1 - Reverse;
      elseif strcmpi(NewStr,'local'),
        Range = [1 10];
      elseif strcmpi(NewStr,'nested'),
        Range = [1 1];
      elseif strcmpi(NewStr,'long-range'),
        Range = [11 Inf];
      elseif strcmpi(NewStr,'LR'),
        Range = [11 Inf];
      end

      EdgeNum = BPequiv{EdgeCode};

      if Near > 0,
        NearEdgeNum = [];
        for j = 1:length(EdgeNum),
          e = EdgeNum(j);
          NearEdgeNum = [NearEdgeNum sign(e) * (100 + abs(e))];
        end
        EdgeNum = NearEdgeNum;

        newBP1 = newBP1 + 100;
        newBP2 = newBP2 + 100;
      end
    end

    if Reverse == 0,
      ReqEdge = [ReqEdge EdgeNum];
      OKPairs = [OKPairs PairCode];
      BP1     = [BP1 newBP1];
      BP2     = [BP2 newBP2];
    else
      ExEdge  = [ExEdge EdgeNum];
      ExPairs = [ExPairs PairCode];
      EBP1    = [EBP1 newBP1];
      EBP2    = [EBP2 newBP2];
    end
end
