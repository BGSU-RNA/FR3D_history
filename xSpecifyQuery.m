% xSpecifyQuery returns the description of a model motif
% The variable Query has several fields, most of which are optional:
%   Query.Description    a useful string, can be long
%   Query.Name           a short string, which will become part of a filename
%
% For geometric searches, these are required:
%   Query.Filename       a string like 1s72, where the query motif is found
%   Query.NTList         a list of nucleotide numbers
%   Query.ChainList      a list of chain specifications, use if needed
%
% For non-geometric or mixed searches, these are optional:
%   Query.Edges          list of required basepairing or stacking interactions
%                        also specify allowed or disallowed pairs here
%   Query.Mask           a mask for which nucleotides to allow (defaults N)
%   Query.AngleWeight    weights to put on the angles (defaults 1)
%   Query.DistanceWeight weights to put on nucleotide distances (defaults 1)
%   Query.DiscCutoff     discrepancy cutoff D_0 (default 0.4)
%   Query.RelCutoff      relaxed cutoff D_1 (default Query.DiscCutoff)
%   Query.MaxDiff        maximum difference between sorted nucleotide indices
%   Query.MinDiff        minimum difference between sorted nucleotide indices
%
%   Query.Geometric      set to 0 to ignore geometry, only use screens. 
%                        Default is 1.
%   Query.ExcludeOverlap set to 1 to eliminate highly redundant motifs with
%     larger discrepancies; often, you get the same candidate with one or two
%     nucleotides different but much higher discrepancy.  Default is 1 when
%     Query.NumNT > 6.
 
function [Query] = xSpecifyQuery(QN);

if nargin > 0,
  Query.Number = QN;
else                        % change the following line to change the query!
  Query.Number =61; 
end

switch Query.Number

case 4
  Query.Description    = 'GRNA hairpin without sequential constraint';
  Query.Name           = 'GNRA4NonSeq';
  Query.Filename       = '1s72';
  Query.NTList         = {'804' '805' '808' '809'};
  Query.ChainList      = {'0' '0' '0' '0'}; 
  Query.Edges{1,4}     = 'cWW';
  Query.Edges{2,3}     = 'tSH';
  Query.DiscCutoff     = 1;      
  Query.ExcludeOverlap = 1;

case 41
  Query.Description    = 'GRNA hairpin with sequential constraint';
  Query.Name           = 'GNRA4Seq';
  Query.Filename       = '1s72';
  Query.NTList         = {'804' '805' '808' '809'};
  Query.ChainList      = {'0' '0' '0' '0'}; 
  Query.Edges{1,4}     = 'cWW';
  Query.Edges{2,3}     = 'tSH';
  Query.DiscCutoff     = 1;      
  Query.MaxDiff(1,4)   = 6;
  Query.ExcludeOverlap = 1;

case 5
  Query.Description    = 'Sarcin five nucleotide geometric';
  Query.Name           = 'Sarcin5Geo';
  Query.Filename       = '1s72';
  Query.NTList         = {'2694' '2701' '2693' '2702' '2692'};
  Query.ChainList      = {'0' '0' '0' '0' '0'};   % all in the 23S
  Query.DiscCutoff     = 0.5;

case 51
  Query.Description    = 'Sarcin five nucleotide symbolic';
  Query.Name           = 'Sarcin5Symb';
  Query.Edges{1,2}     = 'tHS';
  Query.Edges{3,5}     = 'cHS';
  Query.Edges{3,4}     = 'tWH';
  Query.MaxDiff(5,3)   = 2;
  Query.MaxDiff(3,1)   = 2;
  Query.MaxDiff(4,2)   = 2;

case 52
  Query.Description    = 'GRNA hairpin 5 nucleotide';
  Query.Name           = 'GNRA5';
  Query.Filename       = '1s72';
  Query.NTList         = {'804' '805' '807' '808' '809'};
  Query.ChainList      = {'0' '0' '0' '0' '0'}; 
  Query.Edges{1,5}     = 'cWW bif';
  Query.Edges{1,2}     = 's35';
  Query.Edges{3,4}     = 's35';
  Query.DiscCutoff     = 0.8;
  Query.MaxDiff(1,5)   = 6;
  Query.MaxDiff(2,4)   = 4;
  Query.ExcludeOverlap = 1;

case 6
  Query.Description    = 'Kink-turn central base mixed';
  Query.Name           = 'KinkTurnCentral';
  Query.Filename       = '1s72';
  Query.NTList         = {'80' '97' '81' '93' '94' '98'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0'};   % all in the 23S
  Query.Edges{1,2}     = 'tHS';
  Query.DiscCutoff     = 0.7;  
  Query.ExcludeOverlap = 1;

case 61
  Query.Description    = 'Kink-turn closing base pair mixed';
  Query.Name           = 'KinkTurnClosing';
  Query.Filename       = '1s72';
  Query.NTList         = {'80' '97' '81' '93' '100' '77'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0'};   % all in the 23S
  Query.Edges{1,2}     = 'tHS';
  Query.DiscCutoff     = 0.9;    
  Query.ExcludeOverlap = 1;

case 7
  Query.Description    = 'Sarcin seven nucleotide mixed';
  Query.Name           = 'Sarcin7Mixed';
  Query.Filename       = '1s72';
  Query.NTList         = {'2694' '2701' '2693' '2702' '2692' '2691' '2703'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0' '0'};   % all in the 23S
  Query.Edges{3,4}     = 'tWH';
  Query.DiscCutoff     = 0.5;       
  Query.ExcludeOverlap = 0;

case 9
  Query.Description    = 'Sarcin nine nucleotide mixed';
  Query.Name           = 'Sarcin9Mixed';
  Query.Filename       = '1s72';
  Query.NTList         = {'2694' '2701' '2693' '2702' '2692' '2691' '2703' '2690' '2704'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0' '0' '0' '0'};% all in the 23S
  Query.Edges{3,4}     = 'tWH';
  Query.DiscCutoff     = 0.5;
  Query.ExcludeOverlap = 1;


case 2
%  Query.Edges{1,2}  = 'ncWH ~cWH';                 % string
  Query.MaxDiff = [5];
  Query.MinDiff = [5];

case 21
  Query.ReqInter{1,2}  = [101];

case 411
  Query.Description    = 'GRNA hairpin without sequential constraint';
  Query.Filename       = '1s72';
  Query.NTList         = {'804' '809' '805' '808'};
  Query.ChainList      = {'0' '0' '0' '0'}; 
  Query.Mask           = 'NNNN';
  Query.Edges{1,2}     = 'ncWW';
  Query.Edges{3,4}     = 'ntSH';
  Query.DiscCutoff     = 1;         % guaranteed to find all candidates
                                    % with discrepancy below this number
  Query.ExcludeOverlap = 1;

case 42
  Query.Description    = 'Non-geometric search';
  Query.MaxDiff        = [4 Inf 4];
  Query.Edges{1,4}     = 'tHS';
  Query.Edges{2,3}     = 'tHS';

case 43
  Query.Description    = 'What stacks on a cWW?';
  Query.Edges{1,2}     = 'cWW';
  Query.Edges{3,4}     = 'Pair ~tSH ~tHS ~cWW';  % exclude categories we know
  Query.Edges{1,3}     = 'Stack';
  Query.Edges{2,4}     = 'Stack';
  Query.Diff{1,3}      = '=1';
  Query.Diff{2,4}      = '=1';

case 44
  Query.Edges{1,2} = 'tWH ~CA';
  Query.Edges{3,4} = 'tHS';
  Query.Edges{1,3} = 's35 s53 s55 s33';

case 501
  Query.Description    = 'Sarcin five nucleotide quick search';
  Query.Filename       = '1s72';
  Query.NTList         = {'2694' '2701' '2693' '2702' '2692'};
  Query.ChainList      = {'0' '0' '0' '0' '0'};   % all in the 23S
  Query.DiscCutoff     = 0.3;         % guaranteed to find all candidates
                                      % with discrepancy below this number

case 52
  Query.Description    = 'Sarcin five nucleotide non-geometric';
  Query.Mask           = 'NNNNN';
  Query.MaxDiff        = [2 2 Inf 2];
  Query.ReqInter{3,4}  = [10];
  Query.ReqInter{2,5}  = [-4];
  Query.ReqInter{2,1}  = [-9];

case 611
  Query.Description    = 'Kink-turn closing base pair search with specified basepairs';
  Query.Filename       = '1s72';
  Query.NTList         = {'80' '97' '81' '93' '100' '77'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0'};   % all in the 23S
  Query.ReqInter{1,2}  = [10];
  Query.ReqInter{3,4}  = [1];
  Query.ReqInter{5,6}  = [1];
  Query.DiscCutoff     = 0.9;         % guaranteed to find all candidates
                                      % with discrepancy below this number
  Query.ExcludeOverlap = 1;

case 622
  Query.Description    = 'Noncanonical pair between canonical';
  Query.Edges{1,6}     = 'cWW';
  Query.Edges{2,5}     = '~cWW';
  Query.Edges{3,4}     = 'cWW';
  Query.MaxDiff(1,2)   = 1;
  Query.MaxDiff(2,3)   = 1;
  Query.MaxDiff(4,5)   = 1;
  Query.MaxDiff(5,6)   = 1;

case 623
  Query.Description    = 'Non-geometric search';
  Query.MaxDiff(1,2)   = 1;
  Query.MaxDiff(2,3)   = 1;
  Query.MaxDiff(4,5)   = 1;
  Query.MaxDiff(5,6)   = 1;
  Query.Edges{1,6}     = 'cWW AU UA';
  Query.Edges{2,5}     = 'cWW ~CG ~GC';
  Query.Edges{3,4}     = 'cWW AU UA';



end

% Explanation of mask codes:
% A-A
% C-C
% G-G
% U-U
% A,C-M
% A,G-R
% A,U-W
% C,G-S
% C,U-Y
% G,U-K
% A,C,G-V
% A,C,U-H
% A,G,U-D
% C,G,U-B
% A,C,G,U-N

