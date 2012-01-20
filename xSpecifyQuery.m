% xSpecifyQuery returns the description of a model motif
% The variable Query has several fields:
%   Query.Filename       a string like rr0033_23S
%   Query.NTList         a list of nucleotide numbers
%   Query.Mask           a mask for which nucleotides to allow
%   Query.AngleWeight    weights to put on the angles
%   Query.DistanceWeight weights to put on nucleotide distances
%   Query.DiscCutoff     discrepancy cutoff D_0
%   Query.RelCutoff      relaxed cutoff D_1
%   Query.MaxDiff        maximum difference between sorted nucleotide
%     indices; use with Sequential = 1 to screen for sequential constraints
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
  Query.Number = 51; 
end

switch Query.Number

case 4
  Query.Description    = 'GRNA hairpin with sequential constraint';
  Query.Filename       = '1s72';
  Query.NTList         = {'804' '809' '805' '808'};
  Query.ChainList      = {'0' '0' '0' '0'}; 
  Query.Mask           = 'NNNN';
  Query.ReqInter{1,2}  = [1  30];
  Query.ReqInter{3,4}  = [10 -10 30];
  Query.DiscCutoff     = 1;         % guaranteed to find all candidates
                                    % with discrepancy below this number
  Query.MaxDiff        = [1 4 2 1];
  Query.ExcludeOverlap = 1;

case 41
  Query.Description    = 'GRNA hairpin without sequential constraint';
  Query.Filename       = '1s72';
  Query.NTList         = {'804' '809' '805' '808'};
  Query.ChainList      = {'0' '0' '0' '0'}; 
  Query.Mask           = 'NNNN';
  Query.ReqInter{1,2}  = [1  30];
  Query.ReqInter{3,4}  = [10 -10 30];
  Query.DiscCutoff     = 1;         % guaranteed to find all candidates
                                    % with discrepancy below this number
  Query.ExcludeOverlap = 1;

case 5
  Query.Description    = 'Sarcin five nucleotide';
  Query.Filename       = '1s72';
  Query.NTList         = {'2694' '2701' '2693' '2702' '2692'};
  Query.ChainList      = {'0' '0' '0' '0' '0'};   % all in the 23S
  Query.DiscCutoff     = 0.5;         % guaranteed to find all candidates
                                      % with discrepancy below this number

case 51
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

case 6
  Query.Description    = 'Kink-turn central base search';
  Query.Filename       = '1s72';
  Query.NTList         = {'80' '97' '81' '93' '94' '98'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0'};   % all in the 23S
  Query.ReqInter{1,2}  = [10];
  Query.DiscCutoff     = 0.7;         % guaranteed to find all candidates
                                      % with discrepancy below this number
  Query.ExcludeOverlap = 1;

case 61
  Query.Description    = 'Kink-turn closing base pair search';
  Query.Filename       = '1s72';
  Query.NTList         = {'80' '97' '81' '93' '100' '77'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0'};   % all in the 23S
  Query.ReqInter{1,2}  = [10];
  Query.DiscCutoff     = 0.9;         % guaranteed to find all candidates
                                      % with discrepancy below this number
  Query.ExcludeOverlap = 1;

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

case 7
  Query.Description    = 'Sarcin seven nucleotide';
  Query.Filename       = '1s72';
  Query.NTList         = {'2694' '2701' '2693' '2702' '2692' '2691' '2703'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0' '0'};   % all in the 23S
  Query.ReqInter{3,4}  = [-4];
  Query.DiscCutoff     = 0.5;         % guaranteed to find all candidates
                                      % with discrepancy below this number
  Query.ExcludeOverlap = 0;

case 9
  Query.Description    = 'Sarcin 9 nucleotide';
  Query.Filename       = '1s72';
  Query.NTList         = {'2694' '2701' '2693' '2702' '2692' '2691' '2703' '2690' '2704'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0' '0' '0' '0'};   % all in the 23S
  Query.ReqInter{3,4}  = [-4];
  Query.DiscCutoff     = 0.5;         % guaranteed to find all candidates
                                      % with discrepancy below this number
  Query.ExcludeOverlap = 1;


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

