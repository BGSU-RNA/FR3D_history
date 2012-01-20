% zFindExamplars finds/updates the best representative for each category of
% pairs.  The user specifies the paircodes to update.

LMax = 500;                % maximum number of pairs to consider in each class

%load('PairExemplars','Exemplar');

CL = zClassLimits;

                    % 1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
pcodes = [1 5 6];
                    % 1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU
pcodes = [1 5 6 7];

pcodes = [9 11 13 14];

pcodes = [6 7 13 14 15];

pcodes = [7 16];

pcodes = [1];
pcodes = [1 5 6 7 9 11 13 14 15 16];

% load data ----------------------------------------------------------------

if ~exist('File'),
  Filenames = zFileNameList;
  for f=1:length(Filenames),
    File(f) = zGetNTData(Filenames{f},0);
  end
end

% specify parameters for viewing -------------------------------------------

  ViewParam.Mode      = 1; 
  ViewParam.Color     = 1;
  ViewParam.Normal    = 1;
  ViewParam.ColorAxis = [-12 30];
  ViewParam.SortKeys  = [];
  ViewParam.Nearby    = 0;
  ViewParam.Sugar     = 0;
  ViewParam.Hydrogen  = 1;
  ViewParam.Sort      = 0;
  ViewParam.az        = 51;
  ViewParam.el        = 14;

% loop through paircodes and computer classifications ----------------------

for j = 1:length(pcodes),
 pc = pcodes(j);
 CLE = [CL(:,1,pc); [21 22 23]'];    % include stacking, identified differently
 CLE = CLE(find(CLE));               % leave out empty entries
 for row = 1:length(CLE),

  % specify criteria for selection of pairs ----------------------------------

  Param.Paircode = pc;
  Param.Category = CLE(row);        
  Param.Decimal  = 1;        % 1 - use 1.0 only; 0 - round to 1; may not work!
  Param.Group    = 1;        % computer classification matches
  Param.Sequential= 0;
  Param.Inbox    = [-5 5 -5 5 -5 5 -1.1 1.1 -95 275];  % if category = 50,
                           % only pairs in this box will be selected

  fprintf('Paircode %2d Class %5.1f ', pc, CLE(row));

  % select pairs using selection criteria -----------------------------------

  SP = zSelectPairs(File,Param);

  % Here is a brief summary of the format of SP:
  %   SP is an array of structured variables, one for each Selected Pair.
  %   SP(i).Filenum        The number of the file from which this pair comes
  %   SP(i).B1Index        The index of the first base in the pair
  %   SP(i).B2Index        The index of the second base in the pair
  %   SP(i).PairIndex      The index of the pair of bases
  %   SP(i).CI             The hand index of the pair (may be zero)
  %   SP(i).HandClass      The expert classification; 0 if none exists
  %   SP(i).MinDist        The minimum distance between atoms in these bases
  %   SP(i).C1pC1p         The C1' - C1' distance for this pair

  if length(SP) > 0,
    L = min(LMax,length(SP));      % Limit the number of pairs to consider
    PD = zeros(L,L);
    xd = zeros(L,L);
    for k = 1:L,                   % Very slow nested loop
      for m = (k+1):L,
  
        Pair1 = File(SP(k).Filenum).Pair(SP(k).PairIndex);
        Pair2 = File(SP(m).Filenum).Pair(SP(m).PairIndex);
%        PD(k,m) = zPairDiscrepancy(Pair1,Pair2);  % does this match FR3D?
        f1 = SP(k).Filenum;
        f2 = SP(m).Filenum;
        Model = [File(f1).Pair(SP(k).PairIndex).Base1Index ...
                 File(f1).Pair(SP(k).PairIndex).Base2Index];
        Cand  = [File(f2).Pair(SP(m).PairIndex).Base1Index ...
                 File(f2).Pair(SP(m).PairIndex).Base2Index];
        PD(k,m) = xDiscrepancy(File(f1),Model,File(f2),Cand);
      end
    end

%a = nonzeros(PD);
%b = nonzeros(xd);
%[length(a) length(b)]
%clf
%plot(a,b,'.');
%axis([0 max(a) 0 max(b)])
%pause

%    zListPairs(File,SP(1:L),2);
    PD = PD + PD';
%    PD(1:20,1:10)
%    PD

    rs = sum(PD);

    [y,i] = sort(rs);

    for k = 1:min(10,L),
%      zDisplayPair(File(SP(i(k)).Filenum),SP(i(k)),ViewParam);
%      pause
%      [ViewParam.az,ViewParam.el] = view;
    end

    f = SP(i(1)).Filenum;
    p = SP(i(1)).PairIndex;
    E = File(f).Pair(p);
    Exemplar(row,pc).Filename   = File(f).Filename;
    Exemplar(row,pc).Class      = CLE(row);
    Exemplar(row,pc).NT1        = File(f).NT(E.Base1Index);
    Exemplar(row,pc).NT2        = File(f).NT(E.Base2Index);
    Exemplar(row,pc).Pair       = E;
    Exemplar(row,pc).Count      = length(SP);

    % The following is now redundant; can remove these lines later

    Exemplar(row,pc).Displ      = E.Displ;
    Exemplar(row,pc).Rot        = E.Rot;
    Exemplar(row,pc).Base1Index = E.Base1Index;
    Exemplar(row,pc).Base2Index = E.Base2Index;
    Exemplar(row,pc).Base1      = File(f).NT(E.Base1Index).Number;
    Exemplar(row,pc).Base2      = File(f).NT(E.Base2Index).Number;
  
  end

 save('PairExemplarsNew','Exemplar');

 end
end
