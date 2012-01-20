% zFindExemplars finds the best representative for each category of pairs.

% Pair codes:  1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU

pcodes = [6 7 13 14 15];
pcodes = [1 5 6 7 9 11 13 14 15 16];    % pair codes to work on

LMax = 500;                % maximum number of pairs to consider in each class

load('PairExemplars','Exemplar');

CL = zClassLimits;

Pairs{1} = 'AA';
Pairs{5} = 'AC';
Pairs{6} = 'CC';
Pairs{7} = 'GC';
Pairs{9} = 'AG';
Pairs{11} = 'GG';
Pairs{13} = 'AU';
Pairs{14} = 'CU';
Pairs{15} = 'GU';
Pairs{16} = 'UU';

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData('Nonredundant_3_list',0);   % load PDB data
else
  [File,SIndex] = zAddNTData('Nonredundant_3_list',0,File); % add PDB data if needed
end                       
                          % must load full .mat files!
File = File(SIndex);

if ~exist('ExtraFile'),                           % if no molecule data is loaded,
  [ExtraFile,SIndex] = zAddNTData('ExtraInstance_list',0);   % load PDB data
else
  [ExtraFile,SIndex] = zAddNTData('ExtraInstance_list',0,ExtraFile); % add PDB data if needed
end                       

ExtraFile = ExtraFile(SIndex);

if ~exist('ModelFile'),                           % if no molecule data is loaded,
  [ModelFile,SIndex] = zAddNTData('Model_list',0);   % load PDB data
else
  [ModelFile,SIndex] = zAddNTData('Model_list',0,ModelFile); % add PDB data if needed
end                       

ModelFile = ModelFile(SIndex);

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

for j = 1:length(pcodes),            % run through all pair codes specified
 pc = pcodes(j);
 CLE = [CL(:,1,pc); [21 22 23]'];    % include stacking, identified differently
 CLE = CLE(find(CLE));               % leave out empty entries

 fprintf('\n');

 for row = 1:length(CLE),

  % specify criteria for selection of pairs ----------------------------------

  Param.Paircode = pc;
  Param.Category = CLE(row);        
  Param.Decimal  = 1;        % 1 - use 1.0 only; 0 - round to 1; may not work!
  Param.Group    = 1;        % computer classification matches
  Param.Sequential= 0;

%  fprintf('Paircode %2d Class %5.1f ', pc, CLE(row));

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

  if length(SP) > 0,               % instances of this category are found

    L = min(LMax,length(SP));      % Limit the number of pairs to consider
    PD = zeros(L,L);
    for k = 1:L,                   % Very slow nested loop
      for m = (k+1):L,
  
        f1 = SP(k).Filenum;
        f2 = SP(m).Filenum;
        Model = [File(f1).Pair(SP(k).PairIndex).Base1Index ...
                 File(f1).Pair(SP(k).PairIndex).Base2Index];
        Cand  = [File(f2).Pair(SP(m).PairIndex).Base1Index ...
                 File(f2).Pair(SP(m).PairIndex).Base2Index];
        PD(k,m) = xDiscrepancy(File(f1),Model,File(f2),Cand);
      end
    end

    PD = sqrt(PD)/2;            % finish discrepancy calculation

    fprintf('Found %4d instances of %2s %4s class %5.1f', length(SP), Pairs{pc}, zEdgeText(CLE(row),1,pc), CLE(row));

    bigm = max(max(PD));
    if bigm > 1,
      fprintf('%6.2f maximum pair discrepancy\n',bigm);
    else
      fprintf('\n');
    end

    PD = PD + PD';
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

  else                                  % no instances of this category found

    SP = zSelectPairs(ExtraFile,Param); % check extra basepairs

    if length(SP) > 0,

      fprintf('Found %4d instances of %2s %4s class %5.1f outside the non-redundant set\n', length(SP), Pairs{pc}, zEdgeText(CLE(row),1, pc), CLE(row));

      f = SP(1).Filenum;
      p = SP(1).PairIndex;
      E = ExtraFile(f).Pair(p);
      Exemplar(row,pc).Filename   = ExtraFile(f).Filename;
      Exemplar(row,pc).Class      = CLE(row);
      Exemplar(row,pc).NT1        = ExtraFile(f).NT(E.Base1Index);
      Exemplar(row,pc).NT2        = ExtraFile(f).NT(E.Base2Index);
      Exemplar(row,pc).Pair       = E;
      Exemplar(row,pc).Count      = 1;

    else

      SP = zSelectPairs(ModelFile,Param); % check modeled basepairs

      if length(SP) > 0,

        fprintf('Using a model for       %2s %4s class %5.1f\n', Pairs{pc}, zEdgeText(CLE(row),1), CLE(row));

        f = SP(1).Filenum;
        p = SP(1).PairIndex;
        E = ModelFile(f).Pair(p);
        Exemplar(row,pc).Filename   = ModelFile(f).Filename;
        Exemplar(row,pc).Class      = CLE(row);
        Exemplar(row,pc).NT1        = ModelFile(f).NT(E.Base1Index);
        Exemplar(row,pc).NT2        = ModelFile(f).NT(E.Base2Index);
        Exemplar(row,pc).Pair       = E;
        Exemplar(row,pc).Count      = 0;
      else
        fprintf('No instances and no model for %2s %4s %6.1f\n', Pairs{pc}, zEdgeText(CLE(row),1), CLE(row));
      end
    end
  end


  % add information to speed up the discrepancy calculation later

  if ~isempty(Exemplar(row,pc).NT1),
    Exemplar(row,pc).R            = Exemplar(row,pc).NT2.Rot' * Exemplar(row,pc).NT1.Rot;
    Exemplar(row,pc).T1           = (Exemplar(row,pc).NT2.Center - Exemplar(row,pc).NT1.Center) * Exemplar(row,pc).NT1.Rot;
    Exemplar(row,pc).T2           = (Exemplar(row,pc).NT1.Center - Exemplar(row,pc).NT2.Center) * Exemplar(row,pc).NT2.Rot;
    Exemplar(row,pc).AngleWeight  = [1 1];
    Exemplar(row,pc).LDiscCutoff  = Inf;
  end 

  save(['FR3DSource' filesep 'PairExemplars'],'Exemplar'); % Matlab version 7 only
  save PairExemplars_Version_6.mat Exemplar -V6 % for compatibility with older versions

 end
end
