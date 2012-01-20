
% This is a script for reading 3D alignment Excel spreadsheets from Jesse
% Stombaugh

% do this once:
[n,t,r] = xlsread('23S_Ec_Tt_Struct_alignment_12_4_06_JS.xls');

for i = 1:length(n(:,1)),
  A{1,i} = {num2str(n(i, 1)), num2str(n(i, 3))};        % basepairs from file 1
  A{2,i} = {num2str(n(i,14)), num2str(n(i,16))};        % basepairs from file 1
  A{3,i} = {num2str(n(i,26)), num2str(n(i,28))};        % basepairs from file 1
end

Filenames{1} = t{2,2};
Filenames{2} = t{2,15};
Filenames{3} = t{2,27};

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData(Filenames,0);   % load PDB data
else
  [File,SIndex] = zAddNTData(Filenames,0,File); % add PDB data if needed
end

L = length(Filenames);

NumPlots     = (L^2 - L)/2;

stop         = 0;
row          = 300;                          % row of alignment to start at
neighborhood = 1;
sugar        = 1;
fontsize     = 8;
labelbases   = 8;
Write        = 0;
inter        = 1;

% ------------------------------- display menu -----------------------------

while stop == 0,                            

  k=menu('3D alignment','Next element','Previous element', ... % 1,2
         'Change Neighborhood', ...                       % 3
         'Toggle sugar','Toggle display', ...             % 4,5
         'Write to PDB', ...                              % 6
         'Show Interaction Matrices', 'Quit');            % 7,8

  switch k                               % k is the menu choice
    case 1                                      % next row of alignment
      row = row + 1;

    case 2                                      % previous row of alignment
      row = row - 1;

    case 3                                      % enlarge neighborhood
      dn = [2 3 4 1];                           % neighborhood setting list
      neighborhood = dn(neighborhood);

    case 4                                      % toggle sugar
      sugar = 1 - sugar;

    case 5                                      % toggle numbers
      if labelbases == 0,
        labelbases = fontsize;
      elseif labelbases > 0,
        labelbases = 0;
      end

    case 6                                      % write to PDB
      Write = 1;

    case 7                                      % show interaction matrices
      inter = 1 - inter;

    case 8                                      % quit
      stop = 1;

  end  % switch statement for menu

  VP.Sugar = sugar;
  VP.LabelBases = labelbases;

  VP.LineStyle     = '-';
  VP.LineThickness = '1';
  VP.ConnectSugar  = 1;
  VP.Grid          = 0;

  if any([1 2 3 4 5 6] == k),

    % ------------ Make lists of nucleotides shared by all molecules

    startrow = max(1,row-neighborhood);
    endrow   = min(length(A(1,:)), row + neighborhood);
    for k = 1:L,
      NTList{k} = [];
    end

    for r = startrow:endrow,
      OK = 1;
      for k = 1:L,
        if strcmp(A{k,r}{1},'NaN') || strcmp(A{k,r}{2},'NaN'),
          OK = 0;
        end
      end
      if OK > 0,
        for k = 1:L,
          a = length(NTList{k});
          NTList{k}{a+1} = A{k,r}{1};
          NTList{k}{a+2} = A{k,r}{2};
        end
      end
    end

    if inter == 0,

      % ------------ Display nucleotide lists

      for k=1:L,
        fprintf('%s nucleotides', Filenames{k});
        for m = 1:length(NTList{k}),
          fprintf(' %s', NTList{k}{m});
        end
        fprintf('\n');
      end

    else

      % ------------ Show interaction matrices

      for k=1:L,
        zShowInteractionTable(File(k),NTList{k});
      end

    end

    % ------------ Superimpose and display nucleotides
    if length(NTList{1}) > 2,
      f = 1;
      for i= 1:L,
        for j = (i+1):L,    
          figure(f);
          [az,el]=view;                  % get current view (orientation)
          clf
         zSuperimposeNucleotides(File(i),NTList{i},File(j),NTList{j},VP,Write);
          view(az,el);

          f = f + 1;  
        end
      end
    else
      fprintf('Increase size of neighborhood\n');
    end



  end

  rotate3d on
  drawnow

  Write = 0;

end  % end while



