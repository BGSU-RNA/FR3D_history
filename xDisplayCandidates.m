% the function xDisplayCandidates(File,Search) is used to display
% the Model and the best matches to the model

function Search = xDisplayCandidates(File,Search,Level)

if nargin < 3,
  MenuTitle = 'Display options';
  Level     = 0;
  QuitButton = 'Quit';
else
  MenuTitle = ['Level ' num2str(Level)];
  QuitButton = 'Quit level';
end

if isempty(Search.Candidates)
  disp('There are no candidates to display')
  return
end

N = Search.Query.NumNT;

warning off

% if there is no geometric model, use the first candidate to align to

Model = Search.Query;

if Model.Geometric == 0,
  f             = Search.Candidates(1,N+1);
  Model.Indices = double(Search.Candidates(1,1:N));
  Model.NT      = File(f).NT(Model.Indices);
  Model.Centers = cat(1,Model.NT(:).Center);
  Model.Filename = '';
  Model.WeightedCenter = mean(Model.Centers);
  Model.WeightedCenteredCenters = Model.Centers-ones(N,1)*mean(Model.Centers);
end

L = length(Search.Candidates(:,1));  % number of candidates found

if ~isfield(Search,'Marked'),
  Search.Marked = zeros(1,L);         % allow marking certain candidates
end

if ~isfield(Search,'Align'),
  Search.Align = zeros(1,L);         % allow marking certain candidates
end

% -------------------------------------------------------------------------

Display(1).n = 1;               % which candidate is in display window 1
Display(1).sugar = 1;           % display sugars or not
Display(1).neighborhood = 0;    % how large a neighborhood to show
Display(1).superimpose  = 0;    % superimpose the first candidate?
Display(1).supersugar   = 0;    % show sugars of first when superimposing?
Display(1).az           = -37.5;
Display(1).el           = 30;

Numplots = 1;
stop     = 0;
i        = 1;                              % current window
PlotMotif(File,Search,Model,Display,i);    % graph in display window i
rotate3d on

% ------------------------------- display menu -----------------------------

while stop == 0,                            

  k=menu(MenuTitle,'Next plot','Previous Plot','',...
         'Add plot','Larger Neighborhood', ...
         '','Toggle sugar','Toggle superimpose', ...
         'List Marked to screen','Mark/Unmark current','Reverse all marks', ...
         'Write marked to PDB','Display marked','Group marked', ...
         'Align marked',QuitButton);

  ii=gcf;                                 % get current active figure
  if (abs(ii) > 10) | (ii == 0),
    ii = i;
  end
  i = ii;


  [az,el]=view;                          % get current view (orientation)
  Display(i).az = az;
  Display(i).el = el;
  Display(i).x=XLim;                     % current x, y, z limits
  Display(i).y=YLim;
  Display(i).z=ZLim;

  switch k                               % k is the menu choice
    case 1                                      % next plot
      Display(i).n = min(L,Display(i).n+1);     % move to next candidate

    case 2                                      % Previous Plot
      Display(i).n = max(1,Display(i).n-1);

    case 4                                      % Add plot
      Numplots = Numplots + 1;                  % increase number of windows
      Display(Numplots) = Display(i);           % use current settings
      i = Numplots;
      figure(i);
      AlignPlots = 'AlignPlots';

    case 5                                      % toggle Neighborhood
      dn = [4 4 4 4 6 6 8 8 0];                 % neighborhood setting list
      Display(i).neighborhood = dn(1+Display(i).neighborhood);

    case 6                                      % 

    case 7                                      % toggle sugar
      if Display(i).superimpose == 0,
        Display(i).sugar = 1 - Display(i).sugar;
      elseif (Display(i).sugar == 0) & (Display(i).supersugar == 0),
        Display(i).sugar = 1;
      elseif (Display(i).sugar == 1) & (Display(i).supersugar == 0),
        Display(i).supersugar = 1;
      elseif (Display(i).sugar == 1) & (Display(i).supersugar == 1),
        Display(i).sugar = 0;
      elseif (Display(i).sugar == 0) & (Display(i).supersugar == 1),
        Display(i).supersugar = 0;
      end

    case 8                                      % toggle superimpose
      Display(i).superimpose = 1-Display(i).superimpose;

    case 9                                      % list marked on screen
      Search2 = Search;
      j = find(Search.Marked);
      Search2.Candidates  = Search.Candidates(j,:);
      Search2.Discrepancy = Search.Discrepancy(j);
      if length(j) > 0,
        fprintf('Marked candidates:\n');
        xListCandidates(File,Search2);
      end

    case 10                                      % mark current cand
      Search.Marked(Display(i).n) = 1-Search.Marked(Display(i).n);  % toggle mark
      Display(i).n = min(L,Display(i).n+1);             % move to next

    case 11                                    % reverse all marks
      Search.Marked = 1-Search.Marked;
      for j=1:Numplots
        PlotMotif(File,Search,Model,Display,j);
      end

    case 12                                     % write PDB of marked
      Search2 = Search;
      j = find(Search.Marked);
      Search2.Candidates  = Search.Candidates(j,:);
      Search2.Discrepancy = Search.Discrepancy(j);
      Search2.Marked      = Search.Marked(j);
      if length(j) > 0,
        xWriteCandidatePDB(File,Search2);
      end

    case 13                                     % display marked 
      Search2 = Search;
      j = find(Search.Marked);
      if length(j) > 0,
        Search2.Candidates  = Search.Candidates(j,:);
        Search2.Discrepancy = Search.Discrepancy(j);
        Search2.Marked      = Search.Marked(j);
        xDisplayCandidates(File,Search2,Level+1);
      end

    case 14                                     % group and display marked 
      Search2 = Search;
      j = find(Search.Marked);
      if length(j) > 0,
        Search2.Candidates  = Search.Candidates(j,:);
        Search2.Discrepancy = Search.Discrepancy(j);
        Search2.Marked      = Search.Marked(j);
        xGroupCandidates(File,Search2);
      end

    case 15                                     % align marked 
      Search2 = Search;
      j = find(Search.Marked);
      if length(j) > 0,
        Search2.Candidates  = Search.Candidates(j,:);
        Search2.Discrepancy = Search.Discrepancy(j);
        xAlignCandidates(File,Search2,1)
      end

    case 16                                     % quit Display
      if exist('fidOUT','var')
        fclose(fidOUT);
      end
      stop = 1;

    end  % switch statement for menu

  if any([1 2 4 5 6 7 8 9 10 13 14] == k),
    PlotMotif(File,Search,Model,Display,i);
    drawnow
  end

  if Numplots > 1,
      for j=1:Numplots,
        figure(j)
        sh(j) = subplot(1,1,1);
      end
      linkobj = linkprop(sh,...
{'cameraposition',...
 'cameraupvector',...
 'cameratarget',...
 'cameraviewangle'});

      set(gcf, 'UserData', linkobj);
  end

  figure(i)
  rotate3d on

end  % end while

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

function  PlotMotif(File,Search,Model,Display,i)

  N = Model.NumNT;

  figure(i)
  clf

  if (Display(i).superimpose == 1),
    Indices = Model.Indices;
    if (Model.NumNT > 2),
      R = eye(3);
      S = mean(cat(1,Model.NT.Center));
    else
      R = Model.NT(1).Rot;
      S = mean(cat(1,Model.NT.Center));
    end
    MVP.Rotation = R;
    MVP.Shift    = S;
    MVP.LineStyle = '-.';
    MVP.LineThickness = '1';
    MVP.Sugar    = Display(i).supersugar;
    MVP.ConnectSugar = 0;
    MVP.Grid     = 0;
    zDisplayNT(Model,1:N,MVP);
  end

  [s,t] = size(Search.Candidates);
  n       = Display(i).n;
  f       = Search.Candidates(n,N+1);
  Disc    = Search.Discrepancy(n);
  Indices = double(Search.Candidates(n,1:N));
  zShowInteractionTable(File(f),double(Indices),Disc);

  drawnow

  VP.Sugar    = Display(i).sugar;

  if Model.NumNT > 2,
    MC = Model.WeightedCenteredCenters;          % align to the model
    CandiCenters = cat(1,File(f).NT(Indices).Center);
    CC = CandiCenters - ones(N,1)*mean(CandiCenters);

    R = zBestRotation(MC, CC);
    S = mean(CandiCenters);
  else
    R = File(f).NT(Indices(1)).Rot;
    S = mean(cat(1,File(f).NT(Indices).Center));
  end

  VP.Rotation = R;
  VP.Shift    = S;
  VP.Grid     = 0;

  if Display(i).neighborhood > 0,
    a = zeros(1,File(f).NumNT);
    v = Display(i).neighborhood;
    for j=1:length(Indices),
      a = a + (File(f).Distance(Indices(j),:) < v) .* ...
              (File(f).Distance(Indices(j),:) > 0);
    end
    a(Indices) = zeros(1,length(Indices));  % take out ones in Indices
    B = find(a);
    Indices = [Indices B];
  end

  zDisplayNT(File(f),Indices,VP);
 
  xlabel(['Plot ',int2str(n),' of ',int2str(s),'   Disc: ',...
          num2str(Disc)]);
  if Search.Marked(n) == 1;
    yl = 'Marked';
  else
    yl = '';
  end
  ylabel(yl);

  axis equal
  axis vis3d
  view([Display(i).az Display(i).el]);
  drawnow
%  grid off

%    grid on
%    set(gcf,'Renderer','OpenGL')
%    set(gcf,'Renderer','zbuffer')

