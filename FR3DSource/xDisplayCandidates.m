% xDisplayCandidates(File,Search,Level) displays candidates
% graphically and allows various types of analysis

% It may be run directly from Matlab using the command:
%    [Search,File] = xDisplayCandidates([],Search);
% and after this,
%    [Search,File] = xDisplayCandidates(File,Search);

function [Search, File] = xDisplayCandidates(FullFile,Search,Level,UsingFull,Order)

if strcmp(class(Search),'double'),
  S = Search;
  clear Search
  Search.Candidates = S;
end

if isempty(Search.Candidates)
  fprintf('There are no candidates to display\n');
  File = FullFile;
  return
end

[L,N] = size(Search.Candidates);
N = N - 1;                                   % number of nucleotides

Limit = min(L,300);                          % for mutual discrep matrix
p = 1:L;                                     % default permutation for display
q(p) = 1:L;                                  % inverse permutation
                               % p : display position -> real candidate number
                               % q : real candidate number -> display position

if ~isfield(Search,'File'),
  UsingFull = 1;                             % use what was passed in
  File = FullFile;
  Search.Query.Geometric = 0;
  Search.Query.NumNT = N;
  Search.Discrepancy = 1:L;
  FIndex = 1:length(FullFile);
else
  UsingFull = 0;
  File = Search.File;                        % use what was saved w/ Search
  FIndex = 1:length(Search.File);
end

Query = Search.Query;

fontsize = 10;                               % for nucleotide numbers

if nargin < 3,
  MenuTitle  = 'Display options';
  Level      = 0;
  QuitButton = 'Quit display';
else
  MenuTitle  = ['Subset depth ' num2str(Level)];
  QuitButton = 'Return to larger set';
end

if nargin < 5,
  if Query.Geometric == 0,
    Order = 2;
  else
    Order = 1;
  end
end

OrderText = {'by discrepancy from query', 'by file, then sum of nucleotide numbers', 'by similarity', 'by centrality', 'by pair criteria'};

warning off

if ~isfield(Search,'Marked'),
  Search.Marked = zeros(1,L);         % allow marking certain candidates
end

if ~isfield(Search,'Disc'),
  Search.Disc         = sparse(zeros(L,L));
  Search.DiscComputed = sparse(zeros(1,L));
  if Query.Geometric > 0,
    Search.Disc(:,1) = Search.Discrepancy;
    Search.Disc(1,1) = 0;
    Search.DiscComputed(1,1) = 1;
  end
end

NeighMax = 4;

% ----------------------------- find maximum gap between candidate nucleotides

[y,r] = sort(Search.Candidates(1,1:N)); % put nucleotides in increasing order

if isfield(Query,'MaxDiffMat'),
  MaxDiff = diag(Query.MaxDiffMat(r,r),1);
else
  MaxDiff = Inf*ones(1,N-1);
end

maxinsert = zeros(1,N-1);
for c = 1:L,
  maxinsert = max(maxinsert,abs(diff(double(Search.Candidates(c,r))))-1);
end

Display(1).p         = r;
Display(1).MaxDiff   = MaxDiff;
Display(1).MaxInsert = maxinsert;

% ------------------------------------------- Display distance matrix

Search = xMutualDiscrepancy(File,Search,Limit); % calculate some discrepancies

for i=1:L,
  f = Search.Candidates(i,N+1);          % file number
  b = '';
  for j = 1:min(4,N),
    b = [b File(f).NT(Search.Candidates(i,j)).Base];
  end
  n = File(f).NT(Search.Candidates(i,1)).Number;
  n = sprintf('%4s',n);
  if Search.Query.Geometric > 0,
      if isfield(Search,'AvgDisc'),
        d = sprintf('%6.4f',Search.AvgDisc(i));
      else
        d = sprintf('%6.4f',Search.Discrepancy(i));
      end
    else
      d = sprintf('%5d',Search.Discrepancy(i)); % orig candidate number
    end
  Search.Lab{i} = [b n ' ' File(f).Filename];
end

% --------- if there is no geometric model, align to the central candidate

if Query.Geometric == 0,
  [z,j] = sort(sum(Search.Disc));           % sort by average discrepancy
  f              = Search.Candidates(j(1),N+1);
  Query.Indices  = double(Search.Candidates(j(1),1:N));
  Query.NT       = File(f).NT(Query.Indices);
  Query.LocWeight= ones(1,Query.NumNT);
  Query          = xPrecomputeForDiscrepancy(Query);
  Query.Filename = 'Central candidate';
end

% ------------------------------------------- Parameters to display candidates

Display(1).n            = 1;     % which candidate is in display window 1
Display(1).sugar        = 1;     % display sugars or not
Display(1).neighborhood = 0;     % how large a neighborhood to show
Display(1).superimpose  = 0;     % superimpose the first candidate?
Display(1).supersugar   = 0;     % show sugars of first when superimposing?
Display(1).labelbases   = 10;    % show nucleotide numbers
Display(1).az           = -37.5; % standard view
Display(1).el           = 30;

stop     = 0;                              % stop the menu?
i        = 1;                              % current window
nn       = 1;                              % current candidate

PlotMotif(File,Search,Query,Display,i);    % graph in display window i
rotate3d on
DisplayTable(File,Search,Query,Display,i)
drawnow

% ------------------------------- display menu -----------------------------

figure(99)
axis([1 Limit+1 1 Limit+1]);

while stop == 0,                            
  % ---------------------------------------- Display table of discrepancies
  figure(99)
  ax = axis;
  clf
  pp = p(1:Limit);
  zGraphDistanceMatrix(Search.Disc(pp,pp),Search.Lab(pp));
  hold on
  co = {'w*','wo','wd','ws','wv','w<','w>','w^','w+','wx'};
  co = [co co co co co co co co];
  for j = 1:length(Display),
    plot(q(Display(j).n)+0.5,q(Display(j).n)+0.5,co{j});
  end
  m = q(find(Search.Marked));
  plot(m+0.5,m+0.5,'w.');
%  axis(ax);
  title(['Discrepancies between candidates, ordered by ' OrderText{Order}]);
  colormap('default');
  map = colormap;
  map = map((end-8):-1:8,:);
  colormap(map);
  caxis([0 0.8]);
  colorbar('location','eastoutside');



  if (Display(1).neighborhood == NeighMax),
    Neighborhood = 'No Neighborhood';
  else
    Neighborhood = 'Larger Neighborhood';
  end

  Buttons = {'Next candidate','Previous Candidate', ... % 1,2
         'Add plot', Neighborhood, ...                % 3,4
         'Toggle sugar','Toggle display', ...                 % 5,6
         'Mark/Unmark current','Reverse all marks', ...       % 7,8
         'Display marked only', ...                           % 9
         'List to screen','Write to PDB', ...                 % 10,11
         'Sort by centrality', 'Order by Similarity', ...     % 12,13
         'Show Alignment', ...                                % 14
         'Show Scatterplot', 'Navigate with Fig 99', ...      % 15, 16
         QuitButton};                                         % 17

  k=menu(MenuTitle,Buttons);

  ii=gcf;                                 % get current active figure
  if (abs(ii) > length(Display)) | (ii == 0), % other window active?
    ii = i;
  end
  i = ii;                                 % record and save active figure
  i = min(i,length(Display));

  figure(i)
  [az,el]=view;                          % get current view (orientation)
  Display(1).az = az;
  Display(1).el = el;
  Display(1).x=XLim;                     % current x, y, z limits
  Display(1).y=YLim;
  Display(1).z=ZLim;

  switch k                               % k is the menu choice
    case 1                                      % next plot
      n = Display(i).n;                         % actual candidate displayed
      if q(n) + 1 > L,                          % q(n) is display order
        Limit = min(Limit*2,L);                 % increase limit
        Display(i).n = p(1);

Limit
disp('Calculating more discrepancies');
Search = xMutualDiscrepancy(File,Search,Limit); % calculate some discrepancies

p = 1:L;                                     % default permutation for display
q(p) = 1:L;                                  % inverse permutation

      else
        Display(i).n = p(q(n) + 1);
      end

    case 2                                      % Previous Plot
      n = Display(i).n;                         % actual candidate displayed
      if q(n) - 1 < 1,                          % q(n) is display order
        Display(i).n = p(L);
      else
        Display(i).n = p(q(n) - 1);
      end

    case 3                                      % Add plot
      Display(end+1) = Display(i);              % use current settings
      i = length(Display);                      % current figure number
      figure(i);

    case 4                                      % toggle Neighborhood
      Display(1).neighborhood = Display(1).neighborhood + 1;
      if Display(1).neighborhood > NeighMax,
        Display(1).neighborhood = 0; 
      end 

    case 5                                      % toggle sugar
      if Display(1).superimpose == 0,
        Display(1).sugar = 1 - Display(1).sugar;
      elseif (Display(1).sugar == 0) & (Display(1).supersugar == 0),
        Display(1).sugar = 1;
      elseif (Display(1).sugar == 1) & (Display(1).supersugar == 0),
        Display(1).supersugar = 1;
      elseif (Display(1).sugar == 1) & (Display(1).supersugar == 1),
        Display(1).sugar = 0;
      elseif (Display(1).sugar == 0) & (Display(1).supersugar == 1),
        Display(1).supersugar = 0;
      end

    case 6                                      % toggle superimpose/numbers
      if Display(1).superimpose == 0 & Display(1).labelbases == 0,
        Display(1).superimpose = 1;
        fprintf('Superimposing candidate in darker colors');
      elseif Display(1).superimpose == 1 & Display(1).labelbases == 0,
        Display(1).labelbases = fontsize;
      elseif Display(1).superimpose == 1 & Display(1).labelbases > 0,
        Display(1).superimpose = 0;
      elseif Display(1).superimpose == 0 & Display(1).labelbases > 0,
        Display(1).labelbases = 0;
      end

    case 7                                      % mark/unmark current cand
      Search.Marked(Display(i).n) = 1-Search.Marked(Display(i).n); %toggle
      n = Display(i).n;                         % actual candidate displayed
      if q(n) + 1 > L,                          % q(n) is display order
        Display(i).n = p(1);
      else
        Display(i).n = p(q(n) + 1);
      end

    case 8                                      % reverse all marks
      Search.Marked = 1-Search.Marked;

    case 9                                      % display marked only
      j = find(Search.Marked);
      if length(j) > 0,
        [y,m] = sort(q(j));
        j = j(m);                               % put j in display order
        Search2 = SearchSubset(Search,j);
        xDisplayCandidates(File(FIndex),Search2,Level+1);
        figure(99)
        axis(ax);
      end

    case 10                                      % list on screen
      j  = find(Search.Marked);
      jj = find(Search.Marked == 0);
      if (length(j) > 0) && (length(jj) > 0),
        [y,m] = sort(q(j));
        j = j(m);                                % put j in display order
        Search2 = SearchSubset(Search,j);
        fprintf('Marked candidates listed first\n');
        xListCandidates(Search2,Inf);

        [y,m] = sort(q(jj));
        jj = jj(m);                             % put jj in display order
        Search2 = SearchSubset(Search,jj);
        fprintf('Unmarked candidates listed second\n');
        xListCandidates(Search2,Inf);
      else
        Search2 = SearchSubset(Search,p);
        xListCandidates(Search2,Inf);
      end

    case 11                                     % write PDB of all
      Search2 = SearchSubset(Search,p);
      xWriteCandidatePDB(Search2);

    case 12                                     % sort by centrality
%     [z,j] = sort(max(Search.Disc(1:Limit,1:Limit)));% sort by max discrepancy
      [z,j] = sort(sum(Search.Disc));           % sort by average discrepancy
      S.AvgDisc  = z / (Limit - 1);             % average discrep among these
      p(1:Limit) = j;
      p((Limit+1):L) = (Limit+1):L;      
      q(p) = 1:L;
%      Search = xSortByCentrality(File(FIndex),Search,Level,UsingFull);
      Order = 4;

    case 13                                     % group candidates
%      Search = xGroupCandidates(File(FIndex),Search,Level,UsingFull);
%      Search = xOrderCandidates(File(FIndex),Search,Level,UsingFull);

      p(1:Limit) = zOrderbySimilarity(Search.Disc(1:Limit,1:Limit));
      p((Limit+1):L) = (Limit+1):L;
      q(p) = 1:L;
      Order = 3;

    case 14                                     % align
      xAlignCandidates(File(FIndex),Search,1);
%      Text = xFASTACandidates(File(FIndex),Search,1);
%      for t = 1:length(Text),
%        fprintf('%s\n', Text{t});
%      end

    case 15
      ViewParam.Color  = 6;
      ViewParam.FigNum = length(Display)+1;
      ViewParam.Normal = 0;
      ViewParam.ClassLimits = 1;
      ppp = xScatterPairs(Search,1,2,ViewParam);

    case 16
      figure(99)
      pt = get(gca,'CurrentPoint')

      if abs(pt(1,1)-pt(1,2)) > Limit/20,           % clicked off the diagonal
        Search.Marked = 0 * Search.Marked;      % unmark all candidates
        a = sort(pt(1,[1 2]));
        j = p(max(1,floor(a(1))):min(L,floor(a(2))));      % 
        Search.Marked(j) = ones(1,length(j));   % select these candidates
      else                                      % clicked near the diagonal
        newn = max(min(floor(pt(1,1)),L),1);
        Display(i).n = p(newn);
      end

    case 17                                     % quit Display
      if exist('fidOUT','var')
        fclose(fidOUT);
      end
      stop = 1;

  end  % switch statement for menu


  if any([1 2 3 7 16] == k),
      PlotMotif(File(FIndex),Search,Query,Display,i);
  end

  if (k == 4) && (UsingFull == 0),
    fprintf('Loading structure files\n');
    fprintf('If some are not available, Larger Neighborhood will crash\n');
    [File,FIndex] = zAddNTData(Search.CandidateFilenames,2,FullFile);
    for f = 1:length(File),
      if ~isfield(File,'Distance'),
        File(f).Distance = [];
      end
      if isempty(File(f).Distance) && ~isempty(File(f).NumNT),
       if (File(f).NumNT > 0),
        c = cat(1,File(f).NT.Center); % nucleotide centers
        File(f).Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
       end
      end
    end
    FullFile = [];
    UsingFull = 1;
  end

  if (Display(i).n ~= nn) || (k == 4),
    DisplayTable(File(FIndex),Search,Query,Display,i)
    nn = Display(i).n;
  end

  if any([4 5 6 8] == k),
    for j=1:length(Display)
      PlotMotif(File(FIndex),Search,Query,Display,j);
    end
  end

  if length(Display) > 1,
      for j=1:length(Display),
        figure(j)
        sh(j) = subplot(1,1,1);
        rotate3d on
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
  drawnow

end  % end while

if UsingFull == 0,
  File = FullFile;
  FullFile = [];
end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

function  PlotMotif(File,Search,Query,Display,i)

  N = Query.NumNT;

  figure(i)
  clf

  if (Display(1).superimpose == 1),
    Indices = Query.Indices;
    if (Query.NumNT > 2),
      R = eye(3);
      S = mean(cat(1,Query.NT.Center));
    else
      R = Query.NT(1).Rot;
      S = mean(cat(1,Query.NT.Center));
    end
    MVP.Rotation      = R;
    MVP.Shift         = S;
    MVP.LineStyle     = '-.';
    MVP.LineThickness = '1';
    MVP.Sugar         = Display(1).supersugar;
    MVP.ConnectSugar  = 0;
    MVP.Grid          = 0;
    MVP.LabelBases    = Display(1).labelbases;
    zDisplayNT(Query,1:N,MVP);
  end

  [s,t] = size(Search.Candidates);
  n       = Display(i).n;
  f       = Search.Candidates(n,N+1);
  Indices = double(Search.Candidates(n,1:N));

  VP.Sugar    = Display(1).sugar;
  VP.LabelBases = Display(1).labelbases;

  if Query.NumNT > 2,
    MC = Query.WeightedCenteredCenters;          % align to the model
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

  if Display(1).neighborhood > 0,
    v = Display(1).neighborhood;
    Indices = xNeighborhood(File(f),Indices,v,Display(1));
  end

  if exist('amal.txt','file') > 0,
     VP.AtOrigin = 1;
  end

  zDisplayNT(File(f),Indices,VP);

  if isfield(Search,'AvgDisc'),   
    xlabel(['Plot ',int2str(n),' of ',int2str(s),'   Average discrepancy from others ', num2str(Search.AvgDisc(n))]);
  elseif Query.Geometric > 0,
    xlabel(['Plot ',int2str(n),' of ',int2str(s),'   Discrepancy ',...
          num2str(Search.Discrepancy(n))]);
  else
    xlabel(['Plot ',int2str(n),' of ',int2str(s)]);
  end

  if Search.Marked(n) == 1;
    yl = 'Marked';
  else
    yl = '';
  end
  ylabel(yl);

  axis equal
  axis vis3d
  view([Display(1).az Display(1).el]);
  drawnow

  % Commands for Amal's study of base triples:

  if exist('amal.txt','file') > 0,

    N1 = File(f).NT(Indices(1));
    N2 = File(f).NT(Indices(2));
    N3 = File(f).NT(Indices(3));

    ytext = 'Interactions ';
    ytext = [ytext ' ' zEdgeText(File(f).Edge(Indices(1),Indices(2))) ' '];
    ytext = [ytext ' ' zEdgeText(File(f).Edge(Indices(2),Indices(3))) ' '];
    ytext = [ytext ' ' zEdgeText(File(f).Edge(Indices(1),Indices(3))) ' '];

    xlabel(ytext);
    view(2)

    if isfield(Search,'AvgDisc'),   
      ylabel(['Plot ',int2str(n),' of ',int2str(s),'   Average discrepancy from others ', num2str(Search.AvgDisc(n))]);
    elseif Query.Geometric > 0,
      ylabel(['Plot ',int2str(n),' of ',int2str(s),'   Discrepancy ',...
          num2str(Search.Discrepancy(n))]);
    else
      ylabel(['Plot ',int2str(n),' of ',int2str(s)]);
    end


  end
  % end of commands for Amal


%    set(gcf,'Renderer','OpenGL')
%    set(gcf,'Renderer','zbuffer')

% ------------------------------------------------- Display table

function  DisplayTable(File,Search,Query,Display,i)

    N       = Query.NumNT;
    n       = Display(i).n;
    f       = Search.Candidates(n,N+1);
    Indices = double(Search.Candidates(n,1:N));
    
    if isfield(Search,'AvgDisc'),
      fprintf('Average discrepancy from others %6.4f', Search.AvgDisc(n));
    elseif Query.Geometric > 0,
      fprintf('Discrepancy %6.4f', Search.Discrepancy(n));
    else
      fprintf('Candidate #%d', Search.Discrepancy(n));  % integer is cand num
    end

    if Display(1).neighborhood > 0,
      v = Display(1).neighborhood;
      Indices = xNeighborhood(File(f),Indices,v,Display(1));
    end

    zShowInteractionTable(File(f),double(Indices));

    if isfield(File(f),'BasePhosphate'),
      zBasePhosphateTable(File(f),double(Indices));
    end
    drawnow

% -------------------------------------------------- Select subset of search

function [Search2] = SearchSubset(Search,j)

  Search2             = Search;
  Search2.Candidates  = Search.Candidates(j,:);
  Search2.Discrepancy = Search.Discrepancy(j);
  Search2.Marked      = Search.Marked(j);
  Search2.Disc        = Search.Disc(j,j);
  Search2.DiscComputed= Search.DiscComputed(1,j);

