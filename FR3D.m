% FR3D conducts the search given in xSpecifyQuery in the PDB files listed
% Change the list of PDB files to be searched by editing zFileNameList
% Change the query by editing xSpecifyQuery

Search.SaveName = datestr(now,31);  % use date and time to identify this search

if ~exist('GUIactive'),

  % ----------------------------------------- Specify PDB files to search -

  Filenames = zFileNameList;                   % PDB Filenames to search

  % ------------------------------------------- Read search parameters --------

  Query = xSpecifyQuery;                       % get search parameters

end

% ----------------------------------------- Load PDB files if needed --------

if ~exist('File'),
  [File,SIndex] = zAddNTData(Filenames,0);   % load PDB data
else
  [File,SIndex] = zAddNTData(Filenames,0,File); % add PDB data  
end

% ------------------------------------------- Construct details of search ---

if isfield(Query,'Filename'),                 % if query motif is from a file
  [File,QIndex] = zAddNTData(Query.Filename,0,File);  
                                              % load data for Query, if needed
  Query = xConstructQuery(Query,File(QIndex));
                                              % preliminary calculations
else
  Query = xConstructQuery(Query);             % preliminary calculations
end

if isfield(Query,'NumNT'),

% ------------------------------------------- Display query information------

fprintf('Query %s:', Query.Name);

if isfield(Query,'Description'),
  fprintf(' %s\n', Query.Description);
else
  fprintf('\n');
end

% ------------------------------------------- Calc more distances if needed -

if Query.Geometric > 0,                       % if a geometric search
  tic
  CalcFlag = 0;
  for f=1:length(SIndex),
    i = SIndex(f);
    if ceil(Query.DistCutoff) > ceil(max(max(File(i).Distance))),
      c        = cat(1,File(i).NT(1:File(i).NumNT).Center);
      File(i).Distance = zMutualDistance(c,Query.DistCutoff); 
             % sparse matrix of center-center distances, up to Query.DistCutoff
      if length(File(i).NT) > 10,
        zSaveNTData(File(i));
        drawnow;
      end
      CalcFlag = 1;
    end
  end
  if CalcFlag > 0,
    fprintf('Calculated more distances in %5.3f seconds\n', toc);
  end
end

drawnow

% ------------------------------------------- Find candidates ---------------

starttime = cputime;

Candidates = xFindCandidates(File(SIndex),Query);  % find candidates

if ~isempty(Candidates),
 if Query.Geometric > 0,
  [Discrepancy, Candidates] = xRankCandidates(File(SIndex),Query,Candidates);
  fprintf('Found %d candidates in the desired discrepancy range\n',length(Discrepancy));

   if Query.ExcludeOverlap > 0,
     [Candidates, Discrepancy] = xReduceOverlap(Candidates,Discrepancy); 
                                                 % quick reduction in number
     [Candidates, Discrepancy] = xExcludeOverlap(Candidates,Discrepancy,400); 
                                                % find top 400 distinct ones
     fprintf('Removed highly overlapping candidates, kept %d\n', length(Candidates(:,1)));
   end

 else
  A = [Candidates sum(Candidates')'];        % compute sum of indices
  N = Query.NumNT;                           % number of nucleotides
  [y,i] = sortrows(A,[N+1 N+2 1:N]);         % sort by this sum
  Candidates = Candidates(i,:);              % put all permutations together
  Discrepancy = zeros(length(Candidates(:,1)),1);
 end

% -------------------------------------------------- Save results of search

 Search.Query       = Query;
 Search.Filenames   = Filenames;
 Search.TotalTime   = cputime - starttime;
 Search.Date        = Search.SaveName(1:10);
 Search.Time        = Search.SaveName(12:18);
 Search.SaveName    = strrep(Search.SaveName,' ','_');
 Search.SaveName    = strrep(Search.SaveName,':','_');
 Search.Candidates  = Candidates;
 Search.Discrepancy = Discrepancy;

 if ~(exist('SearchSaveFiles') == 7),        % if directory doesn't yet exist
   mkdir('SearchSaveFiles');
 end

 save(['SearchSaveFiles' filesep Search.SaveName], 'Search');

% ------------------------------------------------ Display results

 fprintf('Entire search took %8.4f seconds, or %8.4f minutes\n', (cputime-starttime), (cputime-starttime)/60);

if ~exist('GUIactive'),
  xListCandidates(File(SIndex),Search,1);
  xDisplayCandidates(File(SIndex),Search);
%  xGroupCandidates(File(SIndex),Search);  % doesn't work very well yet!
end

end

end