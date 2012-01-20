% zFileRedundancy(Filenames) explores possible redundancy between PDB files listed in Filenames.  It sorts the files by number of nucleotides, then compares files with similar numbers of nucleotides, performing a Needleman-Wunsch alignment of their bases.  If the percentage of bases which align exceeds the parameter p, the pair is kept for further examination.  It groups together structures connected by chains of greater than p percent sequence identify, then geometrically superimposes all pairs in each group, then prints a report so a human can decide which structures to keep.

% zFileRedundancy('Allfiles_list') % should give a huge report
% zFileRedundancy('NonRedundant_2008_02_21_list') % should show very little possible redundancy

%function [ChosenNames] = zFileRedundancy(t,n,TodayDate)

%if nargin < 1,
%  load PDBInfo
%end 

diary Redundancy_Report_2009-05-13_Criterion_5.txt

Timeline = [];

if ~exist('TodayDate')
  TodayDate = datestr(date,'yyyy-mm-dd');
end

tot = cputime;                    % keep track of total time

p = 0.95;                         % cutoff base match fraction
maxd = 0.5;                       % cutoff discrepancy between structures
NTLimit = 7300;                   % above this limit, do not align sequences
MaxRes  = 4;                      % maximum resolution value to use
SL = 3000;                        % upper limit on # bases to align

Criterion = 5;                   
                                  % 1-earliest release date 
                                  % 2-resolution 
                                  % 3-number of nucleotides
                                  % 4-#pairs
                                  % 5-highest #BP / #nucleotides
                                  % 6-most recent release date

                                  % add 10, use preferred list to override

Preferred = zReadPDBList('Preferred_list',1);

JoinList = {'1FKA','1J5E'};       % structures that are redundant
%JoinList = [JoinList; {'1C2W','2AW4'}];

% --------------------------------- Look up information on these structures

load PDBInfo

i = find(n(:,2) > 0);             % restrict to structures with nucleotides
t = t(i,:);
n = n(i,:);           

fprintf('Found %4d structures containing RNA in the PDB.\n', length(t(:,1)));

i = find(n(:,2) > 1);             % restrict to structures with nucleotides
t = t(i,:);
n = n(i,:);           

fprintf('Found %4d structures with two or more nucleotides.\n', length(t(:,1)));

i = find(n(:,3) > 0);             % restrict to structures with basepairs
t = t(i,:);
n = n(i,:);           

fprintf('Found %4d structures with one or more basepairs.\n', length(t(:,1)));

i = find(n(:,1) > 0);             % omit structures with no resolution, esp NMR
t = t(i,:);
n = n(i,:);           

fprintf('Found %4d structures having a stated resolution.\n', length(t(:,1)));

i = find(n(:,1) <= MaxRes);       % restrict to structures with res <= MaxRes
t = t(i,:);
n = n(i,:);           

fprintf('Found %4d structures with resolution at or better than %6.2f\n', length(t(:,1)), MaxRes);

[y,i] = sort(n(:,2));             % sort by number of nucleotides
t = t(i,:);
n = n(i,:);


%i = find((n(:,2) > 1170) .* (n(:,2) < 1675));
%t = t(i,:);
%n = n(i,:);

%fprintf('Focusing on 16S structures today!\n');


%i = find(n(:,2) > 2000);       % big ribosomes only, for now
%i = find(n(:,2) < 40);       % big ribosomes only, for now
%t = t(i,:);
%n = n(i,:);

F = length(i);                    % number of files

fprintf('Found %4d 3D structures with resolution better than %6.2f and at least one RNA basepair.\n', F, MaxRes);

fprintf('Preparing a redundancy report on %d RNA 3D structures.\n',F);
fprintf('Sequence identity cutoff will be %7.2f.\n', p);
fprintf('Geometric discrepancy cutoff will be %7.2f.\n', maxd);

% ---------------------------------------------- Growth of whole database

for i = 1:F,
  d = datenum(t{i,4}, 'yyyy-mm-dd');
  Timeline = [Timeline; [d 1 n(i,2) n(i,3)]];  % accumulate data
end

[y,i] = sort(Timeline(:,1));                 % sort by date of increase
Timeline = Timeline(i,:);                    % re-order data

Date  = Timeline(:,1);
NumF  = Timeline(:,2);
NumNT = Timeline(:,3);
NumP  = Timeline(:,4);

Year = 1995 + (Date - datenum('01/01/1995','mm/dd/yyyy'))/365;

figure(1)
clf
subplot(3,1,1)
stairs(Year,cumsum(NumF));
title(['Total number of RNA 3D structures with resolution better than ' num2str(MaxRes)]);
axis([1992 2009 0 1.05*sum(NumF)]);

subplot(3,1,2)
stairs(Year,cumsum(NumNT));
title('Number of nucleotides in these structures');
axis([1992 2009 0 1.05*sum(NumNT)]);

subplot(3,1,3)
stairs(Year,cumsum(NumP));
title('Number of basepairs in these structures');
axis([1992 2009 0 1.05*sum(NumP)]);

RTimeline = Timeline;

% --------------------------------- Compare sequences between files

Close = 0.5*sparse(eye(F));       % indicator of whether sequences are close
prop  = 0.5*sparse(eye(F));

clear align
align{F,F} = [];

stim = cputime;
tim = cputime;

Linked = zeros(1,F);              % whether each file is already redundant

for i = 1:(F-1),
 if Linked(i) <= 1 || n(i,2) < 1300,
  if (cputime - tim > 60) || (mod(i,20) == 0) || (i==1),
    fprintf('Comparing structure %4d, %s which has %4d nucleotides, to %2d others of a similar size\n', i, t{i,1}, n(i,2), length(find(p*n(i:F,2) < n(i,2))));
    tim = cputime;
  end

  j = i + 1;                              % index of structure to compare to

  while (j <= F) && (p*length(t{j,11}) < length(t{i,11})),  % while seqs COULD match well enough,

    if n(i,2) <= NTLimit || length(t{i,8}) == 0 || length(t{j,8}) == 0,
      ti = t{i,11};                           % seq from non-red chains
      ti = ti(1:min(length(ti),SL));          % only compare first SL bases
      tj = t{j,11};                           % seq from non-red chains
      tj = tj(1:min(length(tj),SL));          % only compare first SL bases
      [matches,a,b,ss,tt] = zNeedlemanWunsch(ti,tj);
      e = (t{i,11}(a) == t{j,11}(b));         % locations of agreement
      matches = sum(e);
      pro = matches/min(length(ti),length(tj));          % percent agreement

%fprintf('%s %6.4f ', t{j,1}, pro);
drawnow

      if ((n(i,2) - matches < 4) || (pro > p)),
        Close(i,j) = 1;                       % sequences agree well enough
        prop(i,j)  = pro;
        Linked(j)  = Linked(j) + 1;
        if n(i,2) > NTLimit,
          fprintf('Sequence identity with %s is %6.4f\n', t{j,1}, pro);
        end
      end
      if prop(i,j) > 0.2,
        align{i,j} = [a(e); b(e)];        % save alignment data
        align{j,i} = [b(e); a(e)];
      end
    end
    j = j + 1;
  end
 end
%fprintf('\n');
end

fprintf('\nComparing sequences from structures took %7.1f minutes.\n\n', (cputime-stim)/60);

% ------------------------------------------ Join specified structures

for j = 1:length(JoinList(:,1)),
  a = find(ismember(t(:,1),JoinList{j,1}));
  b = find(ismember(t(:,1),JoinList{j,2}));
  Close(a,b) = 1;
  prop(a,b) = 1;
end


Close = Close + Close';           % extend by symmetry
prop  = prop + prop';

closeseq = Close;                 % store this for later

for k = 1:20,                     % Markov transitions for transitivity
  Close = Close * Close;          % extend by transitivity
end

figure(2)
clf
spy(Close)
title(['Structures connected by a chain of more than ' num2str(p*100) '% similarity']);
drawnow

% ---------------------------------------- Superimpose geometrically

fprintf('Superimposing structures with more than %5.2f%% sequence similarity.\n',100*p);

done = zeros(1,F);                % whether each file has been considered

stim = cputime;

for i = 1:(F-1),                  % loop through files
  if done(i) == 0,                % if file i hasn't been considered yet
    j = find(Close(i,:));         % files that are "close" to i
    if length(j) < 2,             % no file is close to i
      done(i) = 1;                % file i has a unique sequence
    elseif n(i,2) < NTLimit,
      File = zAddNTData(t(j,1));     % load nucleotide data
      File = zOrderChains(File);     % put largest chain first

      for m = 1:length(j),
        E  = abs(triu(File(m).Edge));
        np = full(sum(sum((E > 0) .* (E < 16))));
        fprintf('%4s has %4d nucleotides, %4d pairs, ', File(m).Filename, File(m).NumNT, np);
       if isempty(File(m).Info.Resolution),
         fprintf('resolution  ---- ');
       else
         fprintf('resolution %5.2f ', File(m).Info.Resolution);
       end

       Info = File(m).Info;

       fprintf(' %10s %s %s %s\n', Info.ReleaseDate, Info.Source, Info.Descriptor, Info.Author);
      end

      fprintf('Percent agreement of base sequence, using alignment of sequences\n');

      fprintf('           ');
      for m = 1:length(j),
        fprintf(' %4s  ', File(m).Filename);
      end
      fprintf('\n');

      for m = 1:length(j),
        done(j(m)) = 1;
        fprintf('%4s', File(m).Filename);
        if isempty(File(m).Info.Resolution),
          fprintf(' ---- ');
        else
          fprintf('%5.2f ', File(m).Info.Resolution);
        end
        for nn = 1:length(j),
          if nn == m || prop(j(m),j(nn)) == 0,
            fprintf('       ');
          else
            fprintf(' %5.1f ', 100*prop(j(m),j(nn)));
          end
        end
        fprintf('\n');
      end
      fprintf('\n');

      for m = 1:length(j),
        if length(File(m).NT) < 80,
          fprintf('%4s %s\n', File(m).Filename, t{j(m),9});
        end
      end
      fprintf('\n');

      discmat = zeros(length(j));

if length(j) > 30,
  fprintf('Starting to calculate discrepancies\n');
  drawnow
end

      for m = 1:length(j),
        for nn = (m+1) : length(j),
          if ~isempty(align{j(m),j(nn)}),
            malign = align{j(m),j(nn)}(1,:);     % use stored data
            nalign = align{j(m),j(nn)}(2,:);
          else
            malign = [];
            nalign = [];
          end

          if isempty(malign),
            d = Inf;
          else
            d = xDiscrepancy(File(m),malign,File(nn),nalign);
          end
          discmat(m,nn) = d;

          if ~(d <= maxd),                  % allow for d = NaN
            closeseq(j(m),j(nn)) = 0;       % these are not that close!
            closeseq(j(nn),j(m)) = 0;  

%            fprintf('Big discrepancy of %7.4f between %s and %s\n',d,t{j(m),1},t{j(nn),1});

          end     
        end
      end

      discmat = discmat + discmat';

      fprintf('Geometric discrepancy between aligned bases\n');
      fprintf('           ');

      for m = 1:length(j),
        fprintf(' %4s  ', File(m).Filename);
      end
      fprintf('\n');

      for m = 1:length(j),
        fprintf('%4s', File(m).Filename);
        if isempty(File(m).Info.Resolution),
          fprintf(' ---- ');
        else
          fprintf('%5.2f ', File(m).Info.Resolution);
        end
        for nn = 1:length(j),
          if nn == m,
            fprintf('       ');
          elseif discmat(m,nn) < Inf,
            fprintf(' %5.2f ', discmat(m,nn));
          else
            fprintf('       ');
          end
        end
        fprintf('\n');
      end
      fprintf('\n');
    end      
  end
  drawnow
end

fprintf('\nSuperimposing structures took %7.1f minutes.\n\n', (cputime-stim)/60);

% Now repeat the analysis, having removed links between structures that are
% close in sequence but don't superimpose well geometrically.


fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n');
fprintf('Final analysis\n\n');

Timeline = [];
clear ChosenNames
clear Text

Close = closeseq;                 % 

% ------------------------------------------ Join specified structures

for j = 1:length(JoinList(:,1)),
  a = find(ismember(t(:,1),JoinList{j,1}));
  b = find(ismember(t(:,1),JoinList{j,2}));
  Close(a,b) = 1;
  prop(a,b) = 1;
end

for k = 1:20,                     % Markov transitions for transitivity
  Close = Close * Close;
end

figure(3)
clf
spy(Close)
title(['Structures connected by ' num2str(p*100) '% sequence similarity and ' num2str(maxd) ' discrepancy']);
drawnow

fprintf('Listing structures with more than %5.2f%% sequence similarity and less than %5.2f geometric discrepancy.\n',100*p,maxd);

fprintf('Structures are sorted by ');
switch mod(Criterion,10)
  case 1, fprintf('release date, earliest first\n');
  case 2, fprintf('resolution, best resolution first\n');
  case 3, fprintf('number of nucleotides (decreasing)\n');
  case 4, fprintf('number of basepairs (decreasing)\n');
  case 5, fprintf('number of basepairs per nucleotide (decreasing)\n');
  case 6, fprintf('release date, most recent first\n');
end

fprintf('The first structure listed will be chosen as the representative');
if Criterion > 10,
  fprintf(' unless one of these structures is found in the preferred list of structures.\n');
else
  fprintf('.\n');
end

fprintf('\n');

done = zeros(1,F);                % whether each file has been considered
c = 0;                            % counts the number of files selected so far

for i = 1:(F-1),                  % loop through all files
  if done(i) == 0,                % file does not already appear in a report
    j = find(Close(i,:));         % files that are "close" to i
    if length(j) < 2,             % no file is close to i
      done(i) = 1;                % no need to display this one again
      File = zAddNTData(t(i,1)); % load this file
      E  = abs(triu(File(1).Edge));
      np(1) = full(sum(sum((E > 0) .* (E < 16)))); % number of pairs
    else
       
      crit = [];                   % criterion for sorting and choosing
      np   = [];

      for f = 1:length(j),
        done(j(f)) = 1;            % no need to display this one again
%        E  = abs(triu(File(f).Edge));
%        np(f) = full(sum(sum((E > 0) .* (E < 16)))); % number of pairs
        np(f) = n(j(f),3);
        switch mod(Criterion,10)
          case 1, crit(f) = datenum(t{j(f),4}, 'mm/dd/yyyy');
          case 2, crit(f) = n(j(f),1);            % resolution
          case 3, crit(f) = -n(j(f),2);           % number of nucleotides
          case 4, crit(f) = -np(f);               % number of pairs
          case 5, crit(f) = -n(j(f),3)/n(j(f),2); % ratio # bp / # nucleotides
        end
      end

      fprintf('\n');

      File = zAddNTData(t(j,1)); % load these files
      File = zOrderChains(File);

      [y,k] = sort(crit);
      j     = j(k);                         % order by the criterion
      File  = File(k);
      np    = np(k);

      for m = 1:length(j),
        fprintf('%4s has %4d nucleotides, %4d pairs, ', t{j(m),1}, n(j(m),2), n(j(m),3));
       if isempty(File(m).Info.Resolution),
         fprintf('resolution  ---- ');
       else
         fprintf('resolution %5.2f ', File(m).Info.Resolution);
       end

       Info = File(m).Info;

       fprintf(' %10s | %s | %s | %s\n', Info.ReleaseDate, Info.Source, Info.Descriptor, Info.Author);
      end

      if n(i,2) <= NTLimit,

        fprintf('Percent agreement of base sequence, using alignment of sequences\n');
        fprintf('           ');

        for m = 1:length(j),
          fprintf(' %4s  ', File(m).Filename);
        end
        fprintf('\n');

        for m = 1:length(j),
          fprintf('%4s', File(m).Filename);
          if isempty(File(m).Info.Resolution),
            fprintf(' ---- ');
          else
            fprintf('%5.2f ', File(m).Info.Resolution);
          end
          for nn = 1:length(j),
            if nn == m || prop(j(m),j(nn)) == 0,
              fprintf('       ');
            else
              fprintf(' %5.1f ', 100*prop(j(m),j(nn)));
            end
          end
          fprintf('\n');
        end
        fprintf('\n');

        for m = 1:length(j),
          if length(File(m).NT) < 80,
            fprintf('%4s %s\n', File(m).Filename, t{j(m),9});
          end
        end
        fprintf('\n');

        discmat = zeros(length(j));

        for m = 1:length(j),
          for nn = (m+1) : length(j),
            if ~isempty(align{j(m),j(nn)}),
              malign = align{j(m),j(nn)}(1,:);     % use stored data
              nalign = align{j(m),j(nn)}(2,:);
            else
              malign = [];
              nalign = [];
%              [matches,malign,nalign] = zNeedlemanWunsch(t{j(m),9},t{j(nn),9});
            end

            if isempty(malign),
              d = Inf;
            else
              d = xDiscrepancy(File(m),malign,File(nn),nalign);
            end
            discmat(m,nn) = d;

          end
        end

        discmat = discmat + discmat';

        fprintf('Geometric discrepancy between aligned bases\n');
        fprintf('           ');

        for m = 1:length(j),
          fprintf(' %4s  ', File(m).Filename);
        end
        fprintf('\n');

        for m = 1:length(j),
          fprintf('%4s', File(m).Filename);
          if isempty(File(m).Info.Resolution),
            fprintf(' ---- ');
          else
            fprintf('%5.2f ', File(m).Info.Resolution);
          end
          for nn = 1:length(j),
            if nn == m || discmat(m,nn) == Inf,
              fprintf('       ');
            else
              fprintf(' %5.2f ', discmat(m,nn));
            end
          end
          fprintf('\n');
        end
      end
    end

    ff = 1;                               % default file to use

    if Criterion > 10,                    % replace with preferred name
      pp = [];
      for f = 1:length(j),
        p = find(ismember(Preferred,File(f).Filename));
        if ~isempty(p),
          pp = p;
          ff = f;
          fprintf('Found %s in the list of preferred structures\n', File(f).Filename);
          cn = File(f).Filename;
        end
      end
    end

    c = c + 1;
    ChosenNames{c} = File(ff).Filename;

    Equivalents{c,1} = ChosenNames{c};
    for m = 1:length(j),
      Equivalents{c,m+1} = t{j(m),1};        % list these as being equivalent
    end

    if length(File) == 1,
      Text{c} = sprintf('Unique structure is %4s, which ', ChosenNames{c});
    else
      Text{c} = sprintf('Chosen structure is %4s, which ', ChosenNames{c});
    end

    Text{c} = [Text{c} sprintf('has %4d nucleotides, %4d pairs, ', File(ff).NumNT, np(ff))];
    if isempty(File(ff).Info.Resolution),
      Text{c} = [Text{c} sprintf('resolution  ---- ')];
    else
      Text{c} = [Text{c} sprintf('resolution %5.2f ', File(ff).Info.Resolution)];
    end

    Info = File(ff).Info;

    Text{c} = [Text{c} sprintf(' %10s | %s | %s | %s', Info.ReleaseDate, Info.Source, Info.Descriptor, Info.Author)];

    fprintf('%s\n', Text{c});

    if Criterion == 1,                            % date of deposition
      maxNT = length(File(1).NT);
      maxNP = np(1);
      Timeline = [Timeline; [y(1) 1 maxNT maxNP]];
      for m = 1:length(j),
        newNT = length(File(m).NT);
        newNP = np(m);
        if newNT > maxNT,
          Timeline = [Timeline; [y(m) 0 (newNT-maxNT) max(0,newNP-maxNP)]];
          maxNT = newNT;
          maxNP = max(maxNP,newNP);
        end
        if newNP > maxNP,
          Timeline = [Timeline; [y(m) 0 max(0,newNT-maxNT) (newNP-maxNP)]];
          maxNT = newNT;
          maxNP = max(maxNP,newNP);
        end
      end
    end

  end
end

% ----------------------------------------------- List chosen files

fprintf('\nList of chosen files:\n');

for c = 1:length(ChosenNames),
  fprintf('%4s\n', ChosenNames{c});
end

fid = fopen(['PDBFiles' filesep 'Nonredundant_' TodayDate '_list.pdb'],'w');
for c = 1:length(ChosenNames),
  fprintf(fid,'%4s\n', ChosenNames{c});
end
fclose(fid);
fprintf('\nChosen files were written to %s\n', ['PDBFiles' filesep 'Nonredundant_' TodayDate '_list.pdb']);

fprintf('\nInformation about chosen files:\n');

for c = 1:length(ChosenNames),
  fprintf('%4s\n', Text{c});
end

% ---------------------------------------------- Graphs of growth of database

if Criterion == 1,
  [y,i] = sort(Timeline(:,1));                 % sort by date of increase
  Timeline = Timeline(i,:);                    % re-order data

  Date  = Timeline(:,1);
  NumF  = Timeline(:,2);
  NumNT = Timeline(:,3);
  NumP  = Timeline(:,4);

  Year = 1995 + (Date - datenum('01/01/1995','mm/dd/yyyy'))/365;

  figure(4)
  clf
  subplot(3,1,1)
  stairs(Year,cumsum(NumF));
  title('Number of distinct RNA 3D structures');
  axis([1992 2010 0 1.05*sum(NumF)]);

  subplot(3,1,2)
  stairs(Year,cumsum(NumNT));
  title('Number of nucleotides in distinct structures');
  axis([1992 2010 0 1.05*sum(NumNT)]);
  set(gca,'YTick',[0 10000 20000 30000 40000])
  set(gca,'YTickLabel',{'0','10000','20000','30000','40000'});

  subplot(3,1,3)
  stairs(Year,cumsum(NumP));
  title('Number of basepairs in distinct structures');
  axis([1992 2010 0 1.05*sum(NumP)]);
end

% ----------------------------------------- Update PDBInfo with this list

if Criterion == 15,

fprintf('\nUpdating PDB Info with non-redundant list and equivalencies\n');

load PDBInfo

N = lower(t(:,1));

for c = 1:length(ChosenNames),                     % chosen structures
  j = 2;                                           % start in second column
  for j = 1:length(Equivalents(1,:)),
   if ~isempty(Equivalents{c,j}),                   % look for equivalents
    r = find(ismember(N,lower(Equivalents{c,j}))); % find in list t
    if length(r) > 1,
      fprintf('Multiple hits in PDBInfo for %s.\n', Equivalents{c,j});
    end
    for rr = 1:length(r),
      t{r(rr),10} = Equivalents{c,1};               % add a column to repres
    end
   end
  end
end

save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7

end

fprintf('Total elapsed time %8.6f minutes.\n', (cputime-tot)/60);



diary off

break
return

% ----------------------------------------- Read PDB structure names and lists
Filenames = 'Allfiles_list';

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};                % make into a cell array
end

FullList = [];

for j=1:length(Filenames),
  FullList = [FullList; zReadPDBList(Filenames{j})];
end

% FullList = FullList(1:50);

