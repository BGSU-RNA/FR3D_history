% zFileRedundancy(Filenames) explores possible redundancy between PDB files listed in Filenames.  It sorts the files by number of nucleotides, then compares files with similar numbers of nucleotides, performing a Needleman-Wunsch alignment of their bases.  If the percentage of bases which align exceeds the parameter p, the pair is kept for further examination.  It groups together structures connected by chains of greater than p percent sequence identify, then geometrically superimposes all pairs in each group, then prints a report so a human can decide which structures to keep.

% zFileRedundancy('Allfiles_list') % should give a huge report
% zFileRedundancy('NonRedundant_2008_02_21_list') % should show very little possible redundancy

%function [ChosenNames] = zFileRedundancy(t,n,TodayDate)

load PDBInfo

diary Redundancy_Report_2010-05-19_b_Criterion_5.txt

Timeline = [];

if ~exist('TodayDate')
  TodayDate = datestr(date,'yyyy-mm-dd');
end

tot = cputime;                    % keep track of total time

p = 0.95;                         % cutoff base match fraction
maxd = 0.5;                       % cutoff discrepancy between structures
NTLimit = 30000;                  % above this limit, do not align sequences
MaxRes  = 100;                    % maximum resolution value to use
SL = 4000;                        % upper limit on # bases to align

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

for i = 1:length(t(:,1)),
  n(i,4) = length(t{i,11});       % length of longest chain
end

i = find(n(:,2) > 0);             % restrict to structures with nucleotides
t = t(i,:);
n = n(i,:);           

fprintf('Found %4d structures containing RNA in the PDB.\n', length(t(:,1)));

i = find(n(:,4) > 1);             % restrict to structures >= 2 nucleotides
t = t(i,:);
n = n(i,:);           

fprintf('Found %4d structures with two or more nucleotides.\n', length(t(:,1)));

if 0 > 1,
i = find(n(:,3) > 0);             % restrict to structures with basepairs
t = t(i,:);
n = n(i,:);           

fprintf('Found %4d structures with one or more basepairs.\n', length(t(:,1)));
end

if 0 > 1,
  i = find(n(:,1) > 0);           % omit structures with no resolution, esp NMR
  t = t(i,:);
  n = n(i,:);           
 fprintf('Found %4d structures having a stated resolution.\n', length(t(:,1)));
end

if 0 > 1,
  i = find(n(:,1) <= MaxRes);       % restrict to structures with res <= MaxRes
  t = t(i,:);
  n = n(i,:);           
fprintf('Found %4d structures with resolution at or better than %6.2f\n', length(t(:,1)), MaxRes);
end

if max(n(:,2)) > NTLimit,
  i = find(n(:,2) <= NTLimit);
  fprintf('Found %4d structures with %d nucleotides or fewer in their longest chain\n', length(i), NTLimit);
  t = t(i,:);
  n = n(i,:);
end

[y,i] = sort(n(:,4));         % sort by number of nucleotides in longest chain
t = t(i,:);
n = n(i,:);

if 0 > 1,
  i = find((n(:,2) > 1170) .* (n(:,2) < 1675));
  t = t(i,:);
  n = n(i,:);
  fprintf('Focusing on 16S structures today!\n');
end

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

for i = 1:(F-1),                          % loop through files
 if Linked(i) <= 1 || n(i,4) < 1300,      % if already linked 0 or 1 times
                                          % or if a small structure
  J = find( ((1:F)' > i) .* ((p*n(1:F,4) < n(i,4)) + (n(1:F,4) - 4 <= n(i,4))));
  if (cputime - tim > 60) || (mod(i,20) == 0) || (i==1),

    fprintf('Comparing sequence from structure %4d, %s which has %4d nucleotides in its longest chain, to %3d others of a similar size\n', i, t{i,1}, n(i,4), length(J));
    tim = cputime;

  end

  for k = 1:length(J),
    j = J(k);
    if n(i,2) <= NTLimit || length(t{i,8}) == 0 || length(t{j,8}) == 0,
      ti = t{i,11};                           % seq from longest chain
      ti = ti(1:min(length(ti),SL));          % only compare first SL bases
      tj = t{j,11};                           % seq from longest chain
      tj = tj(1:min(length(tj),SL));          % only compare first SL bases
      [matches,a,b,ss,tt] = zNeedlemanWunsch(ti,tj); % sequence alignment
      e = (t{i,11}(a) == t{j,11}(b));         % locations of agreement
      matches = sum(e);                       % number of exact matches
      pro = matches/min(length(ti),length(tj));    % percent identity

%fprintf('%s %6.4f ', t{j,1}, pro);
drawnow

      if (min(length(ti),length(tj)) - matches <= 4) || (pro > p),
        Close(i,j) = 1;                       % sequences agree well enough
        prop(i,j)  = pro;
        Linked(j)  = Linked(j) + 1;           % tally how many times linked
        if n(i,2) > NTLimit,
          fprintf('Sequence identity with %s is %6.4f\n', t{j,1}, pro);
        end
      end
      if prop(i,j) > 0.8,                 % even a little sequence agreement
        align{i,j} = [a(e); b(e)];        % save alignment data
        align{j,i} = [b(e); a(e)];
      end
    end
    j = j + 1;
  end
 end

end

fprintf('\nComparing sequences from structures took %7.1f minutes.\n\n', (cputime-stim)/60);

% ------------------------------------------ Join specified structures

fprintf('\nExtending by transitivity ');

ttime = cputime;

for j = 1:length(JoinList(:,1)),
  a = find(ismember(t(:,1),JoinList{j,1}));
  b = find(ismember(t(:,1),JoinList{j,2}));
  Close(a,b) = 1;
  prop(a,b) = 1;
end


Close = Close + Close';           % extend by symmetry
prop  = prop + prop';

closeseq = Close;                 % store this for later

% Close = closeseq ^ 3;           % structures connected by a short chain 

for k = 1:log2(F),                % Markov transitions for transitivity
  Close = Close * Close;          % extend by transitivity
end

fprintf('took %7.1f minutes.\n\n', (cputime-ttime)/60);

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
    else,
      File = zAddNTData(t(j,1));     % load nucleotide data

      for m = 1:length(j),
        E  = abs(triu(File(m).Edge));
        fprintf('%4s has %4d nucleotides, %4d in longest chain, %4d basepairs, ', File(m).Filename, File(m).NumNT, n(j(m),4), n(j(m),3));
       if isempty(File(m).Info.Resolution),
         fprintf('resolution  ---- ');
       else
         fprintf('resolution %5.2f ', File(m).Info.Resolution);
       end

       Info = File(m).Info;

       fprintf(' %10s | %s | %s | %s\n', Info.ReleaseDate, Info.Source, Info.Descriptor, Info.Author);
      end

      fprintf('Percent agreement of base sequence, using alignment of sequences\n');

      fprintf('           ');
      for m = 1:length(j),
        fprintf(' %4s  ', File(m).Filename);
      end
      fprintf('\n');

      for m = 1:length(j),
        done(j(m)) = 1;
        T = sprintf('%4s', File(m).Filename);
        if isempty(File(m).Info.Resolution),
          T = [T sprintf(' ---- ')];
        else
          T = [T sprintf('%5.2f ', File(m).Info.Resolution)];
        end
        for nn = 1:length(j),
          if nn == m || prop(j(m),j(nn)) == 0,
            T = [T sprintf('       ')];
          else
            T = [T sprintf(' %5.1f ', 100*prop(j(m),j(nn)))];
          end
        end
        fprintf('%s\n',T);
      end
      fprintf('\n');

      fprintf('Sequences of longest chains:\n');
      for m = 1:length(j),                        % print sequences
        if length(t{j(m),11}) < 10000,
          fprintf('%4s %s\n', File(m).Filename, t{j(m),11});
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
            d = xDiscrepancy(File(m),File(m).LongestChain(1)-1+malign,File(nn),File(nn).LongestChain(1)-1+nalign);
          end
          discmat(m,nn) = d;

          if ~(d <= maxd),                  % allow for d = NaN
            closeseq(j(m),j(nn)) = 0;       % these are not that close!
            closeseq(j(nn),j(m)) = 0;  
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
        T = sprintf('%4s', File(m).Filename);
        if isempty(File(m).Info.Resolution),
          T = [T sprintf(' ---- ')];
        else
          T = [T sprintf('%5.2f ', File(m).Info.Resolution)];
        end
        for nn = 1:length(j),
          if nn == m,
            T = [T sprintf('       ')];
          elseif discmat(m,nn) < Inf,
            T = [T sprintf(' %5.2f ', discmat(m,nn))];
          else
            T = [T sprintf('       ')];
          end
        end
        fprintf('%s\n',T);
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

for k = 1:log2(F),                     % Markov transitions for transitivity
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
      File = zAddNTData(t(i,1));  % load this file
      np(1) = n(i,3);      
    else
       
      crit = [];                   % criterion for sorting and choosing
      np   = [];

      for f = 1:length(j),
        done(j(f)) = 1;            % no need to display this one again
        np(f) = n(j(f),3);
                                   % sorting criteria
        crit(f,1) = datenum(t{j(f),4}, 'yyyy-mm-dd');
        crit(f,2) = n(j(f),1);            % resolution
        crit(f,3) = -n(j(f),2);           % number of nucleotides
        crit(f,4) = -np(f);               % number of pairs
        crit(f,5) = -n(j(f),3)/n(j(f),2); % ratio # bp / # nucleotides
%        crit(f,5) = round(crit(f,5)*50);  % round to nearest 0.02
      end

      fprintf('\n');

      File = zAddNTData(t(j,1)); % load these files

      switch mod(Criterion,10)
        case 1, [y,k] = sortrows(crit,[1 5 2 3 4]);
        case 2, [y,k] = sortrows(crit,[2 5]);
        case 3, [y,k] = sortrows(crit,[3 2]);
        case 4, [y,k] = sortrows(crit,[4 5]);
        case 5, [y,k] = sortrows(crit,[5 2 1]);
      end

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
          T = sprintf('%4s', File(m).Filename);
          if isempty(File(m).Info.Resolution),
            T = [T sprintf(' ---- ')];
          else
            T = [T sprintf('%5.2f ', File(m).Info.Resolution)];
          end
          for nn = 1:length(j),
            if nn == m || prop(j(m),j(nn)) == 0,
              T = [T sprintf('       ')];
            else
              T = [T sprintf(' %5.1f ', 100*prop(j(m),j(nn)))];
            end
          end
          fprintf('%s\n',T);
        end
        fprintf('\n');

        fprintf('Sequences of longest chains\n');
        for m = 1:length(j),
          if length(t{j(m),11}) < 10000,
            fprintf('%4s %s\n', File(m).Filename, t{j(m),11});
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
%           [matches,malign,nalign] = zNeedlemanWunsch(t{j(m),11},t{j(nn),11});
            end

            if isempty(malign),
              d = Inf;
            else
              d = xDiscrepancy(File(m),File(m).LongestChain(1)-1+malign,File(nn),File(nn).LongestChain(1)-1+nalign);
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
          T = sprintf('%4s', File(m).Filename);
          if isempty(File(m).Info.Resolution),
            T = [T sprintf(' ---- ')];
          else
            T = [T sprintf('%5.2f ', File(m).Info.Resolution)];
          end
          for nn = 1:length(j),
            if nn == m || discmat(m,nn) == Inf,
              T = [T sprintf('       ')];
            else
              T = [T sprintf(' %5.2f ', discmat(m,nn))];
            end
          end
          fprintf('%s\n',T);
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

    Text{c} = [Text{c} sprintf('has %4d nucleotides, %4d in longest chain, %4d pairs, ', File(ff).NumNT, n(j(ff),4), np(ff))];
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

if Criterion == 5,

fprintf('\nUpdating PDBInfo.mat with non-redundant list and equivalencies\n');

load PDBInfo

N = lower(t(:,1));

for c = 1:length(t(:,1)),
  t{c,10} = '';
end

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

fprintf('Non-redundant list has %d entries.\n', length(ChosenNames));
fprintf('Total elapsed time %8.6f minutes.\n', (cputime-tot)/60);

diary off

% zWriteNRFiles           % put a copy of .mat files in NonRedundantMatFiles
