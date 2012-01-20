% zFileRedundancy(Filenames) explores possible redundancy between PDB files listed in Filenames.  It sorts the files by number of nucleotides, then compares files with similar numbers of nucleotides, performing a Needleman-Wunsch alignment of their bases.  If the percentage of bases which align exceeds the parameter p, the pair is kept for further examination.  It groups together structures connected by chains of greater than p percent sequence identify, then geometrically superimposes all pairs in each group, then prints a report so a human can decide which structures to keep.

% zFileRedundancy('Allfiles_list') % should give a huge report
% zFileRedundancy('NonRedundant_2008_02_21_list') % should show very little possible redundancy

function [ChosenNames] = zFileRedundancy(Filenames)

p = 0.90;                         % cutoff base match fraction
maxd = 0.5;

% ----------------------------------------- Read PDB lists, if any

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};                % make into a cell array
end

FullList = [];

for j=1:length(Filenames),
  FullList = [FullList; zReadPDBList(Filenames{j})];
end



% FullList = FullList(1:50);



load PDBInfo

i = [];

for f = 1:length(FullList),
  pp = find(ismember(upper(t(:,1)),upper(FullList{f})));
  if ~isempty(pp),
    if n(pp(1),1) > 0,            % avoid NMR structures, which have res. 0
      i = [i pp(1)];   
    end
  end
end
t = t(i,:);
n = n(i,:);           

i = find(n(:,2) > 1);             % restrict to structures with nucleotides
t = t(i,:);
n = n(i,:);           

[y,i] = sort(n(:,2));             % sort by number of nucleotides


% i = i(1:50);  % shorten!



t = t(i,:);
n = n(i,:);

F = length(i);                    % number of files

close = 0.5*sparse(eye(F));
prop  = 0.5*sparse(eye(F));

align{F,F} = [];

tim = cputime;

for i = 1:(F-1),
  if (cputime - tim > 60) || (mod(i,20) == 0),
    fprintf('Comparing structure %4d, which has %d nucleotides, to %d others of a similar size\n', i, n(i,2), length(find(p*n(i:F,2) < n(i,2))));
    tim = cputime;
  end
  j = i + 1;
  while (j <= F) && (p*n(j,2) < n(i,2)),  % while seqs COULD match well enough,
    [matches,a,b,ss,tt]   = dNeedlemanWunsch(t{i,9}, t{j,9}, 0.9999, 2);
    matches = sum(t{i,9}(a) == t{j,9}(b));
    pro = matches/min(n([i j],2));
    if ((n(i,2) - matches < 4) || (prop(i,j) > p)),
      close(i,j) = 1;
      prop(i,j) = pro;

%pro
%[ss;tt]

    end
    if prop(i,j) > 0.8,
      align{i,j} = [a; b];        % save alignment data
      align{j,i} = [b; a];
    end

    j = j + 1;
  end
end

close = close + close';
prop  = prop + prop';

closeseq = close;                 % store this for later

for k = 1:20,                     % Markov transitions for transitivity
  close = close * close;
end

figure(1)
clf
spy(close)
title(['Structures connected by a chain of more than ' num2str(p*100) '% similarity']);
drawnow

fprintf('Finding structures with more than %5.2f%% sequence similarity.\n',100*p);

done = zeros(1,F);                % whether each file has been considered

for i = 1:(F-1),
  if done(i) == 0,
    j = find(close(i,:));         % files that are "close" to i
    if length(j) < 2,             % no file is close to i
      done(i) = 1;
    else
%      if length(j) < 6,
%        [a,k] = sort(n(j,1));       % sort by resolution, best first
%      else
%        k = zOrderbySimilarity(-prop(j,j));  % sort by sequence similarity
%      end
%      j = j(k);

      File = zAddNTData(t(j,1),2);

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
          if nn == m,
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

      fprintf('Geometric discrepancy between aligned bases\n');
      fprintf('           ');

      for m = 1:length(j),
        fprintf(' %4s  ', File(m).Filename);
      end
      fprintf('\n');

      discmat = zeros(length(j));

      for m = 1:length(j),
        for nn = (m+1) : length(j),
          if ~isempty(align{j(m),j(nn)}),
            malign = align{j(m),j(nn)}(1,:);     % use stored data
            nalign = align{j(m),j(nn)}(2,:);
          else
            [matches,malign,nalign] = dNeedlemanWunsch(t{j(m),9},t{j(nn),9},0.9999,2);
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
          end     
        end
      end

      discmat = discmat + discmat';

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
          else
            fprintf(' %5.2f ', discmat(m,nn));
          end
        end
        fprintf('\n');
      end
      fprintf('\n');
      

    end
  end
end



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
fprintf('Final analysis\n');

close = closeseq;                 % 

for k = 1:20,                     % Markov transitions for transitivity
  close = close * close;
end

figure(2)
clf
spy(close)
title(['Structures connected by ' num2str(p*100) '% sequence similarity and ' num2str(maxd) ' discrepancy']);
drawnow

fprintf('Finding structures with more than %5.2f%% sequence similarity.\n',100*p);

done = zeros(1,F);                % whether each file has been considered
c = 0;                            % counts the number of files selected so far

for i = 1:(F-1),
  if done(i) == 0,
    j = find(close(i,:));         % files that are "close" to i
    if length(j) < 2,             % no file is close to i
      done(i) = 1;
      ChosenNames{c+1} = t{i,1};
      c = c + 1;
    else
%      if length(j) < 6,
%        [a,k] = sort(n(j,1));       % sort by resolution, best first
%      else
%        k = zOrderbySimilarity(-prop(j,j));  % sort by sequence similarity
%      end
%      j = j(k);

      File = zAddNTData(t(j,1),2);

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
          if nn == m,
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

      fprintf('Geometric discrepancy between aligned bases\n');
      fprintf('           ');

      for m = 1:length(j),
        fprintf(' %4s  ', File(m).Filename);
      end
      fprintf('\n');

      discmat = zeros(length(j));

      for m = 1:length(j),
        for nn = (m+1) : length(j),
          if ~isempty(align{j(m),j(nn)}),
            malign = align{j(m),j(nn)}(1,:);     % use stored data
            nalign = align{j(m),j(nn)}(2,:);
          else
            [matches,malign,nalign] = dNeedlemanWunsch(t{j(m),9},t{j(nn),9},0.9999,2);
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
          else
            fprintf(' %5.2f ', discmat(m,nn));
          end
        end
        fprintf('\n');
      end
      fprintf('\n');
      
      crit = [];
      for f = 1:length(j),
        crit(f) = datenum(File(f).Info.ReleaseDate, 'mm/dd/yyyy');
      end

      [y,k] = min(crit);

      fprintf('Chosen structure is %4s\n\n', File(k).Filename);

      ChosenNames{c+1} = File(k).Filename;
      c = c + 1;

    end
  end
end

for c = 1:length(ChosenNames),
  fprintf('%4s\n', ChosenNames{c});
end