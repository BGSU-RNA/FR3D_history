% zFileRedundancy(Filenames) explores possible redundancy between PDB files listed in Filenames.  It sorts the files by number of nucleotides, then compares files with similar numbers of nucleotides, performing a Needleman-Wunsch alignment of their bases.  If the percentage of bases which align exceeds the parameter p, the pair is kept for further examination.  It groups together structures connected by chains of greater than p percent sequence identify, then geometrically superimposes all pairs in each group, then prints a report so a human can decide which structures to keep.

% zFileRedundancy('Allfiles_list') % should give a huge report
% zFileRedundancy('NonRedundant_2008_02_21_list') % should show very little possible redundancy

function [void] = zFileRedundancy(Filenames)

p = 0.70;                         % cutoff base match fraction

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
    i = [i pp(1)];
  end
end
t = t(i,:);
n = n(i,:);           

i = find(n(:,2) > 1);             % restrict to structures with nucleotides
t = t(i,:);
n = n(i,:);           

[y,i] = sort(n(:,2));             % sort by number of nucleotides
t = t(i,:);
n = n(i,:);

F = length(i);                    % number of files

close = 0.5*sparse(eye(F));
prop  = 0.5*sparse(eye(F));

for i = 1:(F-1),
  if mod(i,20) == 0,
    fprintf('Comparing structure %4d to others of a similar size\n', i);
  end
  j = i + 1;
  while (j <= F) && (p*n(j,2) < n(i,2)),  % while seqs COULD match well enough,
    [matches,a,b,ss,tt]   = dNeedlemanWunsch(t{i,9}, t{j,9}, 0.9999, 2);
    matches = sum(t{i,9}(a) == t{j,9}(b));
    pro = matches/min(n([i j],2));
    if ((n(i,2) - matches < 4) || (prop(i,j) > p)),
      close(i,j) = 1;
      prop(i,j) = pro;
      align{i,j} = [a; b];        % save alignment data
      align{j,i} = [a; b];

%pro
%[ss;tt]

    end
    j = j + 1;
  end
end

close = close + close';
prop  = prop + prop';

for k = 1:20,                     % Markov transitions for transitivity
  close = close * close;
end

spy(close)
drawnow

fprintf('Finding structures with more than %5.2f%% sequence similarity.\n',100*p);

done = zeros(1,F);                % whether each file has been considered

for i = 1:(F-1),
  if done(i) == 0,
    j = find(close(i,:));         % files that are "close" to i
    if length(j) < 2,             % no file is close to i
      done(i) = 1;
    else
      if length(j) < 6,
        [a,k] = sort(n(j,1));       % sort by resolution, best first
      else
        k = zOrderbySimilarity(-prop(j,j));  % sort by sequence similarity
      end
      j = j(k);

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

      fprintf('Geometric discrepancy between aligned bases\n');
      fprintf('           ');

      for m = 1:length(j),
        fprintf(' %4s  ', File(m).Filename);
      end
      fprintf('\n');

      discmat = zeros(length(j));

      for m = 1:length(j),
        for nn = (m+1) : length(j),

          malign = align{j(m),j(nn)}(1,:);     % use stored data
          nalign = align{j(m),j(nn)}(2,:);

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
      

    end
  end
end

