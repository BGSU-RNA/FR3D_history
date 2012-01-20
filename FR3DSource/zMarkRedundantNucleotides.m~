% zMarkRedundantNucleotides(File) identifies redundant chains in File and sets File.Redundant(i,j) = 1 if nucleotides with indices i and j, from different chains, are redundant when these chains are aligned and superimposed

function [File,LongestChain] = zMarkRedundantNucleotides(File, Verbose)

if nargin < 2,
  Verbose = 1;
end

DiscCutoff = 0.4;                               % limit on how dissimilar

for f = 1:length(File),
 File(f).Redundant = sparse(zeros(length(File(f).NT)));
 
 if length(File(f).NT) > 0,
  Chain = cat(2,File(f).NT.Chain);              % all chain identifiers
  U = unique(Chain);                            % unique chain identifiers

  if length(U) > 1,                             % more than one chain
    leng = [];
    clear bases indic

    for u = 1:length(U),                        % loop through chains
      i = find(Chain == U(u));                  % indices of this chain
      bases{u} = cat(2,File(f).NT(i).Base);     % store bases of this chain
      indic{u} = i;                             % store indices of this chain
      leng(u)  = length(i);                     % store length of this chain
    end

    if Verbose > 0,
      for u = 1:length(U),
        fprintf('%s Chain %s:  %s\n', File(f).Filename, U(u), bases{u});
      end
    end

    [y,i] = max(leng);
    LongestChain{f} = U(i(1));                  % chain with greatest length

    % ------------------------------------------- compare chains

    if Verbose > 0,
      fprintf('%s chains with redundancy: ', File(f).Filename);
    end

    for u = 1:length(U),                        % loop through chains
      for v = (u+1):length(U),                  % loop through chains

        % --------------------------------------- align sequences
        [matches,a,b,ss,tt] = dNeedlemanWunsch(bases{u},bases{v});
        e = find(bases{u}(a) == bases{v}(b));   % locations of agreement
        matches = length(e);                    % number of agreements

        if matches > 1,                         % if there are any matches
          a = a(e);                             % focus on the matches
          b = b(e);

          [Disc,R,MM,CM,A] = xDiscrepancy(File(f),indic{u}(a),File(f),indic{v}(b));

          c = find(abs(A) <= 0.8);              % reasonably similar bases

          if Verbose > 2,
            fprintf('\n%s\n%s\n\n', bases{u}(a), bases{v}(b));

            fprintf('%s\n%s\n', bases{u}(a(c)), bases{v}(b(c)));
          end

          if length(c) > 2 && length(c) > length(a)/2, 
            [Disc,R,MM,CM,A] = xDiscrepancy(File(f),indic{u}(a(c)),File(f),indic{v}(b(c)));


            E = eye(length(c),length(c));

            if Disc <= DiscCutoff,
              File(f).Redundant(indic{u}(a(c)),indic{v}(b(c))) = E;
              File(f).Redundant(indic{v}(b(c)),indic{u}(a(c))) = E;

              if Verbose > 0,
                fprintf('%c%c ', U(u), U(v));
              end

            end
          end
        end
      end
    end

    if Verbose > 0,
      fprintf('\n\n');
    end

    if Verbose > 1,
      figure(1)
      clf
      zCircularDiagram(File(f),1);
      figure(2)
      clf
      spy(File(f).Redundant)
      drawnow
    end

  else
    LongestChain{f} = U(1);                     % only one chain
  end
 end
end

