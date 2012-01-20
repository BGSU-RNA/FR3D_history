% zFixPairs identifies nucleotides which have been classified as using the same edge in more than one basepair and chooses the best basepair among them

% function [File] = zFixPairs(File,Verbose)

% File = zAddNTData('2avy');
% File = zAddNTData('Nonredundant_2009-05-14_list');

% Note!  This program does not actually change any classifications!

Verbose = 2;                  % show the nucleotides and wait for key press
Verbose = 1;                  % list the offending pairs

load PairExemplars

for f = 1:length(File),

  [i,j,e] = find(File(f).Edge .* (abs(File(f).Edge) < 13));
  u = zeros(size(e));

  for k = 1:length(i),
    t = zEdgeText(e(k));
    switch upper(t(2))
    case   'W', u(k) = 1;     % WC edge
    case   'H', u(k) = 2;     % Hoogsteen edge
    case   'S', u(k) = 3;     % sugar edge
    end
  end

  w = sparse(i,j,u);

  for b = 1:3,                          % which edge we are checking
    i = find(sum(w'==b) > 1);           % nucleotides using an edge twice
    for a = 1:length(i),
      j = find(w(i(a),:) == b);         % nucleotides that i(a) interacts with
      NT1 = File(f).NT(i(a));
      d = [];
      for c = 1:length(j),
        NT2 = File(f).NT(j(c));
        d(c) = zDistanceToExemplar(Exemplar,NT1,NT2,File(f).Edge(i(a),j(c)));
        if Verbose > 0,
          fprintf('Pair %s %s%5s_%s - %s%5s_%s %s distance %7.4f to exemplar\n', File(f).Filename, NT1.Base,NT1.Number,NT1.Chain,NT2.Base,NT2.Number,NT2.Chain, zEdgeText(File(f).Edge(i(a),j(c))), d(c));
        end
      end
      if Verbose > 1,
        clf
        VP.Sugar = 1;
%        zDisplayNT(File(f),[i(a) j i(a)-1 i(a)+1 min(j)-1 max(j) + 1],VP);
        zDisplayNT(File(f),[i(a) j],VP);
        pause
      end
      for c = 1:length(j),
        if d(c) > min(d),
          if Verbose > 2,
            NT2 = File(f).NT(j(c));
            fprintf('Removing pair %s%5s_%s - %s%5s_%s %s distance %7.4f to exemplar\n', NT1.Base,NT1.Number,NT1.Chain,NT2.Base,NT2.Number,NT2.Chain, zEdgeText(File(f).Edge(i(a),j(c))), d(c));
          end
%          File(f).Edge(i(a),j(c)) = 0;          % remove the worse basepair(s)
%          File(f).Edge(j(c),i(a)) = 0;
        end
      end

      if Verbose > 0,
        fprintf('\n');
      end
    end
  end  
end
