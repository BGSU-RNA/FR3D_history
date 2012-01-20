
function [scores,Name,stacks] = zScoreMotifSequences(codes,best,Edge,BPh,method,File,Motif,stacks)

Name{1} = ['Basepairs using method ' num2str(method)];
Name{2} = ['Basepairs using method ' num2str(method) ' and base-phosphates'];
Name{3} = ['Basepairs using ' num2str(method) ', base-phosphates, and stacking'];

[N,M] = size(Edge);

[a,b] = size(codes);

% -------------------------------- score all possible sequences for basepairs

scores = ones(a,1);                    % place to store all scores

for p = 1:N,
  for q = (p+1):N,                              % loop through interactions
    if abs(Edge(p,q)) > 0 && abs(Edge(p,q)) < 13,  % basepair
      IS = pIsoScore(Edge(p,q),best(p),best(q),method); % use most common seq


      for s = 1:a,                          % loop through all sequences
        scores(s,1) = scores(s,1) * IS(codes(s,p),codes(s,q));
      end
    end
  end
end

% -------------------------------- score all possible sequences for BPh

scores(:,2) = scores(:,1);                  % place to store all scores

for p = 1:N,                                      % loop through bases
  for q = 1:N,
    if p~=q && abs(BPh(p,q)) > 0 && abs(BPh(p,q)) < 20,  % BPh
      [D,Q] = zBasePhosphateGeometry(BPh(p,q));
%      Q = pBPhSpecificity(BPh(p,q));             % use most common seq

      for s = 1:a,                      % loop through all sequences
        scores(s,2) = scores(s,2) * Q(codes(s,p));
      end
    end
  end
end

% -------------------------------- score all possible sequences for stacking

%stacks{1,1} = [];

if 10 > 1,

fprintf('Scoring based on stacking\n');

scores(:,3) = scores(:,2);                  % place to store all scores

for p = 1:N,                                      % loop through bases
  for q = (p+1):N,                              % loop through interactions
    if abs(Edge(p,q)) > 20 && abs(Edge(p,q)) < 24,  % basepair
      if nargin < 8,
        [D,D1] = xPairSubstitutions(File,Edge(p,q),Motif.NT(p),Motif.NT(q));

        L = length(D(:,1));
        S = ones(4,4);                              % start uniform
        for c = 1:L,
          if D(c,3) < 1,
            S(D(c,1),D(c,2)) = S(D(c,1),D(c,2)) + 1 - D(c,3);
          end
        end

        perc = [0.2300    0.2647    0.3230    0.1823];
        M = perc'*perc;
        S = S ./ M;                                 % use relative frequencies
        S = S / sum(sum(S));                        % normalize

        figure(p*N+q)
        clf
        zDisplayNT(Motif,[p q])
        fprintf('In Figure %d, we get this substitution matrix for %s%s\n',p*N+q,Motif.NT(p).Base,Motif.NT(q).Base);
        S
        drawnow

        stacks{p,q} = S;
      else
        S = stacks{p,q};
      end

      for s = 1:a,                      % loop through all sequences
        scores(s,3) = scores(s,3) * S(codes(s,p),codes(s,q));
      end
    end
  end
end

end
