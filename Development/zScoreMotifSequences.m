
function [scores,Name] = zScoreMotifSequences(codes,Search,Edge,BPh,method,Stacks)

Verbose = 1;

Name{1} = ['Basepairs using method ' num2str(method)];
Name{2} = ['Basepairs using method ' num2str(method) ' and base-phosphates'];
Name{3} = ['Basepairs using ' num2str(method) ', base-phosphates, and stacking'];

[N,M] = size(Edge);
L     = length(Search.Candidates(:,1));         % number of candidates
[a,b] = size(codes);

for k = 1:L,
  f = Search.Candidates(k,N+1);                 % file number
  for n = 1:N,
    CandCodes(k,n) = Search.File(f).NT(Search.Candidates(k,n)).Code;
  end
end



% -------------------------------- score all possible sequences for basepairs

scores = ones(a,1);                    % place to store all scores

if 10 > 1,

for p = 1:N,
  for q = (p+1):N,                              % loop through interactions
    if abs(Edge(p,q)) > 0 && abs(Edge(p,q)) < 13,       % basepair
      IS = zeros(4,4);
      for k = 1:L,                              % loop through instances
        IS = IS + pIsoScore(Edge(p,q),CandCodes(k,p),CandCodes(k,q),method);
      end
      IS = IS / sum(sum(IS));

[p q Edge(p,q)]
IS

      for s = 1:a,                              % loop through all sequences
        scores(s,1) = scores(s,1) * IS(codes(s,p),codes(s,q));
      end
    end
  end
end

end

% -------------------------------- score all possible sequences for BPh

scores(:,2) = scores(:,1);                  % place to store all scores

if 10 > 1,

for p = 1:N,                                      % loop through bases
  for q = 1:N,
    if p~=q && abs(BPh(p,q)) > 0 && abs(BPh(p,q)) < 20,  % BPh consensus here
      P = zeros(1,4);                             % BPh substitution probs
      for k = 1:L,                                % loop through instances
        f = Search.Candidates(k,N+1);
        B = Search.File(f).BasePhosphate;
        x = Search.Candidates(k,p);
        y = Search.Candidates(k,q);
        [D,Q] = zBasePhosphateGeometry(mod(B(x,y),100));

        if ~isempty(Q),
          P = P + Q;
        end
      end

      P = P / sum(P);                             % normalize

[p q]
P

      for s = 1:a,                      % loop through all sequences
        scores(s,2) = scores(s,2) * P(codes(s,p));
      end
    end
  end
end

end

% -------------------------------- score all possible sequences for stacking

scores(:,3) = scores(:,2);                  % place to store all scores

for p = 1:N,                                      % loop through bases
  for q = (p+1):N,                                
    if abs(Edge(p,q)) > 20 && abs(Edge(p,q)) < 24,  % base stack

      switch fix(Edge(p,q)),
      case 21,
        S = Stacks{1};                          % s35
      case {22, -22},
        S = Stacks{2};                          % s33
      case {23, -23},
        S = Stacks{3};                          % s55
      case -21,
        S = Stacks{1};                          % s53
        S.Candidates = S.Candidates(:,[2 1 3]);   % swap nucleotide order
      otherwise,
        disp('What case is this?');
        fix(Edge(p,q))
      end

S.SaveName

      NS = length(S.Candidates(:,1));           % number of stacks
      for c = 1:NS,
        f = S.Candidates(c,3);                  % file number
        SC(c,1) = S.File(f).NT(S.Candidates(c,1)).Code;  % code of first base
        SC(c,2) = S.File(f).NT(S.Candidates(c,2)).Code;  % code of second base
      end

      H = zeros(4,4);

      for k = 1:L,                              % loop through instances
        f = Search.Candidates(k,N+1);           % file number
        Model.NT(1) = Search.File(f).NT(Search.Candidates(k,p));
        Model.NT(2) = Search.File(f).NT(Search.Candidates(k,q));

% zDisplayNT(Model)

        [Discrepancy, Candidates, i] = xRankCandidatesIDI(S.File,Model,S.Candidates,0);

        NewSC = SC(i,:);                        % re-order candidate codes

        MinIDI = 40*ones(4,4);
        
        for c = NS:-1:1,
          MinIDI(NewSC(c,1),NewSC(c,2)) = Discrepancy(c);
        end

        fprintf('We get this min IDI matrix for %s%s %s\n',Model.NT(1).Base,Model.NT(2).Base,zEdgeText(fix(Edge(p,q))));

MinIDI

        fprintf('We turn it into this matrix of scores:\n');

        G = 1 ./ (1+MinIDI.^2);
        G = G / sum(sum(G));

        G

        H = H + G;

      end

fprintf('Overall matrix for this stack is\n');

      H = H / sum(sum(H))

      for s = 1:a,                      % loop through all sequences
        scores(s,3) = scores(s,3) * H(codes(s,p),codes(s,q));
      end
    end
  end
end

end
