
function [scores,Name] = zScoreMotifSequences(codes,best,Edge,BPh,method)

Name{1} = ['Basepairs using method ' num2str(method)];
Name{2} = ['Basepairs using method ' num2str(method) ' and base-phosphates'];

[N,M] = size(Edge);

[a,b] = size(codes);

% -------------------------------- score all possible sequences for basepairs

scores = ones(a,1);                    % place to store all scores

for p = 1:N,
  for q = 1:N,                              % loop through interactions
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
    
      Q = pBPhSpecificity(BPh(p,q));             % use most common seq

      for s = 1:a,                      % loop through all sequences
        scores(s,2) = scores(s,2) * Q(codes(s,p));
      end
    end
  end
end

