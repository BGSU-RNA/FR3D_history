% pConsensusInteractions(Search) determines the consensus interactions in the candidates in Search

function [Edge,BPh] = pConsensusInteractions(Search)

[L,N] = size(Search.Candidates);        % L = num instances; N = num NT
N = N - 1;                              % number of nucleotides

f = Search.Candidates(:,N+1);           % file numbers of motifs

% ----------------------------------------------- find consensus basepairs

for a = 1:N,                                    % first NT of possible pair
  for b = 1:N,                                  % second NT of possible pair
    e = [];                                     % record observed edges
    for c = 1:L,                                % run through candidates
      i = Search.Candidates(c,a);               % index of first nucleotide
      j = Search.Candidates(c,b);               % index of second nucleotide
      e = [e Search.File(f(c)).Edge(i,j)];      % append observed interaction
    end

    e = fix(e);                                 % round subcategories

    for d = 1:length(e),
      if any(e(d) == [-1 -2 -7 -8]),            % don't distinguish sign here
        e(d) = -e(d);
      end
    end

    if max(abs(e)) > 0 && min(abs(e)) < 30,     % there was some bp interaction
%  e
      Edge(a,b) = mode(e);                    % use the most common one
  % a more sophisticated method for determining the consensus is needed
  % count how many times each one occurs, giving maybe 0.5 for near
    end
  end
end

% ----------------------------------------------- consensus BPh interactions

for a = 1:N,                                    % first NT of possible BPh
  for b = 1:N,                                  % second NT of possible BPh
    e = [];                                     % record observed edges
    for c = 1:L,                                % run through candidates
      i = Search.Candidates(c,a);               % index of first nucleotide
      j = Search.Candidates(c,b);               % index of second nucleotide
      e = [e Search.File(f(c)).BasePhosphate(i,j)];      % append observed interaction
    end

    e = fix(e);                                 % round subcategories

    if max(abs(e)) > 0 && min(abs(e)) < 30,     % there was some bp interaction
%  e
      BPh(a,b) = mode(e);                    % use the most common one
  % a more sophisticated method for determining the consensus is needed
  % count how many times each one occurs, giving maybe 0.5 for near
    end
  end
end
