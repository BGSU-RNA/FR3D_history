% pConsensusPairSubstitution(a,b,f,File,F,L,Search) looks at the letter pairs corresponding to nucleotides a and b of the query motif in Search, 

function [Score] = pConsensusPairSubstitution(a,b,f,File,F,L,Search,Param)

method          = 4;             % method for assigning pair subst probs

if length(Param) > 1,
  method  = Param(2);
end

Verbose = Param(1);

load PairExemplars

Score = zeros(4,4);              % ready to sum IsoScores for this pair 

for c = 1:L,                            % loop through candidates

  i   = Search.Candidates(c,a);         % index 1 of pair in candidate
  j   = Search.Candidates(c,b);         % index 2 of pair in candidate

if 0 > 1,
Search
c
f
File
File(f(c))
i
j
end

  NT1 = File(f(c)).NT(i);               % retrieve the first nucleotide
  NT2 = File(f(c)).NT(j);               % retrieve the second nucleotide

  if Verbose > 0,
    fprintf('File %4s has %s%4s and %s%4s making %4s; consensus is %4s\n', File(f(c)).Filename, NT1.Base, NT1.Number, NT2.Base, NT2.Number, zEdgeText(File(f(c)).Edge(i,j)), zEdgeText(F.Edge(a,b)));
  end

  newScore = pIsoScore(F.Edge(a,b),NT1.Base,NT2.Base,method,ExemplarIDI,ExemplarFreq);
                                             % use consensus edge
  Score = Score + newScore;

%newScore(1:7)
end

Score = Score / L;                          % average *** need better idea!

%Score