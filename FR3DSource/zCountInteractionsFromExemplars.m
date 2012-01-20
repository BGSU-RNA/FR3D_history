% zCountInteractionsFromExemplars summarizes information on pair interaction frequency from the exemplar file

CountFilename = 'CountsFromExemplars.xls'

MaxCat = 25;
Count = zeros(4,4,MaxCat);

load('PairExemplars','Exemplar');

for c1 = 1:4,                              % loop through code of nucleotide 1
 for c2 = 1:4,                             % loop through code of nucleotide 2
  for r = 1:length(Exemplar(:,1)),         % loop through rows of Exemplar
   pc  = 4*(c2-1)+c1;                      % current paircode

   E = Exemplar(r,pc);                     % current exemplar

   if ~isempty(E.Filename),
     ca    = fix(abs(E.Class));            % current category

     if (E.Class < 0),
       Count(c2,c1,ca) = Count(c2,c1,ca) + E.Count;
     else
       Count(c1,c2,ca) = Count(c1,c2,ca) + E.Count;
     end

    end
  end
 end
end

c1 = 3;
c2 = 2;

for ca = 1:MaxCat,
  t = Count(c1,c2,ca);
  Count(c1,c2,ca) = Count(c2,c1,ca);       % reverse GC and CG counts
  Count(c2,c1,ca) = t;
end

zDisplayPairCounts(Count,CountFilename);

