% zListExamplars lists exemplars for each paircode and category

pcodes = [1 5 6 7 9 11 13 14 15 16];
                    % 1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU

% load exemplars from previous session -------------------------------------

 load('PairExemplars','Exemplar');

% loop through paircodes and computer classifications ----------------------

for j = 1:length(pcodes),
 Paircode = pcodes(j);
 for row = 1:length(Exemplar(:,Paircode)),

  E = Exemplar(row,Paircode);

  if ~isempty(E.Filename),
    n1 = E.NT1;
    fprintf('%10s %2s %5s ',  E.Filename, n1.Base, n1.Number); 
    n2 = E.NT2;
    fprintf('%2s %5s ', n2.Base, n2.Number);
    fprintf('%6.2f \n', E.Class);
  end
 end
end
