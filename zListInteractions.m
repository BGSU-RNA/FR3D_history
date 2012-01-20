% zListInteractions(File,i,j) 

function [void] = zListInteractions(File,i,j)

[y,L] = sort(abs(File.Inter(i,j)));

j = j(L);

fprintf('%10s: %1s%4s with ', File.Filename, File.NT(i).Base, File.NT(i).Number);
for k=1:length(j),
  fprintf('%1s%4s %3.0f | ', File.NT(j(k)).Base, File.NT(j(k)).Number, ...
  File.Inter(i,j(k)));
end
fprintf('\n');