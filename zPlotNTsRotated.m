
function [void] = zPlotNTsRotated(File,Indices,VP,R,S)

for k=1:length(Indices),                 % Loop through all nucleotides
  zPlotOneNTRotated(File.NT(Indices(k)),VP,R,S);
end

% connect the sugars, if sugars are being plotted

if (VP.Sugar > 0) & (VP.ConnectSugar > 0),
  [D,i] = sort(Indices);
  for j=1:(length(Indices)-1),
    if D(j+1)-D(j) == 1,
      A = [File.NT(D(j)).Sugar(5,:); File.NT(D(j+1)).Sugar(10,:)];
      AA = (A - ones(2,1)*S) * R;
      plot3(AA(:,1), AA(:,2), AA(:,3),'k','LineWidth',2,'LineStyle','-');
    end
  end
end
