% zPlotSomeNTAtOrigin puts the first base at the origin and the others in
% their respective locations

function [R,S] = zPlotSomeNTAtOrigin(File,C,ViewParam);

VP = ViewParam;
VP.AtOrigin = 1;

R = File.NT(C(1)).Rot;
S = File.NT(C(1)).Fit(1,:);

zDisplayNT(File,C,VP);

return








% old code

R = File.NT(C(1)).Rot;                     % Rotation matrix for first base
S = File.NT(C(1)).Fit(1,:);               % Location of glycosidic atom

for j=1:length(C),                         % Loop through all nucleotides
  k = C(j);                                % index of current nucleotide
  NT.Code  = File.NT(k).Code;
  L = length(File.NT(k).Fit(:,1));         % Number of base atoms
  NT.Fit   = (File.NT(k).Fit   - ones(L,1) *S) * R; % rotated into position
  s = length(File.NT(k).Sugar(:,1));                % this sometimes varies
  NT.Sugar = (File.NT(k).Sugar - ones(s,1)*S) * R; % rotate sugar too
  zPlotOneNT(NT,ViewParam);
end

if ViewParam.Sugar > 0,
  [D,i] = sort(C);
  for j=1:(length(C)-1),
    if D(j+1)-D(j) == 1,
      A = [File.NT(D(j)).Sugar(5,:); File.NT(D(j+1)).Sugar(10,:)];
      AA = (A - ones(2,1)*S) * R;
      plot3(AA(:,1), AA(:,2), AA(:,3),'k','LineWidth',2,'LineStyle',ViewParam.LineStyle);
    end
  end
end
