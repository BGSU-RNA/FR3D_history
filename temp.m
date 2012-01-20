for v = 1:4,
  figure(v)
  clf
  switch v,
    case 1,     c = PHA(:,4);
                scatter3(PHA(:,1), PHA(:,2), PHA(:,3),4,c,'filled')
    case 2,     c = PHC(:,4);
                scatter3(PHC(:,1), PHC(:,2), PHC(:,3),4,c,'filled')
    case 3,     
  zPlotStandardBase(v,1,1);                % plot base at the origin
  hold on
  rotate3d on
  axis equal
  view(2)
[y,i] = sort(PHG(:,5));
for k = 1:length(i),
  c = PHG(:,4);
  scatter3(PHG(i(k),1), PHG(i(k),2), PHG(i(k),3),4,c(i(k)),'filled')
  drawnow
  rotate3d on
  pause
end
    case 4,     c = PHU(:,4);
                scatter3(PHU(:,1), PHU(:,2), PHU(:,3),4,c,'filled')
  end

  L = {'A','C','G','U'};

  zPlotStandardBase(v,1,1);                % plot base at the origin
  rotate3d on
  axis equal
  view(2)
end
