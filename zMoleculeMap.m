% zMoleculeMap(File) puts at dot at the center of each nucleotide
% in File, colored by the number of interactions the nucleotide has (blue
% means zero, red is the most).  The C1' atom of each nucleotide is
% connected to the C1' atom of the next nucleotide by a black line.
% Interactions are shown between centers of nucleotides.  WC-WC
% interactions (category 1) are shown by red lines, non-canonical planar
% interactions by green, typical sequential stacking (categories 15, 17) by
% dark blue, and typical cross-strand stacking (categories 16, 18) by light
% blue.
% File can either be, for example, 'rr0033_5S' or File(3).

function [void] = zMoleculeMap(File)

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

%-----------------------------------------------------------------------
figure(1)
clf

set(gcf,'Renderer','OpenGL');

for k = 1:length(File.NT),                          % Loop through nucleotides
  ii = find((File.Inter(k,:)<30).*(abs(File.Inter(k,:))>0));
  ni(k) = length(ii);                           % number of interactions
  e = File.NT(k).Center;

  scatter3(e(1),e(2),e(3),12,ni(k),'filled')
  hold on
  for j = 1:length(ii),
    if ii(j) < k,
      f = File.NT(ii(j)).Center;
      inter = fix(File.Inter(k,ii(j)));
      if inter == 1,
        plot3([e(1) f(1)],[e(2) f(2)],[e(3) f(3)],'r');  % WC-WC
      elseif abs(inter) < 15,
        plot3([e(1) f(1)],[e(2) f(2)],[e(3) f(3)],'g');  % planar
      elseif (inter == 15) | (inter == 17),
        plot3([e(1) f(1)],[e(2) f(2)],[e(3) f(3)],'b');  % stacking
      elseif (inter == 16) | (inter == 18),
        plot3([e(1) f(1)],[e(2) f(2)],[e(3) f(3)],'c');  % stacking
      end
    end
  end

  g = File.NT(k).Loc(1,:);                         
  h = File.NT(max(1,k-1)).Loc(1,:);
  plot3([g(1) h(1)],[g(2) h(2)],[g(3) h(3)],'k')  % connect glycosidic atoms


end

%caxis([1 4]);
title('Colored by number of interactions');

return


return

%-----------------------------------------------------------------------
figure(ViewParam.FigNum)
clf

set(gcf,'Renderer','OpenGL');

for k = 1:length(File.NT),                          % Loop through nucleotides
  c = File.NT(k).Code;
  e = File.NT(k).Center;

  scatter3(e(1),e(2),e(3),12,c,'filled')
  hold on
end

caxis([1 4]);
title('Colored by nucleotide type');

return
