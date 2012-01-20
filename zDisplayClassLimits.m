% zDisplayClassLimits(PairCode) uses colored boxes to display the
% displacement, normal vector, and angle of rotation limits for
% classification of PairCode pairs

function [void] = zDisplayClassLimits(paircode)

ClassLimits = zClassLimits;                   % load class limits

B = ClassLimits(:,:,paircode);                % use limits for this paircode

figure(2)
clf

for row = 1:length(B(:,1)),
   hold on
if (B(row,1) < 15) & (abs(B(row,1)) > 0),
   zBox(B(row,[2 4 6]),B(row,[3 5 7]),B(row,1));
   text(B(row,2)+0.2,B(row,4)+0.3,B(row,7)+0.2,num2str(B(row,1)));
end
end

caxis([-12 12]);

zSponer_Locations                 % read in QM locations of atoms in 4 bases

switch paircode,
  case {1, 5, 9, 13}, X = Sponer_Base(1:15,:,1); m.Code = 1;
  case {6, 14},       X = Sponer_Base(1:13,:,2); m.Code = 2;
  case {7, 11, 15},   X = Sponer_Base(1:16,:,3); m.Code = 3;
  case {16},          X = Sponer_Base(1:12,:,4); m.Code = 4;
end

m.Fit = X;
VP.LineStyle = '-';
VP.Sugar = 0;
zPlotOneNT(m,VP);                                 % plot this nucleotide

%title(['Limits on displacements for paircode ' num2str(paircode)]);
%xlabel('Displacement perpendicular to glycosidic bond');
%ylabel('Displacement parallel to glycosidic bond');
%zlabel('Displacement out of the plane of the base at the origin');

return

figure(3)
clf

for row = 1:length(B(:,1)),
   hold on
   zBox(B(row,[2 8 10]),B(row,[3 9 11]),B(row,1));
   text(B(row,2),B(row,8),B(row,11)+20,num2str(B(row,1)));
end

title(['Limits on normal and rotation angle for paircode ' num2str(paircode)]);
xlabel('Displacement perpendicular to glycosidic bond');
ylabel('Third component of normal vector');
zlabel('Appropriate angle of rotation');

figure(4)
clf

for row = 1:length(B(:,1)),
   hold on
   zBox(B(row,[2 4 10]),B(row,[3 5 11]),B(row,1));
   text(B(row,2),B(row,8),B(row,11)+20,num2str(B(row,1)));
end

title(['Limits on displacements for paircode ' num2str(paircode)]);
xlabel('Displacement perpendicular to glycosidic bond');
ylabel('Displacement parallel to glycosidic bond');
zlabel('Appropriate angle of rotation');

figure(5)
clf

for row = 1:length(B(:,1)),
   hold on
   zBox(B(row,[2 4 8]),B(row,[3 5 9]),B(row,1));
   text(B(row,2),B(row,8),B(row,11)+20,num2str(B(row,1)));
end

title(['Limits on displacements for paircode ' num2str(paircode)]);
xlabel('Displacement perpendicular to glycosidic bond');
ylabel('Displacement parallel to glycosidic bond');
zlabel('Third component of normal vector');

