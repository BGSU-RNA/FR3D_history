
File = zAddNTData('1s72',2,File);

VP.AtOrigin = 1;
VP.Sugar = 1;

clf

zDisplayNT(File(1),'A1737 C1738 G1739 U1740 A2039 C2040 G2041 U2042',VP)

M = 1000;
t = 2*pi*(0:M)/M;
hold on

r = 30;

plot3(r*sin(t),r*cos(t),-3*ones(1,M+1));
plot3(r*sin(t),r*cos(t),16*ones(1,M+1));

grid off
title('');

N = 50;

Frames = moviein(N);

figure(1)

for i = 1:N,

  view(i*360/N,0)
  axis vis3d
%  axis off
  drawnow
%  colormap('jet')

  Frames(:,i) = getframe(gca,[0 0 400 400]);
%  axis 

  a = frame2im(Frames(:,1));
  Im(:,:,1,i) = a(:,:,1);
end

% Now save the movie as an mpeg file for use on the Web:
map=colormap    % Uses the previously defined colormap 
mpgwrite(Frames,map,'HelixMovie.mpg')

figure(2)
clf
axis off
movie(Frames,5,5)


