% zBox(X,Y,color) plots a rectangular box with lower left corner X and
% upper right corner Y, using color

function [void] = zBox(X,Y,color)

set(gcf,'Renderer','OpenGL');

a = X(1);
b = X(2);
c = X(3);

x = Y(1);
y = Y(2);
z = Y(3);

%fill3([a,x,x,a],[b,b,y,y],[c,c,c,c],color);
%fill3([a,x,x,a],[b,b,y,y],[z,z,z,z],color);
fill3([a,a,a,a],[b,b,y,y],[c,z,z,c],color);
fill3([x,x,x,x],[b,b,y,y],[c,z,z,c],color);
fill3([a,x,x,a],[b,b,b,b],[c,c,z,z],color);
fill3([a,x,x,a],[y,y,y,y],[c,c,z,z],color);
