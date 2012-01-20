% zPlotExemplar(E,ViewParam)

function [void] = zPlotExemplar(E,ViewParam)

C = [E.Base1Index E.Base2Index];
R = E.NT1.Rot;
S = E.NT1.Fit(1,:);
zPlotOneNTRotated(E.NT1,ViewParam,R,S);
zPlotoneNTRotated(E.NT2,ViewParam,R,S);

grid on
%axis equal

