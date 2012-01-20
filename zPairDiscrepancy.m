% zPairDiscrepancy(Pair1,Pair2) calculates a measure of the discrepancy
% between Pair1 and Pair2.

function [d] = zPairDiscrepancy(Pair1,Pair2)

t = Pair1.Displ - Pair2.Displ;
r = Pair1.Rot * Pair2.Rot';

[ax,ang] = zAxisAngleRadians(r);              % rotation angle without a flip

d = norm(t) + abs(ang)*(57.29577951308232/20); % convert to degrees, combine

%[norm(t) abs(ang)/20]
