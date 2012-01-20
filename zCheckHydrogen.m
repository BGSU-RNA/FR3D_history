% This program is generated by zGenerateCheckHydrogen.m
% based on the Excel file H_bonding_Atoms_from_Isostericity_Table.xls
% created by Jesse Stombaugh.
% zCheckHydrogen(NT1,NT2,Class) computes the angles and distances in the
% hydrogen bonds between two nucleotides assuming their interaction is Class
% The program calls the base that should be at the origin N1, the other N2
%
%
function [Hydrogen] = zCheckHydrogen(NT1,NT2,Class)

Paircode = 4*(NT2.Code-1) + NT1.Code;
switch Paircode
  case {2, 3, 4, 8, 10, 12},
    N1 = NT2;
    N2 = NT1;
  otherwise
    N1 = NT1;
    N2 = NT2;
end
Paircode = 4*(N2.Code-1) + N1.Code;

if (Class == 1) & (Paircode == 1),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(11,:),N2.Fit(9,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(11,:));
elseif (Class == 1) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(3,:),N2.Fit(2,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(3,:));
elseif (Class == 1) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(12,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(12,:));
elseif (Class == 1) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(11,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(11,:),N2.Fit(3,:));
elseif (Class == -1) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(10,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(10,:),N2.Fit(3,:));
elseif (Class == 1) & (Paircode == 6),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(4,:),N2.Fit(13,:));
elseif (Class == 1) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(13,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(12,:),N2.Fit(4,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Fit(3,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(16,:),N2.Fit(3,:));
elseif (Class == 1) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(13,:),N2.Fit(6,:));
elseif (Class == -1) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(12,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(12,:));
elseif (Class == 1) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(13,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(12,:),N2.Fit(4,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Fit(3,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(16,:),N2.Fit(3,:));
elseif (Class == 1) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(11,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
elseif (Class == -1) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(11,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(11,:),N2.Fit(3,:));
elseif (Class == -1) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(13,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(11,:));
elseif (Class == -1) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(11,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
elseif (Class == 1) & (Paircode == 16),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(11,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(11,:),N2.Fit(3,:));
elseif (Class == 2) & (Paircode == 1),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(14,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(14,:));
elseif (Class == 2) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(13,:));
elseif (Class == 2) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(11,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(11,:),N2.Fit(6,:));
elseif (Class == -2) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(4,:),N2.Fit(13,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(14,:),N2.Fit(4,:));
elseif (Class == 2) & (Paircode == 6),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(13,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(13,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(3,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(3,:),N2.Fit(13,:));
elseif (Class == 2) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(16,:),N2.Fit(4,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
elseif (Class == 2) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(13,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(13,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(11,:));
elseif (Class == 2) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(16,:),N2.Fit(4,:));
elseif (Class == 2) & (Paircode == 11),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(12,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(12,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(12,:),N2.Fit(6,:));
elseif (Class == 2) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(11,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(12,:),N2.Fit(6,:));
elseif (Class == -2) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(11,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(14,:),N2.Fit(3,:));
elseif (Class == -2) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(4,:),N2.Fit(11,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(13,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(13,:),N2.Fit(3,:));
elseif (Class == -2) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(6,:),N2.Fit(11,:));
elseif (Class == 2) & (Paircode == 16),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(11,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(11,:),N2.Fit(6,:));
elseif (Class == 3) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(10,:),N2.Fit(10,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(10,:),N2.Fit(10,:));
elseif (Class == 3) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(9,:),N2.Fit(8,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(9,:));
elseif (Class == 3) & (Paircode == 6),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N2.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(4,:),N2.Fit(12,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(3,:),N2.Fit(11,:),N2.Fit(8,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(3,:),N2.Fit(11,:));
elseif (Class == 3) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(13,:));
elseif (Class == 3) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(13,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(9,:),N2.Fit(8,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(9,:));
elseif (Class == -3) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(10,:),N2.Fit(12,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(10,:),N2.Fit(12,:));
elseif (Class == 3) & (Paircode == 11),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Fit(10,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(16,:),N2.Fit(10,:));
elseif (Class == -3) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(10,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(10,:),N2.Fit(11,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(7,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
elseif (Class == -3) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(10,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(10,:),N2.Fit(11,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(7,:),N1.Fit(13,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(13,:),N2.Fit(3,:));
elseif (Class == 3) & (Paircode == 16),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(3,:),N2.Fit(9,:),N2.Fit(8,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(3,:),N2.Fit(9,:));
elseif (Class == 4) & (Paircode == 1),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(10,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(10,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(15,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(15,:));
elseif (Class == 4) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(10,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(10,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(10,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(10,:),N2.Fit(6,:));
elseif (Class == -4) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(10,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(10,:),N2.Fit(13,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(15,:),N2.Fit(4,:));
elseif (Class == 4) & (Paircode == 6),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N2.Fit(11,:),N2.Fit(8,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(4,:),N2.Fit(11,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(3,:),N2.Fit(12,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(3,:),N2.Fit(12,:));
elseif (Class == 4) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(10,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(10,:),N2.Fit(13,:));
elseif (Class == 4) & (Paircode == 11),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(10,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(10,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(16,:),N2.Fit(6,:));
elseif (Class == 4) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(9,:),N2.Fit(8,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(9,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(12,:),N2.Fit(6,:));
elseif (Class == -4) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(10,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(10,:),N2.Fit(11,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
elseif (Class == -4) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(7,:),N1.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(13,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(10,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(10,:),N2.Fit(11,:));
elseif (Class == 4) & (Paircode == 16),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(9,:),N2.Fit(8,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(9,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(11,:),N2.Fit(6,:));
elseif (Class == 5) & (Paircode == 1),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(3,:));
elseif (Class == 5) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(3,:));
elseif (Class == 5) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(3,:));
elseif (Class == 5) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(3,:));
elseif (Class == -5) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(13,:));
elseif (Class == 5) & (Paircode == 6),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(13,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(13,:),N2.Fit(3,:));
elseif (Class == 5) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(13,:));
elseif (Class == 5) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(13,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(13,:),N2.Fit(3,:));
elseif (Class == -5) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(3,:),N2.Fit(12,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(3,:),N2.Fit(12,:));
  Hydrogen(3).Angle    = zAngle(N1.Sugar(3,:),N2.Fit(16,:),N2.Fit(11,:));
  Hydrogen(3).Distance = zDistance(N1.Sugar(3,:),N2.Fit(16,:));
elseif (Class == 5) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Sugar(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(16,:),N2.Sugar(3,:));
elseif (Class == 5) & (Paircode == 11),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(15,:),N2.Fit(11,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(15,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Sugar(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(16,:),N2.Sugar(3,:));
elseif (Class == 5) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Sugar(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(16,:),N2.Sugar(3,:));
elseif (Class == -5) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(3,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(3,:),N2.Fit(11,:));
elseif (Class == -5) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(11,:));
elseif (Class == -5) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(11,:),N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(3,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(3,:),N2.Fit(11,:));
elseif (Class == 5) & (Paircode == 16),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(3,:));
elseif (Class == 6) & (Paircode == 1),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(11,:),N2.Fit(9,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(11,:));
elseif (Class == 6) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Sugar(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Sugar(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(14,:),N2.Fit(3,:));
elseif (Class == 6) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Sugar(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Sugar(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(4,:),N2.Fit(15,:),N2.Fit(11,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(4,:),N2.Fit(15,:));
elseif (Class == 6) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Sugar(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Sugar(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(14,:),N2.Fit(3,:));
elseif (Class == -6) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Sugar(3,:),N2.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Sugar(3,:),N2.Fit(12,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(3,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(3,:),N2.Fit(13,:));
elseif (Class == 6) & (Paircode == 6),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(12,:),N2.Sugar(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Sugar(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(13,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(13,:),N2.Fit(3,:));
elseif (Class == 6) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Sugar(3,:),N2.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Sugar(3,:),N2.Fit(12,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(3,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(3,:),N2.Fit(13,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(11,:),N1.Fit(15,:),N2.Fit(4,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(15,:),N2.Fit(4,:));
elseif (Class == 6) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(12,:),N2.Sugar(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Sugar(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(13,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(13,:),N2.Fit(3,:));
elseif (Class == 6) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
elseif (Class == 6) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
elseif (Class == -6) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(11,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(11,:),N2.Fit(3,:));
elseif (Class == -6) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(11,:));
elseif (Class == -6) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(11,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
elseif (Class == 6) & (Paircode == 16),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(11,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(11,:));
elseif (Class == 7) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(6,:));
elseif (Class == 7) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(10,:),N2.Fit(7,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(10,:));
elseif (Class == -7) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(6,:));
elseif (Class == 7) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(10,:),N2.Fit(7,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(10,:));
elseif (Class == 7) & (Paircode == 11),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(13,:),N2.Fit(7,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(13,:));
elseif (Class == 8) & (Paircode == 1),
  Hydrogen(1).Angle    = zAngle(N1.Fit(10,:),N2.Fit(15,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(10,:),N2.Fit(15,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(10,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(15,:),N2.Fit(10,:));
elseif (Class == 8) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(10,:),N2.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(10,:),N2.Fit(12,:));
elseif (Class == 8) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(6,:));
elseif (Class == 8) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(6,:));
elseif (Class == -8) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(10,:),N2.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(10,:),N2.Fit(12,:));
elseif (Class == 8) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(10,:),N2.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(10,:),N2.Fit(12,:));
elseif (Class == 8) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(6,:));
elseif (Class == -8) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(6,:));
elseif (Class == 8) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(10,:),N2.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(10,:),N2.Fit(12,:));
elseif (Class == 8) & (Paircode == 11),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(13,:),N2.Fit(7,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(13,:));
elseif (Class == -8) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(6,:));
elseif (Class == -8) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(6,:));
elseif (Class == 9) & (Paircode == 1),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
elseif (Class == 9) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
elseif (Class == 9) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(10,:),N2.Fit(15,:),N2.Fit(11,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(10,:),N2.Fit(15,:));
elseif (Class == 9) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
elseif (Class == -9) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(12,:));
elseif (Class == 9) & (Paircode == 6),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(7,:),N1.Fit(10,:),N2.Sugar(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(10,:),N2.Sugar(3,:));
elseif (Class == 9) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(12,:));
elseif (Class == 9) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
elseif (Class == -9) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(6,:));
elseif (Class == 9) & (Paircode == 11),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(15,:),N2.Fit(11,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(15,:));
elseif (Class == -9) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(6,:));
elseif (Class == -9) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(9,:),N2.Fit(8,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(9,:));
elseif (Class == -9) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(11,:),N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(6,:));
elseif (Class == 9) & (Paircode == 16),
  Hydrogen(1).Angle    = zAngle(N1.Fit(8,:),N1.Fit(9,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(9,:),N2.Fit(3,:));
elseif (Class == 10) & (Paircode == 1),
  Hydrogen(1).Angle    = zAngle(N1.Fit(10,:),N2.Fit(11,:),N2.Fit(9,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(10,:),N2.Fit(11,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Sugar(3,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(14,:),N2.Sugar(3,:));
elseif (Class == 10) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Sugar(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(14,:),N2.Sugar(3,:));
elseif (Class == 10) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(10,:),N2.Fit(15,:),N2.Fit(11,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(10,:),N2.Fit(15,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Sugar(3,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(14,:),N2.Sugar(3,:));
elseif (Class == 10) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Sugar(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(14,:),N2.Sugar(3,:));
elseif (Class == -10) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(12,:));
  Hydrogen(2).Angle    = zAngle(N1.Sugar(3,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Sugar(3,:),N2.Fit(13,:));
elseif (Class == 10) & (Paircode == 6),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(13,:),N2.Sugar(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(13,:),N2.Sugar(3,:));
elseif (Class == 10) & (Paircode == 14),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(13,:),N2.Sugar(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(13,:),N2.Sugar(3,:));
elseif (Class == 10) & (Paircode == 11),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N2.Fit(15,:),N2.Fit(11,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(6,:),N2.Fit(15,:));
elseif (Class == -10) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(6,:));
elseif (Class == -10) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(11,:),N1.Fit(15,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(6,:));
elseif (Class == 11) & (Paircode == 1),
  Hydrogen(1).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(3,:));
elseif (Class == 11) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(15,:),N2.Fit(11,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(15,:));
elseif (Class == 11) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(11,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
elseif (Class == -11) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(15,:),N2.Fit(11,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(15,:));
elseif (Class == 11) & (Paircode == 11),
  Hydrogen(1).Angle    = zAngle(N1.Fit(11,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
elseif (Class == -11) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(3,:));
elseif (Class == -11) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(11,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
elseif (Class == 12) & (Paircode == 1),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(11,:),N2.Fit(9,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(11,:));
elseif (Class == 12) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(3,:));
elseif (Class == 12) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(3,:),N2.Fit(15,:),N2.Fit(11,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(3,:),N2.Fit(15,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(11,:),N2.Fit(3,:));
elseif (Class == 12) & (Paircode == 13),
  Hydrogen(1).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(3,:));
elseif (Class == -12) & (Paircode == 9),
  Hydrogen(1).Angle    = zAngle(N1.Fit(9,:),N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(11,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(3,:),N2.Fit(15,:),N2.Fit(11,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(3,:),N2.Fit(15,:));
  Hydrogen(3).Angle    = zAngle(N1.Sugar(3,:),N2.Fit(16,:),N2.Fit(11,:));
  Hydrogen(3).Distance = zDistance(N1.Sugar(3,:),N2.Fit(16,:));
elseif (Class == 12) & (Paircode == 7),
  Hydrogen(1).Angle    = zAngle(N1.Fit(11,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Sugar(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(16,:),N2.Sugar(3,:));
elseif (Class == 12) & (Paircode == 11),
  Hydrogen(1).Angle    = zAngle(N1.Sugar(3,:),N2.Fit(16,:),N2.Fit(11,:));
  Hydrogen(1).Distance = zDistance(N1.Sugar(3,:),N2.Fit(16,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(3,:),N2.Fit(15,:),N2.Fit(11,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(3,:),N2.Fit(15,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(11,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(4).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Sugar(3,:));
  Hydrogen(4).Distance = zDistance(N1.Fit(16,:),N2.Sugar(3,:));
elseif (Class == 12) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Sugar(3,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(16,:),N2.Sugar(3,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(11,:),N1.Fit(15,:),N2.Fit(3,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(15,:),N2.Fit(3,:));
elseif (Class == 13) & (Paircode == 1),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N2.Fit(14,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(4,:),N2.Fit(14,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(15,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(15,:));
elseif (Class == 13) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N2.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(4,:),N2.Fit(12,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(4,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(4,:),N2.Fit(13,:));
elseif (Class == -13) & (Paircode == 5),
  Hydrogen(1).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(4,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(6,:),N1.Fit(15,:),N2.Fit(4,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(15,:),N2.Fit(4,:));
  Hydrogen(3).Angle    = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(3,:));
  Hydrogen(3).Distance = zDistance(N1.Fit(14,:),N2.Fit(3,:));
elseif (Class == 13) & (Paircode == 6),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(4,:),N2.Fit(13,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(3,:),N2.Fit(13,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(3,:),N2.Fit(13,:));
elseif (Class == 13) & (Paircode == 11),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(16,:),N2.Fit(6,:));
elseif (Class == 13) & (Paircode == 15),
  Hydrogen(1).Angle    = zAngle(N1.Fit(4,:),N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(1).Distance = zDistance(N1.Fit(12,:),N2.Fit(6,:));
  Hydrogen(2).Angle    = zAngle(N1.Fit(11,:),N1.Fit(16,:),N2.Fit(6,:));
  Hydrogen(2).Distance = zDistance(N1.Fit(16,:),N2.Fit(6,:));
else
Hydrogen = [];
end
