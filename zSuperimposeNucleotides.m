% zSuperimposeNucleotides(File1,NTList1,File2,NTList2) displays NTList1
% from File1 and NTList2 from File2, rotated to align as well as possible
% At the moment, NTList1 and NTList2 must have the same length

% It can be called in several ways, for example,
% zDisplayNT('1s72',{'27','28','29'},ViewParam), and it will load the
% named datafile and plot the nucleotides by nucleotide number.  Or, if
% data files have already been loaded, one can use zDisplayNT(File(1),[32
% 34 35],ViewParam) to plot the nucleotides in File(1) having indices 32,
% 34, and 35.  Defaults for ViewParam are defined in zDisplayNT; see there
% for the fields of ViewParam.
% One can also use ranges of nucleotide numbers, as in
% zDisplayNT('rr0033_23S',{'2548:2555','2557','2559:2566'},VP);

function [void] = zSuperimposeNucleotides(File1,NTList1,File2,NTList2,ViewParam)

% set default values for the display

VP.Sugar     = 1;
VP.az        = 51;
VP.el        = 14;
VP.LineStyle = '-';              % default - thick solid lines
VP.LineThickness = '2';          
VP.AtOrigin  = 1;                % rotate all so first is at origin
VP.Title     = 1;                % title with nucleotide numbers, filename
VP.Grid      = 1;                % add a grid to the graph
VP.FontSize  = 10;               % will use Matlab's default unless overridden
VP.Rotation  = eye(3);
VP.Shift     = zeros(1,3);
VP.LabelBases = 10;

if nargin < 5,
  ViewParam = VP;
end

% replace defaults with defined values

if isfield(ViewParam,'Sugar'),
  VP.Sugar = ViewParam.Sugar;
end

if isfield(ViewParam,'ConnectSugar'),
  VP.ConnectSugar = ViewParam.ConnectSugar;
else
  if VP.Sugar > 0,
    VP.ConnectSugar = 1;
  else
    VP.ConnectSugar = 0;
  end
end

if isfield(ViewParam,'az'),
  VP.az = ViewParam.az;
end

if isfield(ViewParam,'el'),
  VP.el = ViewParam.el;
end

if isfield(ViewParam,'LineStyle'),
  VP.LineStyle = ViewParam.LineStyle;
end

if isfield(ViewParam,'LineThickness'),
  VP.LineThickness = ViewParam.LineThickness;
end

if isfield(ViewParam,'Color'),
  VP.Color = ViewParam.Color;
end

if isfield(ViewParam,'AtOrigin'),
  VP.AtOrigin = ViewParam.AtOrigin;
end

if isfield(ViewParam,'Title'),
  VP.Title = ViewParam.Title;
end

if isfield(ViewParam,'Grid'),
  VP.Grid = ViewParam.Grid;
end

if isfield(ViewParam,'FontSize'),
  h=gca;
  set(h,'FontSize',ViewParam.FontSize);
end

if isfield(ViewParam,'Rotation'),
  VP.Rotation = ViewParam.Rotation;
end

if isfield(ViewParam,'Shift'),
  VP.Shift = ViewParam.Shift;
end

if isfield(ViewParam,'LabelBases'),
  VP.LabelBases = ViewParam.LabelBases;
end

% --------------------------------------- File1 ----------------------------
% if File is a text string (filename), load the file and display

if strcmp(class(File1),'char'),
  File1name = File1;
  File1 = zGetNTData(File1name,0);
end

if nargin == 1,
  NTList1 = 1:File1.NumNT;                  % display them all
end

% if NTList1 is a cell array of numbers, look up the indices

if strcmp(class(NTList1),'char'),
  NTList1 = {NTList1};
end

if strcmp(class(NTList1),'cell'),
  Indices1 = zIndexLookup(File1,NTList1);
else
  Indices1 = NTList1;
end

% --------------------------------------- File22 ----------------------------
% if File2 is a text string (filename), load the file and display

if strcmp(class(File2),'char'),
  File2name = File2;
  File2 = zGetNTData(File2name,0);
end

if nargin == 1,
  NTList2 = 1:File2.NumNT;                  % display them all
end

% if NTList2 is a cell array of numbers, look up the indices

if strcmp(class(NTList2),'char'),
  NTList2 = {NTList2};
end

if strcmp(class(NTList2),'cell'),
  Indices2 = zIndexLookup(File2,NTList2);
else
  Indices2 = NTList2;
end

% ---------------- Check that the lists are the same length

if (length(Indices1) ~= length(Indices2))
  fprintf('The two lists do not have the same length.  One will be shortened.\n');
  if length(Indices1) > length(Indices2)
    Indices1 = Indices1(1:length(Indices2));
  else 
    Indices2 = Indices2(1:length(Indices1));
  end

  % later, find the best match for one inside the other

end

% ---------------- Check that the lists are long enough

if length(Indices1) < 3
  fprintf('You need more than two nucleotides in each set\n');
end

% ---------------- Calculate the discrepancy between the two structures --

[disc,SuperR] = xDiscrepancy(File1,Indices1,File2,Indices2);

fprintf('Geometric discrepancy is %8.4f between the two sets of nucleotides\n', disc);

% ---------------- Plot the nucleotides ------------------------------------

L        = length(Indices1);

Centers1 = cat(1,File1.NT(Indices1).Center);
CC1      = ones(1,L) * Centers1 / L;

Centers2 = cat(1,File2.NT(Indices2).Center);
CC2      = ones(1,L) * Centers2 / L;

R = File1.NT(Indices1(1)).Rot;             % Rotation matrix for first base
R = eye(3);

clf

zPlotNTsRotated(File1,Indices1,VP,R,CC1);

zPlotNTsRotated(File2,Indices2,VP,R*SuperR',CC2);

Title = strcat(File1.NT(Indices1(1)).Base,File1.NT(Indices1(1)).Number);
for j=2:length(Indices1),
  nt = File1.NT(Indices1(j));
  Title = strcat(Title,'-',nt.Base,nt.Number);
end;

Title = strcat(Title,[' ' strrep(File1.Filename,'_','\_')]);

title(Title);
axis equal
view(VP.az,VP.el);

if VP.Grid == 1,
  grid on
else
  grid off
end

rotate3d on

% ---------------- Write PDB file for first set of nucleotides

Filename = File1.Filename;
Filename = [Filename '_' File1.NT(min(Indices1)).Number '_' File1.NT(max(Indices1)).Number];
Filename = [Filename '.pdb'];

fid = fopen(Filename,'w');                     % open for writing

a = 1;                                         % atom number

for i=1:length(Indices1)
  a = zWriteNucleotidePDB(fid,File1.NT(Indices1(i)),a,0,R,CC1);
end

fclose(fid);
fprintf('Wrote %s\n', Filename);

% ---------------- Write PDB file for second set of nucleotides

Filename = File2.Filename;
Filename = [Filename '_' File2.NT(min(Indices2)).Number '_' File2.NT(max(Indices2)).Number];
Filename = [Filename '.pdb'];

fid = fopen(Filename,'w');                     % open for writing

a = 1;                                         % atom number

for i=1:length(Indices2)
  a = zWriteNucleotidePDB(fid,File2.NT(Indices2(i)),a,0,SuperR,CC2);
end

fclose(fid);
fprintf('Wrote %s\n', Filename);
