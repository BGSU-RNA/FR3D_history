% zSuperimposeNucleotides(File1,NTList1,File2,NTList2,ViewParam,L) 
% superimposes the first L nucleotides in NTList1
% from File1 and NTList2 from File2, rotated to align as well as possible
% It displays all nucleotides from both lists.

% This is not totally relevant:
% It can be called in several ways, for example,
% zDisplayNT('1s72',{'27','28','29'},ViewParam), and it will load the
% named datafile and plot the nucleotides by nucleotide number.  Or, if
% data files have already been loaded, one can use zDisplayNT(File(1),[32
% 34 35],ViewParam) to plot the nucleotides in File(1) having indices 32,
% 34, and 35.  Defaults for ViewParam are defined in zDisplayNT; see there
% for the fields of ViewParam.
% One can also use ranges of nucleotide numbers, as in
% zDisplayNT('rr0033_23S',{'2548:2555','2557','2559:2566'},VP);

function [disc,Shift,SuperR] = zSuperimposeNucleotides(File1,NTList1,File2,NTList2,ViewParam,L)

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
VP.Plot      = 1;
VP.Write     = 0;

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

if isfield(ViewParam,'Plot'),
  VP.Plot = ViewParam.Plot;
end

if isfield(ViewParam,'Write'),
  VP.Write = ViewParam.Write;
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

if nargin < 6,
  L = length(Indices1);
end

L = min([L length(Indices1) length(Indices2)]);

% ---------------- Check that the lists are long enough

if L < 3
  fprintf('You need more than two nucleotides in each set\n');
end

if isfield(ViewParam,'Phosphate'),
  LW = ViewParam.Phosphate;
else
  LW = ones(1,L);
end

% ---------------- Calculate the discrepancy between the two structures --

[disc,SuperR,CC1,CC2] = xDiscrepancy(File1,Indices1(1:L),File2,Indices2(1:L),LW);

Shift = CC2-CC1;

if isfield(ViewParam,'Rotation'),
  R = VP.Rotation;
else
  R = eye(3);
end

% ---------------- Plot the nucleotides ------------------------------------

if VP.Plot > 0

  clf

  VP.Color = [1 0 0];
  zPlotNTsRotated(File1,Indices1(1:L),VP,R,CC1);

  VP.Color = [0.5 0 0];
  zPlotNTsRotated(File1,Indices1((L+1):end),VP,R,CC1);

  VP.Color = [0 1 0];
  zPlotNTsRotated(File2,Indices2(1:L),VP,SuperR*R,CC2);

  VP.Color = [0 0.5 0];
  zPlotNTsRotated(File2,Indices2((L+1):end),VP,SuperR*R,CC2);

  Title = [File1.Filename '(red) versus ' File2.Filename '(green)'];
  title(Title);
  axis equal
  view(VP.az,VP.el);

  if VP.Grid == 1,
    grid on
  else
    grid off
  end

  rotate3d on
end

% ---------------- Write PDB file for first set of nucleotides

if VP.Write > 0,

Filename = File1.Filename;
Filename = [Filename '_' File1.NT(min(Indices1(1:L))).Number '_' File1.NT(max(Indices1(1:L))).Number];
Filename = [Filename '.pdb'];

fid = fopen(Filename,'w');                     % open for writing

a = 1;                                         % atom number

for i=1:length(Indices1(1:L))
  a = zWriteNucleotidePDB(fid,File1.NT(Indices1(i)),a,0,R,CC1);
end

fclose(fid);
fprintf('Wrote %s\n', Filename);

% ---------------- Write PDB file for second set of nucleotides

Filename = File2.Filename;
Filename = [Filename '_' File2.NT(min(Indices2(1:L))).Number '_' File2.NT(max(Indices2(1:L))).Number];
Filename = [Filename '.pdb'];

fid = fopen(Filename,'w');                     % open for writing

a = 1;                                         % atom number

for i=1:length(Indices2)
  a = zWriteNucleotidePDB(fid,File2.NT(Indices2(i)),a,0,SuperR*R,CC2);
  a = a + 1;
end

fclose(fid);
fprintf('Wrote %s\n', Filename);

end