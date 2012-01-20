% zPairViewerGUI loads several pdb files, then starts displaying classification
% data using the GUI written by Ali and Jesse.  Still in development.

if ~exist('File'),                  % If files haven't already been loaded
  Filenames = zFileNameList;        % The master list
  for f=1:length(Filenames),
    File(f) = zGetNTData(Filenames{f},0);
  end
end

Param.Decimal   = 0;        % 1 - use 1.0 only; 0 - round to 1
Param.Group     = 3;        % hand, computer, both, etc.
Param.Expert    = 0;
Param.Inbox     = [-5 5 -5 5 -5 5 -1.1 1.1 -95 275];  % almost everything

ViewParam.Normal    = 1;
ViewParam.ColorAxis = [-12 30];
ViewParam.SortKeys  = [];
ViewParam.Sugar     = 0;
ViewParam.Nearby    = 0;
ViewParam.az        = 51;                  % azimuth
ViewParam.el        = 30;                  % angle

ViewParam.FigNum = 1;

mGUI(File,Param,ViewParam);

