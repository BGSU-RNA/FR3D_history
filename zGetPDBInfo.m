% zGetPDBInfo(File) reads a data file about PDB files and extracts
% information about the current molecule stored in File

function [File] = zGetPDBInfo(File)

[n,t,r] = xlsread('PDB_File_Headers.xls');

PDBNames = lower(t(:,2));              % convert PDB filenames to lowercase

i = strmatch(lower(File.Filename),PDBNames);

if ~isempty(i),
  File.Info.Resolution = n(i(1)-1,1);
  File.Info.Type       = t{i(1),4};
  File.Info.RNA        = t{i(1),5};
  File.Info.Species    = t{i(1),6};
  File.Info.LigandsAndComments = t{i(1),7};
else
  File.Info.Resolution = [];
  File.Info.Type       = '';
  File.Info.RNA        = '';
  File.Info.Species    = '';
  File.Info.LigandsAndComments = '';
end



