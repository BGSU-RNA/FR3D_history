% zGetPDBInfo(File) reads a data file about PDB files and extracts
% information about the current molecule stored in File

function [File] = zGetPDBInfo(File)

load('PDBInfo.mat','n','t','-mat');

PDBNames = lower(t(:,1));              % convert PDB filenames to lowercase

i = strmatch(lower(File.Filename),PDBNames);

if ~isempty(i),
  File.Info.Resolution  = n(i(1),1);
  File.Info.Descriptor  = t{i(1),2};
  File.Info.ReleaseDate = t{i(1),3};
  File.Info.Author      = t{i(1),4};
  File.Info.Keywords    = t{i(1),5};
  File.Info.Source      = t{i(1),7};
  if n(i(1),2) ~= length(File.NT),
    n(i(1),2) = length(File.NT);
    BaseSequence = cat(2,File.NT.Base);
    t{i(1),9} = BaseSequence;
  end

  save('PDBInfo.mat','n','t'); % Matlab version 7
else
  File.Info.Resolution  = [];
  File.Info.Descriptor  = [];
  File.Info.ReleaseDate = [];
  File.Info.Author      = [];
  File.Info.Keywords    = [];
  File.Info.Source      = [];
end


return

[n,t,r] = xlsread('PDB_File_Information.xls');
t = t(2:end,:);                  % remove the header row
t{1,9} = '';                     % make space to save base sequence
n(1,2) = 0;                      % make space to save number of nucleotides
save('PDBInfo.mat','n','t'); % Matlab version 7

%save PDBInfo.mat n t -V6      % for compatibility with earlier versions
