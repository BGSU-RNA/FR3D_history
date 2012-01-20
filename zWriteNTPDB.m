% zWriteNTPDB(File,NTList,Filename) writes NTList from File to 
% (optional) Filename in PDB format, including hydrogen atoms
% NTList can be a character, as in '23' or a cell array
% One can specify a list using ranges of nucleotide numbers, as in
% zWriteNTPDB('1s72',{'2548:2555','2557','2559:2566'});
% One can specify a chain in these ways:
% 2548(A) or 2548:2555(A) or 2548(A):2555(A) or 2548(A):2555

function [void] = zWriteNTPDB(File,NTList,Filename)

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

if nargin == 1,
  NTList = 1:File.NumNT;                  % display them all
end

if strcmp(class(NTList),'char'),
  NTList = {NTList};
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(NTList),'cell'),
  Indices = zIndexLookup(File,NTList);
else
  Indices = NTList;
end

if nargin < 3,
  Filename = File.Filename;
  if strcmp(class(NTList),'cell')
    for i = 1:length(NTList),
      Filename = [Filename '_' NTList{i}];
    end
  else
    Filename = [Filename '_' File.NT(min(Indices)).Number '_' File.NT(max(Indices)).Number];
  end
  Filename = [Filename '.pdb'];
end

% ---------------- write the file

fid = fopen(Filename,'w');                     % open for writing

a = 1;                                         % atom number

for i=1:length(Indices)
  R  = eye(3);
  sh = [0 0 0];

  a = zWriteNucleotidePDB(fid,File.NT(Indices(i)),a,0,R,sh);
end

fclose(fid);

fprintf('Wrote %s\n', Filename);
  