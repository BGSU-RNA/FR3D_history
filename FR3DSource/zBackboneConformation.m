
function [File] = zBackboneConformation(File,Verbose)

if nargin < 2,
  Verbose = 0;
end

File.Backbone = sparse(zeros(length(File.NT)));

if exist('chiropraxis.jar') == 2 ...       % if chiropraxis is here
   && exist(File.PDBFilename) == 2,       % and the PDB file is available

fid = fopen('tempbat.bat','w');

fprintf(fid,'echo off\n');
fprintf(fid,'java -cp chiropraxis.jar chiropraxis.dangle.Dangle rnabb < PDBFiles\\%s.pdb > %s_dangle.txt\n', File.Filename, File.Filename);

fprintf(fid,'suitename.0.3.070628.win.exe < %s_dangle.txt > %s_suitename.txt\n', File.Filename, File.Filename);

%fprintf(fid,'pause\n');

fclose(fid);

!tempbat.bat

fid = fopen([File.Filename '_suitename.txt'],'r');

if fid > 0

  clear T
  L = 1;
  c = 1;

  while L > -1
    L = fgets(fid);
    if L > -1
      if L(1) == ':',
        T{c} = L;
        c = c + 1;
      else
        L = -1;                                   % stop reading file
      end
    end
  end

else

  fprintf('Could not open file %s\n', [File.Filename '_suitename.txt']);

end

fclose(fid);

%123456789012345678901234567890
%:1:A: 388: :  G 33 p 1a 0.630

File.Backbone = sparse(zeros(length(File.NT)));

zBackboneCodes;

for t = 1:length(T),
  a = T{t};
  i = zIndexLookup(File,strrep(a(6:9),' ',''),a(4));
  if isempty(i),
    if Verbose > 0,
      fprintf('%s has no corresponding nucleotide in FR3D\n', a(1:29));
    end
  elseif Verbose > 2,
    fprintf('%s becomes %s%s\n', a(1:29), File.NT(i).Base, File.NT(i).Number);
  end

  bbc = find(ismember(Codes,a(22:23)));       % look up FR3D's internal code

  if isempty(bbc),
    fprintf('Found a new code, %s\n', a(22:23));
  elseif bbc < length(Codes) && ~isempty(i),
    if i > 1,
      if File.Covalent(i-1,i) == 1,             % linked to previous NT
        File.Backbone(i-1,i) = bbc;
      end
    end
  end
end

delete('tempbat.bat');
delete([File.Filename '_dangle.txt']);
delete([File.Filename '_suitename.txt']);

else
  fprintf('Cannot determine backbone conformations for %s\n', File.Filename);
end

return

% determine necessary codes

for t = 1:length(T),
  g{t} = T{t}(22:23);
end
G = unique(g);


% now read the successful output of suitename, convert the codes to numbers, store them, match the codes to numbers in zEdgeText, and implement this in xPairwiseScreen.  Also make it possible to see the backbone conformations for candidates.  Work, but not so bad.


