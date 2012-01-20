% zExtractAtomsPDB(Filename) reads Filename.pdb and writes out lines 
% containing ATOM to Filename_Atoms.pdb

function [void] = zExtractAtomsPDB(Filename,Outputname)

fid = fopen([Filename '.pdb'],'r');

if fid > 0

  out = fopen([Outputname '.pdb'],'w');

  L = 1;

  while L > -1
    L = fgets(fid);
    if L > -1
      if strcmp(L(1:min(4,length(L))),'ATOM'),
        fprintf(out,'%s',L);
      end
    end
  end

  fclose(fid);
  fclose(out);

  fprintf('Read  %s.pdb\n', Filename)

else

  fprintf('Could not open file %s.pdb\n', Filename);

end
