% zWriteExemplarPDB reads the file of pair exemplars and writes them out to a single PDB file, spaced 20 Angstroms apart in a plane (using Separate = 0) or to individual PDB files, using Separate = 1.

function [void] = zWriteExemplarPDB(Separate)

if nargin < 1,
  Separate = 0;
end

% load exemplars -------------------------------------

load('PairExemplars','Exemplar');

if Separate == 0,

  % loop through pairs and classifications ----------------------

  fid = fopen('Isostericity\PairExemplarPDB.pdb','w');       % open for writing

  a = 1;                                         % atom number
  c = 1;                                         % pair number

  for pc = 1:length(Exemplar(1,:)),
    for row = 1:length(Exemplar(:,pc)),

      E = Exemplar(row,pc);
  
      if ~isempty(E.NT1),

        R = E.NT1.Rot;
        sh = E.NT1.Center;

        a = zWriteNucleotidePDB(fid,E.NT1,a,c,R,sh);
        a = zWriteNucleotidePDB(fid,E.NT2,a,c,R,sh);
        c = c + 1;

      end

    end
  end

  fclose(fid);

  fprintf('Wrote one pair exemplar PDB file.pdb\n');

else                                       % write to separate PDB files

 for c1 = 1:4,
  for c2 = 1:4,
    pc  = 4*(c2-1)+c1;                     % current paircode
    pc2 = 4*(c1-1)+c2;                     % paircode when reversed
    for row = 1:length(Exemplar(:,pc)),

      E = Exemplar(row,pc);
  
      if ~isempty(E.NT1),
        WriteOneExemplar(E,pc);
      end

      if (any(pc == [2 3 4 8 10 12])),     % need to produce symmetric version
        E = Exemplar(row,pc2);             % get the right exemplar
        if ~isempty(E.Class),
          if any(fix(E.Class) == [1 2 7 8])  % but only for these families
            E.Class = -E.Class;            % other side of diagonal
            WriteOneExemplar(E,pc2);
          end
        end
      end

    end


  end
 end

 fprintf('Wrote individual pair exemplar PDB files.pdb\n');

end

% ---------------------------------------------------------------------------

function [void] = WriteOneExemplar(E,pc)

if (E.Class < 0),
  Tem   = E.NT1;                       % reverse order of nucleotides
  E.NT1 = E.NT2;
  E.NT2 = Tem;
  E.Class = -E.Class;
end

if (pc == 7),
  Tem   = E.NT1;                       % reverse order of GC pairs
  E.NT1 = E.NT2;
  E.NT2 = Tem;
end

R  = E.NT1.Rot;
sh = E.NT1.Fit(1,:);

% simple form for Jesse's online database of exemplars

ET = zEdgeText(E.Class,1,E.NT1.Code,E.NT2.Code);
ET = strrep(ET,'WW','Ww');
ET = strrep(ET,'HH','Hh');
ET = strrep(ET,' ','');

FN = [ET '_' E.NT1.Base E.NT2.Base '_Exemplar.pdb'];

CP = norm(E.NT1.Sugar(1,:) - E.NT2.Sugar(1,:));

if ~isempty(strfind(E.Filename, 'Model')),
  E.NT1.Number = '1';
  E.NT2.Number = '2';
end

fprintf('%5s %s%s %s%5s %s%5s %5s %4.1f %5d %4.1f %s\n', ET, E.NT1.Base, E.NT2.Base, E.NT1.Base, E.NT1.Number, E.NT2.Base, E.NT2.Number, E.Filename(1:min(end,5)), CP, E.Count, abs(E.Class), FN);


fid = fopen(['Isostericity\' ET '_' E.NT1.Base E.NT2.Base '_Exemplar.pdb'],'w');

a = 1;                                 % atom number

a = zWriteNucleotidePDB(fid,E.NT1,a,0,R,sh);
a = zWriteNucleotidePDB(fid,E.NT2,a,0,R,sh);

fclose(fid);
