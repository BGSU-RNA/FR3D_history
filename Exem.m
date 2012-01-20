
% load exemplars -------------------------------------

  load('PairExemplars','Exemplar');

% loop through pairs and classifications ----------------------

fid = fopen('PairExemplarPDB.pdb','w');       % open for writing

a = 1;                                         % atom number
c = 1;                                         % shift

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

fprintf('Wrote PairExemplarPDB.pdb\n');
