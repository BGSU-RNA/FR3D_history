
load PairExemplars

Filenames = zReadPDBList('NonRedundant_2008_02_21_list');

for f = 1:length(Filenames),
  Filenames{f} = upper(Filenames{f});
end

[s,t] = size(Exemplar);

c = 1;

for i = 1:s,
  for j = 1:t,
    fn = upper(Exemplar(i,j).Filename);
    if ~isempty(fn),
      p  = find(ismember(Filenames,fn));
      if isempty(p),
        fprintf('Exemplar (%d,%d) was found in %s\n', i,j,fn);
      end
    end
  end
end
